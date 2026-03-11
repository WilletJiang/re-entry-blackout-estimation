use crate::data::{AeroCoeff, AeroTable};
use crate::types::DynamicsParams;
use nalgebra::{UnitQuaternion, Vector3};

pub fn isa_atmosphere(h_m: f64) -> (f64, f64) {
    let h_m = h_m.max(0.0);
    let (t, rho) = if h_m < 11_000.0 {
        let t = 288.15 - 0.0065 * h_m;
        let rho = 1.225 * (t / 288.15).powf(4.2561);
        (t, rho)
    } else if h_m < 20_000.0 {
        let t = 216.65;
        let exponent = -(h_m - 11_000.0) / 6341.62;
        let rho = 0.36391 * exponent.exp();
        (t, rho)
    } else if h_m < 32_000.0 {
        let t = 196.65 + 0.001 * (h_m - 20_000.0);
        let rho = 0.08803 * (t / 216.65).powf(-34.1632);
        (t, rho)
    } else if h_m < 47_000.0 {
        let t = 228.65 + 0.0028 * (h_m - 32_000.0);
        let rho = 0.01322 * (t / 228.65).powf(-12.2011);
        (t, rho)
    } else {
        let t = 270.65;
        let exponent = -(h_m - 47_000.0) / 6341.62;
        let rho = 0.00143 * exponent.exp();
        (t, rho)
    };
    let a = (1.4 * 287.05287 * t).sqrt();
    (rho.max(1e-6), a)
}

pub struct AeroResult {
    pub force_b: Vector3<f64>,
    pub torque_b: Vector3<f64>,
    pub coeff: AeroCoeff,
    pub q_dyn: f64,
}

pub struct NominalBodyDynamics {
    pub aero: AeroResult,
    pub command_force_b: Vector3<f64>,
    pub command_torque_b: Vector3<f64>,
    pub total_force_b: Vector3<f64>,
    pub total_torque_b: Vector3<f64>,
    pub specific_force_b: Vector3<f64>,
}

pub fn aero_force_torque(
    table: &AeroTable,
    dyn_params: &DynamicsParams,
    r_i: Vector3<f64>,
    v_i: Vector3<f64>,
    q_ib: &UnitQuaternion<f64>,
    omega_b: Vector3<f64>,
) -> AeroResult {
    let speed = v_i.norm().max(1e-3);
    let alt = r_i.norm() - dyn_params.radius_earth;
    let (rho, a) = isa_atmosphere(alt);
    let mut mach = speed / a;

    if !table.mach.is_empty() {
        let mach_min = table.mach[0];
        let mach_max = *table.mach.last().unwrap();
        if mach < mach_min {
            mach = mach_min;
        } else if mach > mach_max {
            mach = mach_max;
        }
    }

    let c_bi = q_ib.to_rotation_matrix();
    let v_b = c_bi.transpose() * v_i;
    let v_b_norm = v_b.norm().max(1.0e-6);
    let mut alpha = v_b.z.atan2(v_b.x).to_degrees();
    let mut beta = (v_b.y / v_b_norm).asin().to_degrees();

    if !table.alpha.is_empty() {
        let alpha_min = table.alpha[0];
        let alpha_max = *table.alpha.last().unwrap();
        if alpha < alpha_min {
            alpha = alpha_min;
        } else if alpha > alpha_max {
            alpha = alpha_max;
        }
    }
    if !table.beta.is_empty() {
        let beta_min = table.beta[0];
        let beta_max = *table.beta.last().unwrap();
        if beta < beta_min {
            beta = beta_min;
        } else if beta > beta_max {
            beta = beta_max;
        }
    }

    let coeff = table.interp(mach, alpha, beta);
    let q_dyn = 0.5 * rho * speed * speed;
    let force_coeff = Vector3::new(coeff.cx, coeff.cy, coeff.cz);
    let mut moment_coeff = Vector3::new(coeff.cl, coeff.cm, coeff.cn);
    let inv_two_speed = 0.5 / speed;
    let p_hat = omega_b.x * dyn_params.b_ref * inv_two_speed;
    let q_hat = omega_b.y * dyn_params.c_ref * inv_two_speed;
    let r_hat = omega_b.z * dyn_params.b_ref * inv_two_speed;
    moment_coeff.x += coeff.cl_p * p_hat;
    moment_coeff.y += coeff.cm_q * q_hat;
    moment_coeff.z += coeff.cn_r * r_hat;
    let force_b = q_dyn * dyn_params.s_ref * force_coeff;
    let torque_b = q_dyn
        * dyn_params.s_ref
        * Vector3::new(
            dyn_params.b_ref * moment_coeff.x,
            dyn_params.c_ref * moment_coeff.y,
            dyn_params.b_ref * moment_coeff.z,
        );
    AeroResult {
        force_b,
        torque_b,
        coeff,
        q_dyn,
    }
}

pub fn nominal_body_dynamics(
    table: &AeroTable,
    dyn_params: &DynamicsParams,
    r_i: Vector3<f64>,
    v_i: Vector3<f64>,
    q_ib: &UnitQuaternion<f64>,
    omega_b: Vector3<f64>,
    force_cmd_b: Vector3<f64>,
    torque_cmd_b: Vector3<f64>,
) -> NominalBodyDynamics {
    nominal_body_dynamics_with_residual(
        table,
        dyn_params,
        r_i,
        v_i,
        q_ib,
        omega_b,
        force_cmd_b,
        torque_cmd_b,
        Vector3::zeros(),
        Vector3::zeros(),
    )
}

pub fn nominal_body_dynamics_with_residual(
    table: &AeroTable,
    dyn_params: &DynamicsParams,
    r_i: Vector3<f64>,
    v_i: Vector3<f64>,
    q_ib: &UnitQuaternion<f64>,
    omega_b: Vector3<f64>,
    force_cmd_b: Vector3<f64>,
    torque_cmd_b: Vector3<f64>,
    aero_res_force_b: Vector3<f64>,
    aero_res_torque_b: Vector3<f64>,
) -> NominalBodyDynamics {
    let aero = aero_force_torque(table, dyn_params, r_i, v_i, q_ib, omega_b);
    let total_force_b = aero.force_b + force_cmd_b + aero_res_force_b;
    let total_torque_b = aero.torque_b + torque_cmd_b + aero_res_torque_b;
    let specific_force_b = total_force_b / dyn_params.mass;

    NominalBodyDynamics {
        aero,
        command_force_b: force_cmd_b,
        command_torque_b: torque_cmd_b,
        total_force_b,
        total_torque_b,
        specific_force_b,
    }
}

#[cfg(test)]
mod tests {
    use super::{nominal_body_dynamics, nominal_body_dynamics_with_residual, AeroCoeff, AeroTable};
    use crate::types::DynamicsParams;
    use nalgebra::{UnitQuaternion, Vector3};

    fn single_coeff_table() -> AeroTable {
        AeroTable {
            mach: vec![1.0],
            alpha: vec![0.0],
            beta: vec![0.0],
            coeffs: vec![AeroCoeff {
                cx: 1.0,
                cy: 0.0,
                cz: 0.0,
                cl: 0.0,
                cm: 0.0,
                cn: 0.0,
                cl_p: 0.0,
                cm_q: 0.0,
                cn_r: 0.0,
            }],
            n_a: 1,
            n_b: 1,
        }
    }

    #[test]
    fn nominal_dynamics_combines_aero_and_rcs_terms() {
        let table = single_coeff_table();
        let dyn_params = DynamicsParams::default();
        let command_force_b = Vector3::new(10.0, -5.0, 3.0);
        let command_torque_b = Vector3::new(2.0, -1.0, 0.5);
        let nominal = nominal_body_dynamics(
            &table,
            &dyn_params,
            Vector3::new(dyn_params.radius_earth + 20_000.0, 0.0, 0.0),
            Vector3::new(1_500.0, 0.0, 0.0),
            &UnitQuaternion::identity(),
            Vector3::zeros(),
            command_force_b,
            command_torque_b,
        );

        assert!((nominal.total_force_b - (nominal.aero.force_b + command_force_b)).norm() < 1.0e-9);
        assert!(
            (nominal.total_torque_b - (nominal.aero.torque_b + command_torque_b)).norm() < 1.0e-9
        );
    }

    #[test]
    fn residual_terms_flow_through_shared_nominal_helper() {
        let table = single_coeff_table();
        let dyn_params = DynamicsParams::default();
        let aero_res_force_b = Vector3::new(4.0, -2.0, 1.0);
        let aero_res_torque_b = Vector3::new(0.5, -0.25, 0.125);
        let nominal = nominal_body_dynamics_with_residual(
            &table,
            &dyn_params,
            Vector3::new(dyn_params.radius_earth + 20_000.0, 0.0, 0.0),
            Vector3::new(1_500.0, 0.0, 0.0),
            &UnitQuaternion::identity(),
            Vector3::zeros(),
            Vector3::zeros(),
            Vector3::zeros(),
            aero_res_force_b,
            aero_res_torque_b,
        );

        assert!(
            (nominal.total_force_b - (nominal.aero.force_b + aero_res_force_b)).norm() < 1.0e-9
        );
        assert!(
            (nominal.total_torque_b - (nominal.aero.torque_b + aero_res_torque_b)).norm() < 1.0e-9
        );
    }
}
