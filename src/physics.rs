use crate::data::{AeroCoeff, AeroTable};
use crate::types::DynamicsParams;
use nalgebra::{UnitQuaternion, Vector3};

pub fn isa_atmosphere(h_m: f64) -> (f64, f64) {
    // Clamp to sea level: negative altitudes are non-physical and can explode
    // density/linearizations when the estimator state is uninitialized or
    // temporarily inconsistent.
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

    // Clamp Mach, alpha, beta to table ranges to avoid evaluating far
    // outside the aerodynamic grid while preserving saturation behavior.
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
    let v_b_norm = v_b.norm().max(1e-6);
    let mut alpha = v_b.z.atan2(v_b.x).to_degrees();
    let mut beta = (v_b.y / v_b_norm).asin().to_degrees();

    if !table.alpha.is_empty() {
        let a_min = table.alpha[0];
        let a_max = *table.alpha.last().unwrap();
        if alpha < a_min {
            alpha = a_min;
        } else if alpha > a_max {
            alpha = a_max;
        }
    }
    if !table.beta.is_empty() {
        let b_min = table.beta[0];
        let b_max = *table.beta.last().unwrap();
        if beta < b_min {
            beta = b_min;
        } else if beta > b_max {
            beta = b_max;
        }
    }

    let coeff = table.interp(mach, alpha, beta);
    let q_dyn = 0.5 * rho * speed * speed;
    let f_coeff = Vector3::new(coeff.cx, coeff.cy, coeff.cz);
    let mut m_coeff = Vector3::new(coeff.cl, coeff.cm, coeff.cn);
    let inv_2_speed = 0.5 / speed;
    let p_hat = omega_b.x * dyn_params.b_ref * inv_2_speed;
    let q_hat = omega_b.y * dyn_params.c_ref * inv_2_speed;
    let r_hat = omega_b.z * dyn_params.b_ref * inv_2_speed;
    m_coeff.x += coeff.cl_p * p_hat;
    m_coeff.y += coeff.cm_q * q_hat;
    m_coeff.z += coeff.cn_r * r_hat;
    let force_b = q_dyn * dyn_params.s_ref * f_coeff;
    let torque_b = q_dyn
        * dyn_params.s_ref
        * Vector3::new(
            dyn_params.b_ref * m_coeff.x,
            dyn_params.c_ref * m_coeff.y,
            dyn_params.b_ref * m_coeff.z,
        );
    AeroResult {
        force_b,
        torque_b,
        coeff,
        q_dyn,
    }
}
