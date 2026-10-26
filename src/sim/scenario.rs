//! Explicit synthetic reference case. These coefficients are a dissipative
//! surrogate, not an extension or calibration of the ten-row preview table.
use crate::{
    data::{AeroCoeff, AeroTable},
    types::{DynamicsParams, State},
};
use nalgebra::{UnitQuaternion, Vector3};

pub fn reference_entry(params: &DynamicsParams) -> (State, AeroTable) {
    let mut state = State::new();
    state.r = Vector3::new(params.radius_earth + 141_000.0, 0.0, 0.0);
    let gamma = (-2.0_f64).to_radians();
    state.v = Vector3::new(3_000.0 * gamma.sin(), 0.0, 3_000.0 * gamma.cos());
    state.q = UnitQuaternion::rotation_between(&-Vector3::z(), &state.v.normalize()).unwrap();
    let mach = vec![1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0];
    let alpha: Vec<f64> = (-180..=180).step_by(5).map(f64::from).collect();
    let beta: Vec<f64> = (-90..=90).step_by(5).map(f64::from).collect();
    let mut coeffs = Vec::with_capacity(mach.len() * alpha.len() * beta.len());
    for _ in &mach {
        for &angle in &alpha {
            for &sideslip in &beta {
                let a = angle.to_radians();
                let b = sideslip.to_radians();
                let drag = 1.2 + 0.8 * a.sin().powi(2);
                coeffs.push(AeroCoeff {
                    cx: -drag * a.cos() * b.cos(),
                    cy: -drag * b.sin(),
                    cz: -drag * a.sin() * b.cos(),
                    cl: 0.05 * b.sin(),
                    cm: -0.05 * (a + std::f64::consts::FRAC_PI_2).sin(),
                    cn: 0.0,
                    cl_p: -0.5,
                    cm_q: -11.0,
                    cn_r: -1.4,
                });
            }
        }
    }
    let table = AeroTable {
        n_a: alpha.len(),
        n_b: beta.len(),
        mach,
        alpha,
        beta,
        coeffs,
    };
    (state, table)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reference_case_descends_and_aerodynamic_force_dissipates_energy() {
        let params = DynamicsParams::default();
        let (state, table) = reference_entry(&params);
        assert!(state.r.dot(&state.v) < 0.0);
        for alpha in [-170.0_f64, -91.0, -45.0, 0.0, 70.0, 173.0] {
            for beta in [-50.0_f64, 0.0, 51.0] {
                let a = alpha.to_radians();
                let b = beta.to_radians();
                let velocity = Vector3::new(a.cos() * b.cos(), b.sin(), a.sin() * b.cos());
                let coeff = table.interp(7.0, alpha, beta);
                assert!(Vector3::new(coeff.cx, coeff.cy, coeff.cz).dot(&velocity) < 0.0);
            }
        }
    }
}
