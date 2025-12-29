use entry_blackout_estimation::{
    data::{AeroCoeff, AeroTable, load_vehicle_props},
    estimator::{Estimator, GnssMeas, RadarMeas, StarMeas},
    types::{DynamicsParams, NoiseParams, ImuSample},
    physics::isa_atmosphere,
};
use nalgebra::{linalg::Cholesky, Matrix3, Matrix6, UnitQuaternion, Vector3};

#[test]
fn load_data_and_interp() {
    let _dyn = load_vehicle_props("data/Vehicle Properties Preview.csv").unwrap();
    let aero = AeroTable::load("data/Aero Coeff Table.csv").unwrap();
    let c = aero.interp(5.0, 10.0, 0.0);
    assert!(c.cx > 0.0);
    let (rho, a) = isa_atmosphere(20_000.0);
    assert!(rho > 0.0);
    assert!(a > 200.0);
}

#[test]
fn ekf_propagate_update_and_smooth() {
    let noise = NoiseParams::default();
    let dyn_params = DynamicsParams::default();
    let aero = AeroTable::load("data/Aero Coeff Table.csv").unwrap();
    let mut est = Estimator::new(noise, dyn_params, 5.0);
    est.set_aero_table(aero);

    let dt = 0.01;
    for k in 1..=100 {
        est.propagate(ImuSample {
            t: k as f64 * dt,
            dt,
            accel_meas: Vector3::zeros(),
            gyro_meas: Vector3::zeros(),
            torque_cmd: Vector3::zeros(),
            force_cmd: Vector3::zeros(),
        });
    }
    let r_cov = Matrix6::identity() * 4.0;
    est.update_gnss(GnssMeas {
        t: 0.5,
        pos: Vector3::zeros(),
        vel: Vector3::zeros(),
        r_cov,
    });
    est.update_star(StarMeas {
        t: 0.4,
        q_bi: UnitQuaternion::identity(),
        r_cov: Matrix3::identity() * 1e-4,
    });
    est.update_radar(RadarMeas {
        t: 0.6,
        alt: 0.0,
        sigma: 2.0,
    });

    // covariance should stay PD
    let p = est.cov();
    let chol = Cholesky::new(p);
    if chol.is_none() {
        let mut nan = 0usize;
        let mut inf = 0usize;
        let mut min_diag = p[(0, 0)];
        let mut max_diag = p[(0, 0)];
        for i in 0..entry_blackout_estimation::types::STATE_DIM {
            for j in 0..entry_blackout_estimation::types::STATE_DIM {
                let v = p[(i, j)];
                if v.is_nan() {
                    nan += 1;
                } else if !v.is_finite() {
                    inf += 1;
                }
            }
            min_diag = min_diag.min(p[(i, i)]);
            max_diag = max_diag.max(p[(i, i)]);
        }
        eprintln!(
            "[smoke] covariance not PD: nan={} inf={} min_diag={:.3e} max_diag={:.3e}",
            nan, inf, min_diag, max_diag
        );
    }
    assert!(chol.is_some());

    // smoother should return something
    let smoothed = est.rts_smooth(0.2);
    assert!(smoothed.is_some());
}

#[test]
fn aero_interp_clamped_no_double_count() {
    let coeff = AeroCoeff {
        cx: 1.0,
        cy: 2.0,
        cz: 3.0,
        cl: 4.0,
        cm: 5.0,
        cn: 6.0,
        cl_p: 7.0,
        cm_q: 8.0,
        cn_r: 9.0,
    };
    let table = AeroTable {
        mach: vec![0.0],
        alpha: vec![0.0],
        beta: vec![0.0],
        coeffs: vec![coeff.clone()],
        n_a: 1,
        n_b: 1,
    };

    let c = table.interp(10.0, 20.0, -30.0);
    assert!((c.cx - coeff.cx).abs() < 1e-12);
    assert!((c.cy - coeff.cy).abs() < 1e-12);
    assert!((c.cz - coeff.cz).abs() < 1e-12);
    assert!((c.cl - coeff.cl).abs() < 1e-12);
    assert!((c.cm - coeff.cm).abs() < 1e-12);
    assert!((c.cn - coeff.cn).abs() < 1e-12);
    assert!((c.cl_p - coeff.cl_p).abs() < 1e-12);
    assert!((c.cm_q - coeff.cm_q).abs() < 1e-12);
    assert!((c.cn_r - coeff.cn_r).abs() < 1e-12);
}
