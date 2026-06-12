#[cfg(feature = "diagnostics")]
use entry_blackout_estimation::{
    data::RcsTorqueSeries,
    sim::{simulate_single_with_trace, BlackoutTraceConfig, MonteCarloTraceResult, SimConfig},
};
use entry_blackout_estimation::{
    data::{load_vehicle_props, AeroCoeff, AeroTable},
    estimator::{Estimator, GnssMeas, RadarMeas, StarMeas},
    physics::isa_atmosphere,
    types::{DynamicsParams, ImuSample, NoiseParams},
};
use nalgebra::{linalg::Cholesky, Matrix3, Matrix6, UnitQuaternion, Vector3};
#[cfg(feature = "diagnostics")]
use std::sync::Arc;

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

#[cfg(feature = "diagnostics")]
#[test]
fn pvt_recovery_window_velocity_sigma_stays_inside_acceptance() {
    let trace_cfg = recovery_trace_config();
    let trace_out = run_trace_case(0, false, trace_cfg);
    let vel_sigma3_max_mps = max_trace_value_in_window(
        &trace_out.trace.vel_sigma3_mps,
        trace_cfg,
        trace_cfg.blackout_end,
    );

    assert!(trace_out.result.passed);
    assert!(vel_sigma3_max_mps <= 0.5 + 1.0e-3);
}

#[cfg(feature = "diagnostics")]
#[test]
fn tight_seed_zero_stays_bounded_through_blackout_exit() {
    let trace_cfg = recovery_trace_config();
    let trace_out = run_trace_case(0, true, trace_cfg);
    let max_pos_err_m = max_trace_value_in_window(
        &trace_out.trace.pos_err_m,
        trace_cfg,
        trace_cfg.blackout_end,
    );

    assert!(trace_out.result.passed);
    assert!(max_pos_err_m <= 10.0);
}

#[cfg(feature = "diagnostics")]
#[test]
fn tight_seed_five_stays_bounded_through_recovery_window() {
    let trace_cfg = recovery_trace_config();
    let trace_out = run_trace_case(5, true, trace_cfg);
    let max_pos_err_m = max_trace_value_in_window(
        &trace_out.trace.pos_err_m,
        trace_cfg,
        trace_cfg.blackout_end,
    );

    assert!(trace_out.result.passed);
    assert!(max_pos_err_m <= 10.0);
}

#[cfg(feature = "diagnostics")]
fn recovery_trace_config() -> BlackoutTraceConfig {
    BlackoutTraceConfig {
        t_start: 159.995,
        t_end: 162.0,
        dt: 0.0025,
        blackout_start: 120.0,
        blackout_end: 160.0,
    }
}

#[cfg(feature = "diagnostics")]
fn run_trace_case(
    seed: u64,
    use_tight_coupling: bool,
    trace_cfg: BlackoutTraceConfig,
) -> MonteCarloTraceResult {
    let dyn_params = load_vehicle_props("data/Vehicle Properties Preview.csv").unwrap();
    let aero = Arc::new(AeroTable::load("data/Aero Coeff Table.csv").unwrap());
    let rcs = Arc::new(RcsTorqueSeries::load("data/RCS Torque Series 100Hz Preview.csv").unwrap());
    let noise = NoiseParams::tuned();
    let mut cfg = SimConfig::default();
    cfg.use_tight_coupling = use_tight_coupling;
    simulate_single_with_trace(seed, dyn_params, noise, aero, rcs, cfg, trace_cfg)
}

#[cfg(feature = "diagnostics")]
fn max_trace_value_in_window(
    samples: &[f32],
    trace_cfg: BlackoutTraceConfig,
    window_start_s: f64,
) -> f32 {
    let mut max_value = f32::NEG_INFINITY;
    for (sample_idx, &sample) in samples.iter().enumerate() {
        if trace_cfg.time_at(sample_idx) + trace_cfg.dt < window_start_s {
            continue;
        }
        if sample > max_value {
            max_value = sample;
        }
    }
    max_value
}
