use super::*;
use rand::Rng;

fn make_run(seed: u64, noise: NoiseParams, cfg: SimConfig) -> SimRun {
    let params = crate::data::load_vehicle_props("data/Vehicle Properties Preview.csv").unwrap();
    let aero = Arc::new(AeroTable::load("data/Aero Coeff Table.csv").unwrap());
    let rcs = Arc::new(RcsTorqueSeries::load("data/RCS Torque Series 100Hz Preview.csv").unwrap());
    SimRun::new(seed, params, noise, aero, rcs, cfg, None)
}

#[test]
fn truth_and_imu_stream_are_identical_across_modes_tuning_and_dropouts() {
    let cfg = SimConfig {
        t_final: 3.0,
        history_sec: 1.0,
        blackout_start: 1.0,
        blackout_end: 2.0,
        star_convergence_count: 0,
        ..SimConfig::default()
    };
    let mut pvt = make_run(11, NoiseParams::tuned(), cfg.clone());
    let mut tight = make_run(
        11,
        NoiseParams::tuned(),
        SimConfig {
            use_tight_coupling: true,
            ..cfg.clone()
        },
    );
    let mut other = make_run(
        11,
        NoiseParams::default(),
        SimConfig {
            dropout_rate: 1.0,
            ..cfg
        },
    );
    other.no_star = true;
    for step in 0..pvt.total_steps {
        pvt.step(step);
        tight.step(step);
        other.step(step);
        for run in [&tight, &other] {
            assert_eq!(pvt.truth.r, run.truth.r);
            assert_eq!(pvt.truth.v, run.truth.v);
            assert_eq!(pvt.truth.q, run.truth.q);
            assert_eq!(pvt.truth.omega, run.truth.omega);
            assert_eq!(pvt.truth.bg, run.truth.bg);
            assert_eq!(pvt.truth.ba, run.truth.ba);
            assert_eq!(pvt.truth.aero_res_force, run.truth.aero_res_force);
            assert_eq!(pvt.truth.aero_res_torque, run.truth.aero_res_torque);
            assert_eq!(pvt.truth.t, run.truth.t);
            assert_eq!(
                pvt.imu_rng.clone().gen::<u64>(),
                run.imu_rng.clone().gen::<u64>()
            );
        }
    }
    assert_eq!(pvt.trajectory, tight.trajectory);
    assert_eq!(pvt.trajectory, other.trajectory);
    assert_eq!(pvt.star_rng.gen::<u64>(), tight.star_rng.gen::<u64>());
    assert_eq!(pvt.radar_rng.gen::<u64>(), tight.radar_rng.gen::<u64>());
    assert_ne!(pvt.est.state.r, tight.est.state.r);
}

#[test]
fn supplied_initial_truth_is_used_and_legacy_ascent_is_reported() {
    let cfg = SimConfig {
        t_final: 0.0,
        ..SimConfig::default()
    };
    let legacy = make_run(0, NoiseParams::tuned(), cfg.clone());
    assert!((legacy.trajectory.initial_radial_velocity_mps - 5821.6).abs() < 1e-9);
    assert_eq!(legacy.trajectory.radar_eligible_samples, 0);
    assert_eq!(legacy.trajectory.outside_aero_table_samples, 1);
    let mut initial = legacy.truth.clone();
    initial.v.x = -100.0;
    let run = make_run(
        0,
        NoiseParams::tuned(),
        SimConfig {
            initial_truth: Some(initial.clone()),
            ..cfg
        },
    );
    assert_eq!(run.truth.r, initial.r);
    assert_eq!(run.truth.v, initial.v);
    assert_eq!(run.truth.bg, initial.bg);
    assert_eq!(run.trajectory.initial_radial_velocity_mps, -100.0);
    assert!(!run.finish().result.passed);
}

#[test]
fn nonfinite_recovery_state_fails_even_when_covariance_is_small() {
    let mut run = make_run(0, NoiseParams::tuned(), SimConfig::default());
    run.est.state.t = 160.0;
    run.est.state.r.x = f64::NAN;
    run.est.l_cov.fill(0.0);
    run.accumulate_recovery_window_stats(64000);
    let result = run.finish().result;
    assert!(!result.recovery_numerically_valid);
    assert!(!result.passed);
}
