use super::{
    history::{HistoryBuffer, PropRecord, Snapshot},
    math::qr_r_upper_in_place,
    Estimator, GnssMeas, Mat2NByN, MatN, RadarMeas, StarMeas,
};
use crate::types::{DynamicsParams, ImuSample, NoiseParams, STATE_DIM};
use nalgebra::{linalg::Cholesky, Matrix3, Matrix6, UnitQuaternion, Vector3};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

#[test]
fn rts_smoother_corrects_attitude() {
    let noise = NoiseParams::default();
    let dyn_params = DynamicsParams::default();
    let mut est = Estimator::new(noise, dyn_params, 5.0);

    let dt = 0.1;
    for k in 1..=10 {
        est.propagate(ImuSample {
            t: k as f64 * dt,
            dt,
            accel_meas: Vector3::zeros(),
            gyro_meas: Vector3::zeros(),
            torque_cmd: Vector3::zeros(),
            force_cmd: Vector3::zeros(),
        });
    }

    let (_xs, ps) = est.rts_smooth(0.3).expect("smoother output");
    let det = ps.determinant();
    assert!(det.is_finite() && det > 0.0);
}

#[test]
fn out_of_order_measurements_rewind_and_replay() {
    let noise = NoiseParams::default();
    let dyn_params = DynamicsParams::default();
    let mut est = Estimator::new(noise, dyn_params, 2.0);

    let dt = 0.0025;
    let mut t = 0.0;
    for _ in 0..32 {
        t += dt;
        est.propagate(ImuSample {
            t,
            dt,
            accel_meas: Vector3::zeros(),
            gyro_meas: Vector3::zeros(),
            torque_cmd: Vector3::zeros(),
            force_cmd: Vector3::zeros(),
        });
    }
    let latest_t = est.state.t;
    assert!((latest_t - 0.08).abs() < 1e-9);

    est.update_gnss(GnssMeas {
        t: 0.04,
        pos: Vector3::new(100.0, -50.0, 25.0),
        vel: Vector3::new(1.0, -0.4, 0.2),
        r_cov: Matrix6::identity() * 4.0,
    });
    assert!((est.state.t - latest_t).abs() < 1e-9);

    est.update_star(StarMeas {
        t: 0.06,
        q_bi: UnitQuaternion::from_scaled_axis(Vector3::new(1e-4, -2e-4, 3e-4)),
        r_cov: Matrix3::identity() * (20e-6 * 20e-6),
    });
    assert!((est.state.t - latest_t).abs() < 1e-9);

    est.update_radar(RadarMeas {
        t: 0.05,
        alt: 1200.0,
        sigma: 1.0,
    });
    assert!((est.state.t - latest_t).abs() < 1e-9);

    est.propagate(ImuSample {
        t: latest_t + dt,
        dt,
        accel_meas: Vector3::zeros(),
        gyro_meas: Vector3::zeros(),
        torque_cmd: Vector3::zeros(),
        force_cmd: Vector3::zeros(),
    });
    assert!((est.state.t - (latest_t + dt)).abs() < 1e-9);
    let _ = Cholesky::new(est.cov()).expect("covariance should remain PD");
}

#[test]
fn out_of_order_raw_gps_rewind_and_replay() {
    let noise = NoiseParams::default();
    let dyn_params = DynamicsParams::default();
    let mut replayed = Estimator::new(noise.clone(), dyn_params.clone(), 2.0);
    let mut chronological = Estimator::new(noise, dyn_params.clone(), 2.0);

    let dt = 0.0025;
    let initial_r = Vector3::new(dyn_params.radius_earth + 60_000.0, 0.0, 0.0);
    let initial_v = Vector3::new(7_600.0, 0.0, 0.0);
    replayed.state.r = initial_r;
    replayed.state.v = initial_v;
    chronological.state.r = initial_r;
    chronological.state.v = initial_v;

    for k in 1..=16 {
        let t = k as f64 * dt;
        chronological.propagate(ImuSample {
            t,
            dt,
            accel_meas: Vector3::zeros(),
            gyro_meas: Vector3::zeros(),
            torque_cmd: Vector3::zeros(),
            force_cmd: Vector3::zeros(),
        });
    }
    let raw_meas = make_raw_gps_meas(chronological.state.r, chronological.state.v, 0.04);
    chronological.update_gps_raw(raw_meas.clone());
    for k in 17..=32 {
        let t = k as f64 * dt;
        chronological.propagate(ImuSample {
            t,
            dt,
            accel_meas: Vector3::zeros(),
            gyro_meas: Vector3::zeros(),
            torque_cmd: Vector3::zeros(),
            force_cmd: Vector3::zeros(),
        });
    }

    for k in 1..=32 {
        let t = k as f64 * dt;
        replayed.propagate(ImuSample {
            t,
            dt,
            accel_meas: Vector3::zeros(),
            gyro_meas: Vector3::zeros(),
            torque_cmd: Vector3::zeros(),
            force_cmd: Vector3::zeros(),
        });
    }
    let latest_t = replayed.state.t;
    replayed.update_gps_raw(raw_meas);

    assert!((replayed.state.t - latest_t).abs() < 1e-9);
    assert!((replayed.state.r - chronological.state.r).norm() < 1.0e-6);
    assert!((replayed.state.v - chronological.state.v).norm() < 1.0e-9);
    let cov_diff = replayed.cov() - chronological.cov();
    assert!(
        cov_diff.norm() < 1.0e-6,
        "covariance mismatch norm={:.3e}",
        cov_diff.norm()
    );
    let _ = Cholesky::new(replayed.cov()).expect("covariance should remain PD after raw replay");
}

#[test]
fn robust_tuning_changes_gnss_gate() {
    let noise = NoiseParams::default();
    let dyn_params = DynamicsParams::default();
    let meas = GnssMeas {
        t: 0.1,
        pos: Vector3::new(0.6, 0.0, 0.0),
        vel: Vector3::zeros(),
        r_cov: Matrix6::identity() * 0.01,
    };

    let mut baseline = Estimator::new(noise.clone(), dyn_params.clone(), 1.0);
    baseline.l_cov = Cholesky::new(MatN::identity() * 0.01).unwrap().l();
    baseline.update_gnss(meas.clone());
    assert_eq!(baseline.last_gnss_accepted, Some((meas.t, true)));

    let mut tightened = Estimator::new(noise, dyn_params, 1.0);
    tightened.l_cov = Cholesky::new(MatN::identity() * 0.01).unwrap().l();
    tightened.set_robust_tuning(0.99, 10.0, 0.0);
    tightened.update_gnss(meas);
    // Tightening the position gate still permits the independently valid
    // velocity observation; the rejected position must not move the state.
    assert_eq!(tightened.last_gnss_accepted, Some((0.1, true)));
    assert!(baseline.state.r.x > 0.2);
    assert_eq!(tightened.state.r.x, 0.0);
}

#[test]
fn iae_enable_updates_and_disable_resets_scales() {
    let noise = NoiseParams::default();
    let dyn_params = DynamicsParams::default();
    let mut est = Estimator::new(noise, dyn_params, 1.0);
    est.l_cov = Cholesky::new(MatN::identity() * 0.01).unwrap().l();
    est.enable_iae(0.35);
    assert!(est.iae_gamma.is_some());

    est.update_gnss(GnssMeas {
        t: 0.1,
        pos: Vector3::new(0.6, 0.0, 0.0),
        vel: Vector3::zeros(),
        r_cov: Matrix6::identity() * 0.01,
    });

    let scale = est.iae_scales[super::robust::ROBUST_SENSOR_GNSS_PVT];
    assert!(scale > 1.0, "expected adaptive scale > 1.0, got {}", scale);

    est.disable_iae();
    assert!(est.iae_gamma.is_none());
    for scale in est.iae_scales {
        assert!((scale - 1.0).abs() < 1.0e-12);
    }
}

#[test]
fn qr_r_upper_in_place_matches_ata() {
    let mut rng = StdRng::seed_from_u64(0x5151_5151);
    let mut a = Mat2NByN::zeros();
    for i in 0..super::STATE_DIM2 {
        for j in 0..STATE_DIM {
            a[(i, j)] = rng.gen_range(-1.0..1.0);
        }
    }
    let a_ref = a;
    let r = qr_r_upper_in_place(&mut a);
    let ata = a_ref.transpose() * a_ref;
    let rtr = r.transpose() * r;

    let mut max_abs = 0.0f64;
    for i in 0..STATE_DIM {
        for j in 0..STATE_DIM {
            max_abs = max_abs.max((ata[(i, j)] - rtr[(i, j)]).abs());
        }
    }
    assert!(max_abs < 1.0e-8, "max_abs={:.3e}", max_abs);
}

#[test]
fn time_update_sqrt_matches_dense_cov_update() {
    let mut rng = StdRng::seed_from_u64(0xDEAD_BEEF);

    let mut a = MatN::zeros();
    let mut b = MatN::zeros();
    let mut phi = MatN::identity();
    for i in 0..STATE_DIM {
        for j in 0..STATE_DIM {
            a[(i, j)] = rng.gen_range(-1.0..1.0);
            b[(i, j)] = rng.gen_range(-1.0..1.0);
            phi[(i, j)] += 0.01 * rng.gen_range(-1.0..1.0);
        }
    }

    let p0 = a * a.transpose() + MatN::identity() * 1.0e-3;
    let q = b * b.transpose() + MatN::identity() * 1.0e-6;
    let l0 = Cholesky::new(p0).unwrap().l();

    let noise = NoiseParams::default();
    let dyn_params = DynamicsParams::default();
    let mut est = Estimator::new(noise, dyn_params, 1.0);
    est.l_cov = l0;

    est.time_update_sqrt(&phi, &q);

    let p_ref = phi * (l0 * l0.transpose()) * phi.transpose() + q;
    let p_new = est.l_cov * est.l_cov.transpose();

    let mut max_abs = 0.0f64;
    for i in 0..STATE_DIM {
        for j in 0..STATE_DIM {
            max_abs = max_abs.max((p_ref[(i, j)] - p_new[(i, j)]).abs());
        }
    }
    assert!(max_abs < 5.0e-8, "max_abs={:.3e}", max_abs);
}

#[test]
fn process_noise_factor_avoids_artificial_bias_growth() {
    let noise = NoiseParams::tuned();
    let dyn_params = DynamicsParams::default();
    let mut est = Estimator::new(noise, dyn_params, 10.0);
    est.state.r = Vector3::new(est.dyn_params.radius_earth + 60_000.0, 0.0, 0.0);
    est.state.v = Vector3::new(7_600.0, 0.0, 0.0);

    let dt = 0.0025;
    let steps = 8_000usize;
    for k in 1..=steps {
        let t = 0.0025 * f64::from(u32::try_from(k).expect("test step count fits in u32"));
        est.propagate(ImuSample {
            t,
            dt,
            accel_meas: Vector3::zeros(),
            gyro_meas: Vector3::zeros(),
            torque_cmd: Vector3::zeros(),
            force_cmd: Vector3::zeros(),
        });
    }

    let p = est.cov();
    let bg_sigma_dps = p[(12, 12)].abs().sqrt().to_degrees();
    let ba_sigma_ug = p[(15, 15)].abs().sqrt() / 9.80665 * 1.0e6;

    assert!(
        bg_sigma_dps < 0.01,
        "gyro-bias sigma exploded: actual={bg_sigma_dps:.6} deg/s expected<0.01 deg/s"
    );
    assert!(
        ba_sigma_ug < 1_000.0,
        "accel-bias sigma exploded: actual={ba_sigma_ug:.3} ug expected<1000 ug"
    );
}

#[test]
fn rts_smoother_uses_predicted_next_state_not_filtered_next_state() {
    let noise = NoiseParams::default();
    let dyn_params = DynamicsParams::default();
    let mut est = Estimator::new(noise, dyn_params, 1.0);

    let l_cov = match Cholesky::new(MatN::identity()) {
        Some(cholesky) => cholesky.l(),
        None => panic!("failed to factor identity covariance in RTS regression setup"),
    };

    let mut snapshot_state_0 = crate::types::State::new();
    snapshot_state_0.t = 0.0;
    snapshot_state_0.r.x = 0.0;

    let mut predicted_state_1 = crate::types::State::new();
    predicted_state_1.t = 1.0;
    predicted_state_1.r.x = 1.0;

    let mut filtered_state_1 = crate::types::State::new();
    filtered_state_1.t = 1.0;
    filtered_state_1.r.x = 11.0;

    let internal_0 = est.snapshot_internal();
    let internal_1 = est.snapshot_internal();

    let mut history = HistoryBuffer::new(1.0, 0.01);
    history.push_snapshot(Snapshot {
        state: snapshot_state_0,
        l_cov,
        internal: internal_0,
    });
    history.push_prop(PropRecord {
        phi: MatN::identity(),
        q: MatN::zeros(),
        predicted_end_state: predicted_state_1,
    });
    history.push_snapshot(Snapshot {
        state: filtered_state_1,
        l_cov,
        internal: internal_1,
    });
    est.history = history;

    let smoothed = match est.rts_smooth(1.0) {
        Some(output) => output,
        None => panic!("RTS smoother unexpectedly returned None on a valid synthetic interval"),
    };
    let smoothed_state = smoothed.0;

    assert!(
        (smoothed_state.r.x - 10.0).abs() < 1.0e-9,
        "RTS used the wrong next-state reference: actual={} expected=10.0",
        smoothed_state.r.x
    );
}

#[test]
fn ten_second_history_stays_within_memory_budget() {
    const MEMORY_BUDGET_BYTES: usize = 20 * 1024 * 1024;

    let noise = NoiseParams::tuned();
    let dyn_params = DynamicsParams::default();
    let mut est = Estimator::new(noise, dyn_params, 10.0);
    est.state.r = Vector3::new(est.dyn_params.radius_earth + 60_000.0, 0.0, 0.0);
    est.state.v = Vector3::new(7_600.0, 0.0, 0.0);

    let dt = 0.0025;
    let gnss_cov = Matrix6::identity() * 4.0;
    let steps = 10usize * 400usize;
    for k in 1..=steps {
        let t = dt * f64::from(u32::try_from(k).expect("test step count fits in u32"));
        est.propagate(ImuSample {
            t,
            dt,
            accel_meas: Vector3::zeros(),
            gyro_meas: Vector3::zeros(),
            torque_cmd: Vector3::zeros(),
            force_cmd: Vector3::zeros(),
        });
        if k % 40 == 1 {
            est.update_gnss(GnssMeas {
                t,
                pos: est.state.r,
                vel: est.state.v,
                r_cov: gnss_cov,
            });
        }
    }

    let health = est.health();
    assert!(
        health.history_capacity_bytes < MEMORY_BUDGET_BYTES,
        "history memory exceeded budget: actual={} bytes expected<{} bytes",
        health.history_capacity_bytes,
        MEMORY_BUDGET_BYTES
    );
}

fn make_raw_gps_meas(r_user: Vector3<f64>, v_user: Vector3<f64>, t: f64) -> super::GpsRawMeas {
    let gps_alt = 20_200_000.0_f64;
    let radius_earth = 6_378_137.0_f64;
    let mu_earth = 3.986004418e14_f64;
    let r_orbit = radius_earth + gps_alt;
    let v_orbit = (mu_earth / r_orbit).sqrt();

    let mut sats = Vec::new();
    for i in 0..8 {
        let phase = (i as f64) * std::f64::consts::TAU / 8.0;
        let sat_pos = Vector3::new(r_orbit * phase.cos(), r_orbit * phase.sin(), 0.0);
        let sat_vel = Vector3::new(-v_orbit * phase.sin(), v_orbit * phase.cos(), 0.0);
        let dr = sat_pos - r_user;
        let rho = dr.norm().max(1.0e-3);
        let los = dr / rho;
        let rr = (sat_vel - v_user).dot(&los);
        sats.push(super::GpsSatMeas {
            pos_ecef: sat_pos,
            vel_ecef: sat_vel,
            pr_m: rho,
            pr_sigma: 1.5,
            range_rate_mps: Some(rr),
            rr_sigma: Some(0.02),
        });
    }

    super::GpsRawMeas { t, sats }
}

#[test]
fn angular_rate_observer_covariance_includes_memory_bias_and_shared_noise() {
    let mut est = Estimator::new(NoiseParams::default(), DynamicsParams::default(), 1.0);
    est.state.r = Vector3::new(est.dyn_params.radius_earth + 60_000.0, 0.0, 0.0);
    let dt = 0.0025;
    let bias_var = 2e-8;
    let mut prior = MatN::identity() * 1e-6;
    for axis in 0..3 {
        prior[(9 + axis, 9 + axis)] = 100.0;
        prior[(12 + axis, 12 + axis)] = bias_var;
    }
    est.l_cov = prior.cholesky().unwrap().l();
    let gyro_density_sq = est.noise.gyro_noise.powi(2);
    let (gain, _) = est.rate_observer_gains(dt);
    est.propagate(ImuSample {
        t: dt,
        dt,
        accel_meas: Vector3::zeros(),
        gyro_meas: Vector3::zeros(),
        torque_cmd: Vector3::zeros(),
        force_cmd: Vector3::zeros(),
    });
    let cov = est.cov();
    for axis in 0..3 {
        let inertia = est.dyn_params.inertia[(axis, axis)];
        let expected_rate_variance = (1.0 - gain).powi(2) * 100.0
            + gain.powi(2) * (bias_var + gyro_density_sq / dt)
            + ((dt * (1.0 - 0.5 * gain)).powi(2) * 1e-6 + est.noise.torque_noise.powi(2) * dt)
                / inertia.powi(2);
        assert!((cov[(9 + axis, 9 + axis)] - expected_rate_variance).abs() < 1e-10);
        assert!((cov[(9 + axis, 12 + axis)] + gain * bias_var).abs() < 1e-12);
        assert!(
            (cov[(6 + axis, 9 + axis)] - gain * (dt * bias_var + gyro_density_sq)).abs() < 1e-12
        );
    }
}

#[test]
fn angular_noise_factor_preserves_all_cross_covariance_blocks() {
    let mut est = Estimator::new(NoiseParams::default(), DynamicsParams::default(), 1.0);
    est.set_blackout_residual(120.0, 160.0, 1000.0, Some(300.0));
    for t in [0.0, 140.0, 161.0] {
        est.state.t = t;
        let q = est.process_noise(0.0025);
        let factor = est.process_noise_factor(&q);
        assert!((q - factor * factor.transpose()).norm() < 1e-10);
    }
}

#[test]
fn star_tracker_recovers_from_a_contaminated_first_fix() {
    let mut est = Estimator::new(NoiseParams::tuned(), DynamicsParams::default(), 10.0);
    est.set_precalibration(false, 120.0);
    est.state.r = Vector3::new(est.dyn_params.radius_earth + 100_000.0, 0.0, 0.0);
    let r_cov = Matrix3::identity() * (20e-6_f64).powi(2);
    est.update_star(super::StarMeas {
        t: 0.0,
        q_bi: UnitQuaternion::from_scaled_axis(Vector3::new(0.05, 0.0, 0.0)),
        r_cov,
    });
    assert!(est.state.q.angle() > 0.04);
    for k in 1..=1000 {
        let t = k as f64 * 0.0025;
        est.propagate(ImuSample {
            t,
            dt: 0.0025,
            gyro_meas: Vector3::zeros(),
            accel_meas: Vector3::zeros(),
            force_cmd: Vector3::zeros(),
            torque_cmd: Vector3::zeros(),
        });
        if k % 200 == 0 {
            est.update_star(super::StarMeas {
                t,
                q_bi: UnitQuaternion::identity(),
                r_cov,
            });
        }
    }
    assert!(
        est.state.q.angle() < 1e-5,
        "attitude was locked to its first outlier"
    );
    assert_eq!(est.star_reject_streak, 0);
    assert!(est.cov().cholesky().is_some());
}

#[test]
fn history_retains_anchor_for_the_entire_ten_second_lag() {
    let mut est = Estimator::new(NoiseParams::tuned(), DynamicsParams::default(), 10.0);
    est.state.r.x = est.dyn_params.radius_earth;
    for k in 1..=4801 {
        est.propagate(ImuSample {
            t: k as f64 * 0.0025,
            dt: 0.0025,
            gyro_meas: Vector3::zeros(),
            accel_meas: Vector3::zeros(),
            force_cmd: Vector3::zeros(),
            torque_cmd: Vector3::zeros(),
        });
    }
    let newest = est.history.snapshots.back().unwrap().state.t;
    let oldest = est.history.snapshots.front().unwrap().state.t;
    assert!(newest - oldest >= 10.0 - 1e-9);
    let (smoothed, _) = est.rts_smooth(10.0).expect("full lag must be retained");
    assert!((newest - smoothed.t - 10.0).abs() < 0.011);
}

#[test]
fn angular_linearization_matches_perturbed_observer_states() {
    fn estimator() -> Estimator {
        let mut est = Estimator::new(NoiseParams::tuned(), DynamicsParams::default(), 1.0);
        est.state.r.x = est.dyn_params.radius_earth + 60_000.0;
        est.state.q = UnitQuaternion::from_scaled_axis(Vector3::new(0.1, -0.2, 0.03));
        est.state.omega = Vector3::new(0.02, -0.01, 0.03);
        est.state.aero_res_torque = Vector3::new(1.0, 2.0, 3.0);
        est
    }
    let dt = 0.0025;
    let imu = ImuSample {
        t: dt,
        dt,
        gyro_meas: Vector3::new(0.02, -0.01, 0.03),
        accel_meas: Vector3::zeros(),
        force_cmd: Vector3::zeros(),
        torque_cmd: Vector3::zeros(),
    };
    let base = estimator();
    let (phi, _) = base.linearize(dt, Vector3::zeros());
    let epsilon = 1e-6;
    for column in (6..15).chain(21..24) {
        let mut plus = estimator();
        let mut minus = estimator();
        let mut error = super::VecN::zeros();
        error[column] = epsilon;
        plus.apply_state_error(error);
        minus.apply_state_error(-error);
        plus.propagate_internal(imu.clone(), false);
        minus.propagate_internal(imu.clone(), false);
        let attitude = (plus.state.q * minus.state.q.inverse()).scaled_axis() / (2.0 * epsilon);
        let rate = (plus.state.omega - minus.state.omega) / (2.0 * epsilon);
        for axis in 0..3 {
            assert!(
                (attitude[axis] - phi[(6 + axis, column)]).abs() < 1e-6,
                "attitude Jacobian column {column}"
            );
            assert!(
                (rate[axis] - phi[(9 + axis, column)]).abs() < 1e-6,
                "rate Jacobian column {column}"
            );
        }
    }
}

#[test]
fn gnss_position_outlier_keeps_an_independently_valid_velocity() {
    let mut est = Estimator::new(NoiseParams::tuned(), DynamicsParams::default(), 1.0);
    est.state.r.x = est.dyn_params.radius_earth + 60_000.0;
    let original_r = est.state.r;
    let mut p = MatN::identity() * 1e-6;
    for i in 0..3 {
        p[(i, i)] = 1.0;
        p[(i + 3, i + 3)] = 0.01;
    }
    est.l_cov = p.cholesky().unwrap().l();
    let mut r_cov = Matrix6::identity() * 3.61;
    for i in 3..6 {
        r_cov[(i, i)] = 0.0064;
    }
    est.update_gnss(GnssMeas {
        t: 0.0,
        pos: original_r + Vector3::new(1000.0, 0.0, 0.0),
        vel: Vector3::new(0.1, 0.0, 0.0),
        r_cov,
    });
    assert!(est.last_gnss_base_nis.unwrap().1 > 1000.0);
    assert_eq!(est.last_gnss_accepted, Some((0.0, true)));
    assert_eq!(est.state.r, original_r);
    assert!(est.state.v.x > 0.05 && est.state.v.x < 0.1);
    assert!(est.cov()[(3, 3)] < 0.01);
}

#[test]
fn raw_acquisition_recovers_a_prior_tail_but_not_inconsistent_ranges() {
    let truth_r = Vector3::new(6_438_137.0, 0.0, 3_000_000.0);
    let truth_v = Vector3::new(0.0, 3_000.0, 0.0);
    let mut est = Estimator::new(NoiseParams::tuned(), DynamicsParams::default(), 1.0);
    est.state.r = truth_r + Vector3::new(30_000.0, 0.0, 0.0);
    est.state.v = truth_v;
    let mut prior = MatN::identity() * 1e-6;
    for i in 0..3 {
        prior[(i, i)] = 5_000.0_f64.powi(2);
        prior[(i + 3, i + 3)] = 50.0_f64.powi(2);
    }
    est.l_cov = prior.cholesky().unwrap().l();
    let clean = make_raw_gps_meas(truth_r, truth_v, 0.0);
    est.update_gps_raw(clean.clone());
    assert_eq!(est.last_gnss_accepted, Some((0.0, false)));
    let mut bad = clean.clone();
    bad.sats[0].pr_m += 10_000.0;
    for _ in 0..12 {
        est.update_gps_raw(bad.clone());
    }
    assert_eq!(est.state.r, truth_r + Vector3::new(30_000.0, 0.0, 0.0));
    assert!(
        (est.cov() - prior).norm() < 1e-7,
        "inconsistent ranges widened the prior"
    );
    for _ in 0..9 {
        est.update_gps_raw(clean.clone());
    }
    assert_eq!(est.last_gnss_accepted, Some((0.0, true)));
    assert!((est.state.r - truth_r).norm() < 1.0);
    assert!((est.state.v - truth_v).norm() < 0.01);
    assert!(est.cov().cholesky().is_some());
}
