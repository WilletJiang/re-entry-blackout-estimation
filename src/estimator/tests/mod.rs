use super::{math::qr_r_upper_in_place, Estimator, GnssMeas, Mat2NByN, MatN, RadarMeas, StarMeas};
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
    assert_eq!(tightened.last_gnss_accepted, Some((0.1, false)));
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
