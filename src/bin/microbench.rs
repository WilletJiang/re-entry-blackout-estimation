#[cfg(feature = "diagnostics")]
use entry_blackout_estimation::estimator::GpsSatMeas;
use entry_blackout_estimation::estimator::{Estimator, GnssMeas};
use entry_blackout_estimation::types::{DynamicsParams, ImuSample, NoiseParams};
use nalgebra::{Matrix6, UnitQuaternion, Vector3};
use std::time::Instant;

const PR_SIGMA_M: f64 = 1.5;
const RR_SIGMA_MPS: f64 = 0.02;

fn ns_per_iter(elapsed: std::time::Duration, iters: usize) -> f64 {
    (elapsed.as_secs_f64() * 1.0e9) / (iters as f64).max(1.0)
}

fn main() -> anyhow::Result<()> {
    let noise = NoiseParams::tuned();
    let dyn_params = DynamicsParams::default();
    let dt = 0.0025;

    bench_propagate(&noise, &dyn_params, dt);
    #[cfg(feature = "diagnostics")]
    bench_time_update(&noise, &dyn_params, dt);
    bench_update_gnss(&noise, &dyn_params);
    #[cfg(feature = "diagnostics")]
    bench_update_gps_raw(&noise, &dyn_params);
    Ok(())
}

fn bench_propagate(noise: &NoiseParams, dyn_params: &DynamicsParams, dt: f64) {
    let mut est = Estimator::new(noise.clone(), dyn_params.clone(), 1.0);
    est.state.r = Vector3::new(dyn_params.radius_earth + 60_000.0, 0.0, 0.0);
    est.state.v = Vector3::new(7_600.0, 0.0, 0.0);
    est.state.q = UnitQuaternion::identity();

    let warm = 2_000usize;
    let iters = 20_000usize;
    for k in 1..=warm {
        est.propagate(ImuSample {
            t: k as f64 * dt,
            dt,
            accel_meas: Vector3::zeros(),
            gyro_meas: Vector3::zeros(),
            torque_cmd: Vector3::zeros(),
            force_cmd: Vector3::zeros(),
        });
    }
    let t0 = est.state.t;
    let start = Instant::now();
    for k in 1..=iters {
        est.propagate(ImuSample {
            t: t0 + k as f64 * dt,
            dt,
            accel_meas: Vector3::zeros(),
            gyro_meas: Vector3::zeros(),
            torque_cmd: Vector3::zeros(),
            force_cmd: Vector3::zeros(),
        });
    }
    let elapsed = start.elapsed();
    println!(
        "propagate(): iters={} ns/iter={:.1}",
        iters,
        ns_per_iter(elapsed, iters)
    );
}

#[cfg(feature = "diagnostics")]
fn bench_time_update(noise: &NoiseParams, dyn_params: &DynamicsParams, dt: f64) {
    let mut est = Estimator::new(noise.clone(), dyn_params.clone(), 1.0);
    let mut phi = entry_blackout_estimation::estimator::MatN::identity();
    for i in 0..3 {
        phi[(i, i + 3)] = dt;
    }
    let qk = entry_blackout_estimation::estimator::MatN::identity() * 1.0e-3;

    let warm = 5_000usize;
    let iters = 50_000usize;
    for _ in 0..warm {
        est.bench_time_update_sqrt(&phi, &qk);
    }
    let start = Instant::now();
    for _ in 0..iters {
        est.bench_time_update_sqrt(&phi, &qk);
    }
    let elapsed = start.elapsed();
    println!(
        "time_update_sqrt(): iters={} ns/iter={:.1}",
        iters,
        ns_per_iter(elapsed, iters)
    );
}

fn bench_update_gnss(noise: &NoiseParams, dyn_params: &DynamicsParams) {
    let mut est = Estimator::new(noise.clone(), dyn_params.clone(), 1.0);
    est.state.r = Vector3::new(dyn_params.radius_earth + 60_000.0, 0.0, 0.0);
    est.state.v = Vector3::new(7_600.0, 0.0, 0.0);
    est.state.q = UnitQuaternion::identity();
    let r_cov = Matrix6::identity() * 4.0;
    let warm = 2_000usize;
    let iters = 20_000usize;
    for k in 0..warm {
        est.update_gnss(GnssMeas {
            t: k as f64 * 0.1,
            pos: est.state.r,
            vel: est.state.v,
            r_cov: r_cov.clone(),
        });
    }
    let start = Instant::now();
    for k in 0..iters {
        est.update_gnss(GnssMeas {
            t: (warm + k) as f64 * 0.1,
            pos: est.state.r,
            vel: est.state.v,
            r_cov: r_cov.clone(),
        });
    }
    let elapsed = start.elapsed();
    println!(
        "update_gnss(): iters={} ns/iter={:.1}",
        iters,
        ns_per_iter(elapsed, iters)
    );
}

#[cfg(feature = "diagnostics")]
fn bench_update_gps_raw(noise: &NoiseParams, dyn_params: &DynamicsParams) {
    let mut est = Estimator::new(noise.clone(), dyn_params.clone(), 1.0);
    est.state.r = Vector3::new(dyn_params.radius_earth + 60_000.0, 0.0, 0.0);
    est.state.v = Vector3::new(7_600.0, 0.0, 0.0);
    est.state.q = UnitQuaternion::identity();

    let gps_alt = 20_200_000.0;
    let r_orbit = dyn_params.radius_earth + gps_alt;
    let v_orbit = (dyn_params.mu_earth / r_orbit).sqrt();

    let mut sats = Vec::new();
    for i in 0..8 {
        let phase = (i as f64) * std::f64::consts::TAU / 8.0;
        let sat_pos = Vector3::new(r_orbit * phase.cos(), r_orbit * phase.sin(), 0.0);
        let sat_vel = Vector3::new(-v_orbit * phase.sin(), v_orbit * phase.cos(), 0.0);
        let dr = sat_pos - est.state.r;
        let rho = dr.norm().max(1e-3);
        let los = dr / rho;
        let pr = rho;
        let rr = (sat_vel - est.state.v).dot(&los);
        sats.push(GpsSatMeas {
            pos_ecef: sat_pos,
            vel_ecef: sat_vel,
            pr_m: pr,
            pr_sigma: PR_SIGMA_M,
            range_rate_mps: Some(rr),
            rr_sigma: Some(RR_SIGMA_MPS),
        });
    }
    let warm = 1_000usize;
    let iters = 5_000usize;
    for k in 0..warm {
        est.update_gps_raw_ref(k as f64 * 0.1, &sats);
    }
    let start = Instant::now();
    for k in 0..iters {
        est.update_gps_raw_ref((warm + k) as f64 * 0.1, &sats);
    }
    let elapsed = start.elapsed();
    println!(
        "update_gps_raw(): iters={} ns/iter={:.1}",
        iters,
        ns_per_iter(elapsed, iters)
    );
}
