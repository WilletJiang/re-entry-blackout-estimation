# Interface Specification

This note defines the flight style integration contract for the estimator. It complements `tasks.md` and `design.md` by focusing on call semantics, timestamp discipline, and the small set of rules that must be respected if rewind and replay are to remain valid.

## Timestamp Discipline

Every measurement passed to the estimator shall carry the physical measurement time `t` in seconds. Arrival time shall not replace it. If a star tracker sample reaches the estimator one hundred and twenty milliseconds after exposure, the caller shall submit `meas.t = t_now - 0.12`. Submitting `meas.t = t_now` would be incorrect. The same rule applies to radar and any delayed raw GNSS packet. The estimator will rewind within its history horizon and apply the update at `meas.t`.

## Construction and Configuration

An estimator instance is created with the standard constructor shown below.

```rust
let mut est = Estimator::new(noise, dyn_params, history_sec);
```

`history_sec` sets the fixed lag horizon and should be chosen long enough to cover the maximum expected delay together with the smoothing window. Optional runtime configuration is then applied through the existing setter methods. The most relevant controls are `set_blackout_residual`, `set_precalibration`, `set_robust_tuning`, `enable_iae`, `disable_iae`, `enable_aero_learning`, and `disable_aero_learning`. These calls modify estimator behaviour prospectively. They do not retroactively alter the already recorded history.

## Main Loop Contract

The estimator is driven by a chronological IMU propagation loop. Each IMU sample should be submitted at the IMU rate through `propagate(ImuSample { ... })`. Measurement updates may then be injected whenever data become available, including delayed arrivals, provided that each payload carries the correct measurement timestamp. A minimal integration pattern is shown below.

```rust
let mut est = Estimator::new(noise, dyn_params, history_sec);

est.set_precalibration(true, blackout_start);
est.set_robust_tuning(gate_p, r_scale_cap, student_nu);

for imu in imu_stream {
    est.propagate(imu);

    if let Some(gnss) = maybe_gnss() {
        est.update_gnss(gnss);
    }
    if let Some(star) = maybe_star() {
        est.update_star(star);
    }
    if let Some(radar) = maybe_radar() {
        est.update_radar(radar);
    }
}
```

When an aerodynamic table is configured, `torque_cmd` and `force_cmd` are consumed by the same nominal aero and RCS composition path used by truth propagation. This keeps the nominal dynamics seen by simulation and estimation consistent.

## Measurement Submission Contract

The supported primary calls are `update_gnss`, `update_star`, `update_radar`, and `update_gps_raw_ref`. `update_gnss` expects a `GnssMeas` carrying timestamp, position, velocity, and covariance. `update_star` expects a `StarMeas` carrying timestamp, body to inertial quaternion, and covariance. `update_radar` expects a `RadarMeas` carrying timestamp, altitude, and standard deviation. These calls are replay aware and may be issued with delayed timestamps as long as the referenced time remains inside the history horizon.

## Tight Coupled Raw GNSS Contract

The raw GNSS path is exposed through `update_gps_raw_ref(t, sats)`. Each satellite element shall provide ECEF position, ECEF velocity, pseudorange, pseudorange standard deviation, range rate, and range rate standard deviation. A typical call shape is shown below.

```rust
est.update_gps_raw_ref(t, &[
    GpsSatMeas {
        pos_ecef,
        vel_ecef,
        pr_m,
        pr_sigma,
        range_rate_mps,
        rr_sigma,
    },
]);
```

This interface is part of the core estimator and is available in `flight` builds as well as diagnostics builds. Raw packets follow the same rewind and replay semantics as GNSS PVT, star tracker, and radar. A delayed raw packet shall be submitted with its original measurement timestamp.

## Replay Safe External Interventions

Some supervisory actions change state or covariance in ways that must survive history rewind. These actions shall use replay safe interfaces. Use `inflate_rv_record` instead of `inflate_rv`, and use `apply_bias_calibration_record` instead of `apply_bias_calibration`. The replay safe forms both mutate the live state and record an event so that deterministic replay reconstructs the same result after a rewind. Calling the non recorded forms from an external supervisor creates a hidden divergence between the visible current state and the recoverable history state.

## Health Reporting

`Estimator::health()` exposes the most useful operational counters. It reports the last GNSS NIS value and accept decision, the recent GNSS rejection streak, and the approximate history buffer usage in bytes. This interface serves runtime monitoring and does not alter the estimation state.

The integration rule is therefore simple. Propagate strictly in time, submit every measurement with its physical timestamp, and route any externally imposed covariance or bias intervention through a replay safe API. If those three conditions are met, the estimator can preserve both fixed lag smoothing and delayed measurement consistency without any special handling in the caller.
