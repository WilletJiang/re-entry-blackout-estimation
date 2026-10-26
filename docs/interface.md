# Interface Specification

This note defines the flight style integration contract for the estimator. It complements `tasks.md` and `design.md` by focusing on call semantics, timestamp discipline, and the small set of rules that must be respected if rewind and replay are to remain valid.

## Timestamp Discipline

Every measurement passed to the estimator shall carry the physical measurement time `t` in seconds. Arrival time shall not replace it. If a star tracker sample reaches the estimator one hundred and twenty milliseconds after exposure, the caller shall submit `meas.t = t_now - 0.12`. Submitting `meas.t = t_now` would be incorrect. The same rule applies to radar and any delayed raw GNSS packet. The estimator will rewind within its history horizon and apply the update at `meas.t`.

## Construction and Configuration

An estimator instance is created with the standard constructor shown below.

```rust
let mut est = Estimator::new(noise, dyn_params, history_sec);
```

`history_sec` sets the fixed lag horizon and should be chosen long enough to cover the maximum expected delay together with the smoothing window. For the task in `tasks.md`, the diagnostics harness uses ten seconds by default so that delayed measurement handling, smoothing, and memory accounting are evaluated on the required horizon. Optional runtime configuration is then applied through the existing setter methods. The most relevant controls are `set_blackout_residual`, `set_precalibration`, `set_robust_tuning`, `enable_iae`, `disable_iae`, `enable_aero_learning`, and `disable_aero_learning`. These calls modify estimator behaviour prospectively. They do not retroactively alter the already recorded history.

## Main Loop Contract

The estimator is driven by a chronological IMU propagation loop. Each IMU sample should be submitted at the IMU rate through `propagate(ImuSample { ... })`. Measurement updates may then be injected whenever data become available, including delayed arrivals, provided that each payload carries the correct measurement timestamp. A minimal integration pattern is shown below.

```rust
let mut est = Estimator::new(noise, dyn_params, history_sec);

// Continuous random-walk bias estimation is the default.
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

The raw GNSS path is exposed through `update_gps_raw_ref(t, sats)`. Each satellite element shall provide position and velocity in the same inertial frame as the navigation state, clock-corrected pseudorange/range rate, and their standard deviations. The public fields retain the legacy names `pos_ecef` and `vel_ecef`; those names do not perform a coordinate transformation. Convert Earth-fixed ephemerides before submitting them to an inertial estimator. The reference simulator uses an idealized inertial constellation and does not model receiver clock states or atmospheric propagation corrections. A typical call shape is shown below.

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

`Estimator::health()` exposes the most useful operational counters. It reports the last GNSS NIS value and accept decision, the recent GNSS rejection streak, and the approximate history buffer usage in bytes. This interface serves runtime monitoring and does not alter the estimation state. The diagnostics harness records the peak history counters and memory estimate from this interface into its summary outputs.

The integration rule is therefore simple. Propagate strictly in time, submit every measurement with its physical timestamp, and route any externally imposed covariance or bias intervention through a replay safe API. If those three conditions are met, the estimator can preserve both fixed lag smoothing and delayed measurement consistency without any special handling in the caller.

## Noise Units and Angular-Rate Output

`gyro_noise` and `accel_noise` are continuous white-noise densities defined by
`Var(integrated noise over dt) = density^2 * dt`. They are not per-sample standard
deviations. Bias random-walk fields likewise specify a variance rate.
`State::omega` is a causal rate/torque-observer output, with nominal rigid-body
feedforward and configurable `rate_observer_bandwidth` (rad/s; default 2).
Its covariance includes rate memory, bias uncertainty and the shared gyro-noise
correlation with attitude and torque. Attitude still integrates the gyro directly.
See [design.md](design.md) for the observer equations and bandwidth convention.

`last_gnss_nis` retains the robust diagnostic. `last_gnss_base_nis`,
`last_star_base_nis`, and `last_radar_base_nis` expose pre-robustification NIS.
Read diagnostics immediately after an update; delayed replay may replace later
diagnostics. These values are monitoring outputs, not replay state.

The reference batch disables constant-bias precalibration. Explicitly enabling
that legacy path assumes a constant calibration bias and does not establish
consistency for a random-walk bias. The old `RobustTuning::reacq_*_scale` fields
remain source-compatible but no longer cause an artificial blackout-exit reset.

## Streaming retrospective estimates

`visit_smoothed(lag_sec, callback)` visits stored checkpoints newest first using
only the observations already assimilated. The callback receives references to a
state and covariance valid for that call. It allocates no history-sized workspace.
On `None`, discard any partial callback results. `rts_smooth(lag_sec)` retains its
single-oldest-state interface. These calls are synchronous: the target scheduler
must budget their measured execution time separately from ordinary propagation.
The deadline harness uses one two-second pass at simulation time 162 s; it does
not insert this analysis pass into the flight-style 100 Hz resource cycle.
