## Interface Notes (Flight-style usage)

### Timestamp contract

All measurements must carry their **measurement timestamp** `t` (seconds), not “arrival time”.

If a sensor is delayed by a fixed amount (e.g. star tracker 120 ms), the correct call is:

- `meas.t = t_now - delay`

The estimator will **rewind/replay** within its history horizon to apply the update at `meas.t`.

### Main loop

Typical usage:

1. Construct:

- `Estimator::new(noise, dyn_params, history_sec)`

2. Configure (optional):

- `set_blackout_residual(start, end, ...)`
- `set_precalibration(enabled, blackout_start)`
- `enable_aero_learning(rate)` / `disable_aero_learning()`

3. Run:

- IMU: call `propagate(ImuSample { t, dt, accel_meas, gyro_meas, torque_cmd, force_cmd })` at IMU rate
- GNSS PVT: call `update_gnss(GnssMeas { t, pos, vel, r_cov })`
- Star: call `update_star(StarMeas { t, q_bi, r_cov })`
- Radar: call `update_radar(RadarMeas { t, alt, sigma })`

### Tight-coupled GNSS

For raw pseudorange + Doppler:

- `update_gps_raw_ref(t, sats)` where each sat provides `{pos_ecef, vel_ecef, pr_m, pr_sigma, range_rate_mps, rr_sigma}`

This path is available in `flight` builds (not diagnostics-gated).

### Replay-safe “external interventions”

If an external supervisor needs to change estimator state/covariance in a way that must survive rewind/replay, use:

- `inflate_rv_record(pos_scale, vel_scale)` instead of `inflate_rv(...)`
- `apply_bias_calibration_record()` instead of `apply_bias_calibration()`

These are recorded as events and checkpointed so deterministic replay remains valid.

### Health / counters

Use `Estimator::health()` to get:

- last GNSS NIS and accept flag
- GNSS reject streak
- history buffer sizes and approximate capacity bytes

