## Design Notes — Entry Blackout Estimation (SE-2317)

This repo targets the acceptance criteria in `tasks.md`: survive a GNSS blackout and recover within the **post-blackout 2 s** window, while keeping **NEES/NIS** statistically consistent and meeting strict runtime/memory constraints.

### 1) State & filter form

The estimator is a **24-state error-state SR-EKF**:

- `r (3), v (3), δθ (3), ω (3), bg (3), ba (3), aero_res_force (3), aero_res_torque (3)`

Propagation is strapdown-style using IMU specific force and gyro rate, plus gravity (and optional learned aero correction). The nominal aerodynamic / RCS composition now comes from the same helper on both the truth and filter sides; the filter still uses IMU specific force for state integration so nominal terms are shared without double-counting measured acceleration. Covariance is propagated in **square-root form** (QR) for numerical stability.

### 2) Process noise (covariance realism)

`process_noise()` in `src/estimator/support/noise.rs` is phase-aware:

- phase-dependent `q_scale_rv_{pre,blk,post}` and `q_scale_ba_{pre,blk,post}`
- OU-style aero residual force/torque with optional strengthened blackout values
- additional attitude/PV growth during star-tracker saturation (matches the task spec window)

### 3) Measurement updates & robustness

Measurement updates are gated using χ² (NIS) and use Joseph-form covariance updates with an explicit SPD repair (`enforce_pd` + Cholesky/jitter/spectral fallback). Runtime `RobustTuning` now connects the sweep knobs to live behavior:

- `gate_p` overrides the effective χ² gate level
- `student_nu` enables Student-t-style covariance inflation
- `r_scale_cap` bounds covariance inflation
- `iae_gamma` enables adaptive innovation-based covariance scaling (IAE) with replay-safe internal state

Key choices:

- **GNSS PVT** updates use an “active-set” update (only `r,v`) to prevent cross-covariance from yanking attitude/bias states.
- **Star tracker** is attitude-only (conservative by design).
- **Radar altimeter** is radial height only.

### 4) Time alignment (delays / out-of-order)

The spec contains delayed sensors (star/radar). To keep the estimator time-consistent, the implementation uses:

**rewind-to-measurement-time + deterministic replay** (within the fixed-lag history buffer).

Mechanism (conceptually):

1. Restore the closest snapshot `≤ t_meas`.
2. Replay logged IMU/measurement events up to `t_meas`.
3. Apply the update at `t_meas`.
4. Replay the remaining events back to “now”.

Critical invariant: **every state/covariance change that must survive a rewind must be either**

- produced by replayed events, or
- checkpointed into a snapshot at the same timestamp.

For external callers, use replay-safe APIs:

- `Estimator::inflate_rv_record(...)`
- `Estimator::apply_bias_calibration_record()`

### 5) Tight coupling (raw pseudorange + Doppler)

Tight coupling is implemented in `Estimator::update_gps_raw_ref()`:

- pseudorange rows are robustified (heavy-tail protection)
- pseudorange is prevented from directly updating `v` (avoids RV cross-cov noise injection)
- Doppler-only velocity geometry is checked; the weakest direction is projected out when ill-conditioned
- raw GPS updates are recorded as replay events, so delayed / out-of-order tight-coupled measurements use the same rewind/replay path as GNSS PVT, star, and radar

This is the root fix behind the tight-mode tail failures (velocity “noise injection” under degenerate Doppler geometry).

### 6) Bias strategy

Accelerometer bias is handled with a dedicated **pre-calibration** (`BiasCal`) rather than relying on a single EKF to do “trajectory + bias calibration” under blackout-level model mismatch.

The calibrator accumulates a small normal equation over the pre-blackout window and applies a bounded correction at blackout entry.

### 7) Memory budget & history configuration

History is stored as:

- snapshots at a **100 Hz** default rate (`history_snapshot_dt = 0.01`)
- IMU events at the IMU rate (for replay)
- aggregate Φ/Q blocks between snapshots (for RTS smoothing)

For a 10 s fixed-lag window, this is designed to stay within the `tasks.md` **≤ 20 MB** target.

At runtime, you can inspect the current buffer sizes using `Estimator::health()` (includes `history_capacity_bytes`).
