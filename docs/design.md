# Design Note: Entry Blackout Estimation

This note records the estimator design currently used to satisfy the acceptance basis defined in `tasks.md`. The estimator must cross a forty second GNSS blackout, recover within the two second post blackout window, preserve covariance credibility, and remain inside the runtime and memory envelope of the task.

## Design Objective

The design is shaped by four constraints. First, blackout removes the most informative translational observation for a prolonged interval. Second, star tracker and radar measurements are delayed and may arrive out of order. Third, entry dynamics are strongly nonlinear and imperfectly modelled, so a numerically fragile covariance form is unacceptable. Fourth, a filter with optimistic covariance creates downstream risk, so NEES and NIS are part of the design target and the validation basis.

## State Form and Propagation

The estimator uses a 24 state error state square root EKF. The nominal state contains position, velocity, attitude, angular rate, gyroscope bias, accelerometer bias, residual aerodynamic force, and residual aerodynamic torque. Attitude is corrected multiplicatively, which keeps the quaternion normalized and localizes covariance to a three component attitude error. Covariance is propagated in square root form through QR based time updates so that positive definiteness and conditioning remain controlled during long inertial propagation.

Translational state propagation follows a strapdown formulation driven by debiased IMU specific force and body rate. At the same time, truth propagation and filter side nominal linearization use the same aerodynamic and RCS composition helper from `physics.rs`. This arrangement gives the simulator and the filter a shared deterministic dynamics model without double counting measured acceleration inside the state integrator.

## Uncertainty Model

The process model is phase aware. The translational noise block and the accelerometer bias random walk can be scaled independently before blackout, during blackout, and after blackout through the `q_scale_rv_*` and `q_scale_ba_*` parameters. Residual aerodynamic force and torque are represented as first order Gauss Markov states with optional blackout specific excitation levels. Additional attitude and translational growth is introduced during the star tracker saturation interval so that covariance reflects the temporary loss of attitude information.

This structure keeps covariance realistic across regime changes. A uniform process noise model would understate uncertainty during blackout and overstate it during nominal observation periods.

## Measurement Update Policy

The update architecture is intentionally conservative. GNSS PVT updates operate on the position and velocity subspace only, which prevents spurious cross covariance from pulling attitude and bias states in poorly supported directions. Star tracker updates remain attitude only. Radar updates remain scalar altitude updates. Raw GNSS updates fuse pseudorange and Doppler directly when tight coupling is enabled.

All updates are gated by NIS based chi square tests. The covariance update uses Joseph form together with explicit positive definiteness enforcement. When numerical roundoff or aggressive tuning threatens symmetry or definiteness, the code falls back through Cholesky repair, diagonal jitter, and spectral regularization.

## Temporal Alignment and History Semantics

Delayed data are handled with a fixed lag history buffer. Every measurement is interpreted at its measurement timestamp. Arrival time is not used as the update time. When a delayed or out of order measurement is received, the estimator restores the latest snapshot not later than that timestamp, replays the recorded event stream to the measurement time, applies the update there, and then replays the remaining events forward to the present state.

This contract yields two important invariants. The first is temporal consistency: state and covariance are updated at the time where the measurement model is valid. The second is replay determinism: any state or covariance modification that must survive a rewind must either be reproduced by replayed events or be checkpointed in a snapshot at the same timestamp. The replay safe interfaces `inflate_rv_record` and `apply_bias_calibration_record` exist for exactly this reason.

## Tight Coupling Strategy

The tight coupled path is implemented in `update_gps_raw_ref`. Pseudorange rows are robustified against heavy tailed residuals. Pseudorange is also prevented from directly updating velocity, which reduces noise injection through position and velocity cross covariance. Doppler rows carry the main velocity information, but that information becomes unreliable when the satellite geometry loses rank in the weakest direction. The implementation therefore checks the Doppler geometry and removes the weakest velocity direction when conditioning deteriorates.

This conditioning step removes the rare and severe recovery phase failures previously observed in the tight coupled tail.

## Bias Management

Accelerometer bias is the most delicate hidden state in this task because blackout removes the measurements that best constrain translational drift. The repository therefore separates bias preparation from the main filter. A dedicated pre blackout calibration stage accumulates a small normal equation over the observation rich interval before blackout and applies a bounded correction at blackout entry. The EKF then carries the corrected bias state through the unobservable segment. Bias identification during blackout is intentionally limited because the available information is weak.

## Runtime and Memory Envelope

History snapshots are stored at a default cadence of one hundred hertz. IMU events are stored at IMU rate for replay, and the propagation records needed for RTS smoothing are retained between snapshots. With a ten second fixed lag horizon, this design is intended to remain within the twenty megabyte task limit. Runtime health can be inspected through `Estimator::health()`, which exposes recent GNSS gate status together with history usage statistics.

In summary, the estimator combines four mechanisms that address distinct blackout recovery risks. The square root covariance form protects numerical stability, the history buffer preserves time consistency, tight coupled raw GNSS improves recovery observability, and pre blackout bias calibration reduces hidden drift before the least observable phase.
