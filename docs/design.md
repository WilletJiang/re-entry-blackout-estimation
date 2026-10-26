# Estimator design

The core is a 24-dimensional error-state square-root EKF. The nominal state is
`r, v, q, omega, bg, ba, residual_force, residual_torque`. Attitude errors are
left-multiplicative, in the inertial frame: `q_true = Exp(delta_theta) q_est`.
The reference harness preserves the causal IMU-rate assessment over `[160,162] s`.
It separately evaluates the task's permitted retrospective checkpoint estimates,
conditioned only on observations received by the recovery deadline at 162 s.

## Propagation

Translation uses debiased measured specific force. Aerodynamic residual force
must not also be injected into its covariance as though it were an unmeasured
acceleration. Attitude uses the debiased interval-average gyro reading. Its
bias Jacobian is `-R(q) dt`, not a body-frame Jacobian attached to an inertial
attitude error.

Angular rate is estimated by a causal second-order rate/torque observer. With
natural frequency `w = 2 rad/s`, `a = sqrt(2) w dt`, `b = w² dt`, and nominal
rigid-body angular acceleration `alpha`, its update is

```
innovation = gyro - bg - omega_old - alpha * dt/2
omega_new = omega_old + alpha * dt + a * innovation
torque_residual_new = exp(-dt/tau) * torque_residual_old + J * b * innovation
```

Nominal acceleration includes aerodynamic moment, commanded torque, residual
torque, and `-omega × J omega`. This estimates the current rate; it introduces no
fixed output delay. Its response bandwidth is finite. Attitude propagation
continues to use the gyro directly, avoiding dependence of attitude integration
on aerodynamic moment fidelity.

The transition retains rate memory and the coupling to gyro bias and residual
torque. Gyro noise enters attitude, rate, and torque through the **same** random
input, with gains `R dt`, `a I`, and `J b`. The process covariance includes their
cross terms. Its angular 9×9 block is factored without adding arbitrary noise to
zero-variance dimensions. QR propagation preserves the covariance square root.
The local angular model is a linearization, not a guarantee under arbitrary
unmodelled aerodynamic moments; independent physical model validation remains
necessary.

## Observations and bias recovery

The default estimator estimates random-walk biases continuously, using the full
cross covariance in GNSS and star updates. Radar uses the same gain for the state
and Joseph covariance updates. The earlier constant-bias batch calibration is
still available through `set_precalibration(true, ...)`, but is disabled by
default and is not used by the reference campaign.

GNSS PVT and raw pseudorange/Doppler share the timestamped replay mechanism.
There is no special clean GNSS fix or scheduled covariance reset at blackout
exit. Ordinary 10 Hz observations remain subject to dropout and contamination.
The raw path checks its original innovation before residual-dependent covariance
inflation, so a severe outlier cannot make itself appear consistent by inflating
its own measurement covariance. The default sequential mode uses all observable
state directions, without the former recovery-only geometry projection. With a
large uncertain raw-GNSS prior, iterative WLS first handles the nonlinear position
acquisition. Its post-fit residual and its innovation against the prior are both
gated; its full joint position/velocity covariance is assimilated once. The
original raw NIS is retained. Repeated raw rejections also enter acquisition.
Eight internally valid acquisition solutions rejected by the prior trigger RV
covariance widening, as in the PVT recovery path; inconsistent WLS residuals break
the streak. This prevents a rare Gaussian prior-tail draw from rejecting all later
good observations. The next epoch must pass the ordinary gates. A rejected PVT position no longer discards an
independently gated valid velocity measurement.

Eight consecutive rejected PVT epochs trigger covariance widening for a lost
navigation solution. Three isolated outliers were too short a trigger under the
specified 5% contamination. For stars, three **consistent** rejected attitude
innovations trigger reacquisition; unrelated outliers do not. This repairs an
incorrect first star fix that would otherwise lock out subsequent good fixes.
Reacquisition counters and previous rejected star innovations are checkpointed.

## History and smoothing

Every attempted observation, including a rejected one, is retained. A delayed
observation restores the last snapshot at its measurement time, then re-evaluates
all later observations during replay. Discarding rejected observations would
make gate decisions and reacquisition counters depend on arrival order.

Pruning preserves an anchor at or just before the full ten-second cutoff and all
IMU events needed after that anchor. RTS smoothing traverses backwards with one
state/covariance pair; it does not allocate three arrays proportional to history
length. `visit_smoothed` streams these checkpoints newest-first and uses a
Cholesky solve rather than an explicit predicted-covariance inverse. The deadline
assessment calls it once at 162 s, retains exactly the 100 Hz grid in `[160,162]`,
and never revises those results using later arrivals. The causal 400 Hz results
remain separately available. This is a first-order EKF/RTS approximation, not a
continuous-time guarantee between checkpoints. Propagation after warmup and the smoother are checked for zero heap
allocation. Delayed replay and raw GNSS still allocate temporary workspaces;
`resourcebench` includes them in the peak measurement.

## Reference-model assumptions

The reference scenario is synthetic. PVT position/velocity standard deviations
are 1.9 m and 0.08 m/s per axis; raw pseudorange/range-rate standard deviations
are 1.5 m and 0.02 m/s. These modes do not have identical GNSS information quality.
The simulator uses idealized inertial satellite orbits and omits receiver clock
states, atmospheric propagation, and sensor scale-factor errors. Physical-model
validation and target execution measurements are separate from the numerical
screen; see [verification.md](verification.md).
