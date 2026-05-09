# GN&C Task Specification: SE 2317 Entry Blackout State Estimation

## 1. Scope

This specification defines the state estimation task for a reusable orbital class stage during atmospheric entry. The estimator shall maintain a navigation solution through a plasma induced GNSS blackout and shall recover a guidance usable state estimate within two seconds after blackout exit. Compliance is defined by quantitative accuracy, consistency, robustness, runtime, and memory requirements. Unless stated otherwise, the requirements in this document are normative.

## 2. Mission Scenario

The vehicle is a six degree of freedom rigid body that performs a belly first blunt body entry. Attitude control is provided by RCS actuation together with movable aerodynamic surfaces. The aerodynamic regime is strongly nonlinear over approximately Mach two through Mach twelve and is accompanied by high dynamic pressure, strong thermal loading, IMU drift growth, and scale factor nonlinearity.

The operational objective is precise. The vehicle shall enter a controllable guidance corridor no later than two seconds after the end of the blackout interval. State error alone is insufficient evidence of success. The reported covariance shall remain statistically credible over the same recovery segment.

## 3. Sensor Suite and Timing Assumptions

The estimator shall support the asynchronous sensor suite described below.

| Sensor | Rate | Observation content | Timing and availability assumptions |
| --- | --- | --- | --- |
| IMU | 400 Hz | Three axis gyroscope and accelerometer | Thermal drift shall be modelled as random walk |
| GNSS | 10 Hz | Position and velocity solution | Complete signal loss from `t0 + 120 s` through `t0 + 160 s` |
| Star tracker | 2 Hz | Quaternion attitude | Fixed delay of `120 ms`; unavailable from `t0 + 130 s` through `t0 + 150 s` |
| Radar altimeter | 20 Hz | Terrain normal altitude | Fixed delay of `70 ms`; valid only below `25 km` |

Measurement anomalies are part of the task assumptions. Observations may arrive out of order. Heavy tailed contamination shall occur with probability five percent, represented by a Gaussian mixture with ninety five percent nominal samples and five percent Cauchy outliers. Random dropout shall also be considered during robustness evaluation.

## 4. Estimation Objective

The estimator shall propagate the vehicle state continuously through the blackout interval and shall assimilate asynchronous measurements before blackout and after reacquisition. Delayed or out of order measurements shall be applied at their measurement timestamp. Receipt time shall not be used as the update time. The estimator shall support a fixed lag smoothing mechanism over a horizon no shorter than ten seconds if smoothing is required to satisfy the recovery objective.

The minimum state content shall cover position, velocity, attitude, angular rate, gyroscope bias, and accelerometer bias. Additional residual states may be introduced if they improve model mismatch absorption without violating the runtime and memory limits.

## 5. Functional Requirements

The estimator shall satisfy the following functional requirements.

1. It shall propagate inertial state at the IMU rate and maintain numerical stability over the full mission segment.
2. It shall process GNSS, star tracker, and radar measurements despite asynchronous arrival and fixed sensor delays.
3. It shall remain operational through complete GNSS loss over the blackout interval.
4. It shall reject or down weight non nominal measurements so that sparse outliers do not corrupt either the state or the covariance.
5. It shall expose sufficient diagnostic information to evaluate innovation consistency and estimation consistency over Monte Carlo campaigns.
6. It shall support replay safe handling of any externally imposed state or covariance intervention that must survive delayed measurement rewinding.

## 6. Reference State and Model Formulation

This section is informative. It defines the reference formulation used by the current repository realization.

The augmented state is written as

$$
x = \big[r,\ v,\ q,\ \omega,\ b_g,\ b_a,\ f_{\mathrm{res}},\ \tau_{\mathrm{res}}\big]^\top
$$

where $r\in\mathbb R^3$ is position, $v\in\mathbb R^3$ is velocity, $q\in\mathbb H$ is the unit quaternion attitude, $\omega\in\mathbb R^3$ is angular rate, $b_g$ and $b_a$ are gyroscope and accelerometer biases, and $f_{\mathrm{res}}$ and $\tau_{\mathrm{res}}$ are residual aerodynamic force and torque states used to absorb model mismatch during blackout and reacquisition.

A representative continuous time model is

$$
\begin{aligned}
\dot r &= v \\
\dot v &= R(q)\big(a_\text{meas}-b_a-n_a\big) + g(r) \\
\dot q &= \tfrac{1}{2}\,\Omega(\omega)\,q \\
\dot \omega &= J^{-1}\!\left(\tau_\text{cmd} + \tau_\text{aero}(r,v) + \tau_{\mathrm{res}} - \omega\times J\omega\right) \\
\dot b_g &= w_g,\quad \dot b_a = w_a \\
\dot f_{\mathrm{res}} &= -\tau_f^{-1} f_{\mathrm{res}} + w_{f,\mathrm{res}},\quad
\dot \tau_{\mathrm{res}} = -\tau_\tau^{-1} \tau_{\mathrm{res}} + w_{\tau,\mathrm{res}}
\end{aligned}
$$

where $R(q)$ is the body to inertial direction cosine matrix, $\Omega(\omega)$ is the quaternion rate matrix, and the $w$ terms are zero mean process noises.

A representative measurement model is

$$
\begin{aligned}
z_\mathrm{GNSS}(t_k) &= h_\mathrm{GNSS}(x(t_k)) + \nu_g \\
z_\mathrm{raw}(t_k) &= h_\mathrm{raw}(x(t_k)) + \nu_\mathrm{raw} \\
z_\star(t_k-\delta t_\star) &= h_\star(x(t_k-\delta t_\star)) + \nu_\star \\
z_h(t_k-\delta t_h) &= h_h(x(t_k-\delta t_h)) + \nu_h
\end{aligned}
$$

Any delayed star tracker, radar, or raw GNSS measurement should be interpreted at its physical measurement time. In the current repository this is realized through rewind, update, and deterministic replay inside a fixed lag history buffer.

## 7. Quantitative Acceptance Criteria

Let `t_exit = t0 + 160 s` denote blackout exit. The primary acceptance window is the interval from `t_exit` through `t_exit + 2 s`. Unless otherwise stated, the quantitative targets below apply to that interval after any permitted fixed lag smoothing has been applied.

### 7.1 Accuracy and Bias Recovery

| Quantity | Requirement |
| --- | --- |
| Position error | `3σ ≤ 40 m` |
| Velocity error | `3σ ≤ 0.5 m s^-1` |
| Attitude error | `3σ ≤ 0.15 deg` equivalent cone angle |
| Angular rate error | `3σ ≤ 0.05 deg s^-1` |
| Gyroscope bias norm | `3σ ≤ 0.005 deg s^-1` |
| Accelerometer bias norm | `3σ ≤ 120 µg` |

### 7.2 Consistency

Over a Monte Carlo campaign of at least one hundred runs, NEES and NIS ninety five percent confidence interval coverage shall be no lower than ninety percent. A solution that satisfies the accuracy targets while systematically understating covariance shall not be accepted.

### 7.3 Robustness

The estimator shall not diverge under a stress campaign that combines five percent outlier contamination with twenty percent random dropout. Under the same campaign, covariance underestimation shall not exceed thirty percent in the aggregate consistency assessment.

### 7.4 Runtime and Memory

The estimator shall support a one hundred hertz main processing loop with a single cycle compute budget no greater than five milliseconds on an ARM A78AE class single core reference target. Runtime memory attributable to the estimator, including the fixed lag history needed for delayed measurement handling and smoothing, shall not exceed twenty megabytes. The fixed lag horizon shall support at least ten seconds of history.

## 8. Repository Realization Note

This section is informative and describes the current repository realization of the task. The repository implements a 24 state error state square root EKF that augments the nominal translational, rotational, and bias states with residual aerodynamic force and torque states. Delayed sensors are handled through a fixed lag rewind, update, and replay buffer, and raw GNSS tight coupled updates participate in the same event replay path. Truth propagation and filter side nominal linearization call the same aerodynamic and RCS composition helper so that both paths share the same deterministic nominal dynamics.
