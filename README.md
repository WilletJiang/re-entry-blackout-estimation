# Entry Blackout Estimation

Chinese documentation: [README_CN.md](README_CN.md)

## Overview

This repository addresses the state estimation problem defined by task `SE 2317` for a reusable orbital stage during atmospheric entry. The task is defined by quantitative and operational criteria. During plasma blackout, GNSS becomes unavailable for forty seconds, inertial drift grows under thermal stress, and the vehicle must recover a controllable navigation solution within two seconds after signal reacquisition. Success depends on nominal accuracy, delay handling, observability management, covariance credibility, and bounded runtime.

The repository is a constrained estimator study with explicit dynamics, sensor timing, anomaly assumptions, and acceptance criteria. The vehicle follows six degree of freedom rigid body dynamics with aerodynamic control and RCS actuation across a strongly nonlinear entry regime from Mach two to Mach twelve. The IMU runs at four hundred hertz. GNSS runs at ten hertz and is fully unavailable from `t0 + 120 s` to `t0 + 160 s`. The star tracker runs at two hertz, arrives with a fixed delay of one hundred and twenty milliseconds, and is unavailable from `t0 + 130 s` to `t0 + 150 s`. The radar altimeter runs at twenty hertz, is delayed by seventy milliseconds, and is valid only below twenty five kilometres. Measurements may drop out, may contain heavy tailed outliers, and may arrive out of order.

## Estimator Design

The current implementation uses a 24 state error state square root EKF. The state carries position, velocity, attitude, angular rate, gyroscope bias, accelerometer bias, and residual aerodynamic force and torque terms. Covariance propagation remains in square root form through QR based updates so that numerical conditioning stays controlled during long blackout propagation and during aggressive measurement reacquisition. Translational propagation is strapdown IMU based. Truth propagation and filter side nominal linearization call the same nominal aerodynamic and RCS composition helper, which keeps the deterministic part of the dynamics aligned across simulation and estimation.

Delay handling is a first class requirement in this implementation. Delayed or out of order measurements are applied through a fixed lag history buffer that performs rewind, update, and deterministic replay at the measurement timestamp. This mechanism is used for GNSS PVT, raw GNSS, star tracker, and radar updates. The result is a time aligned state and covariance history even when the arrival sequence differs from the physical sequence. For the blackout interval itself, the repository also includes fixed lag RTS smoothing so that information recovered after signal return can improve the estimate over the trailing segment of the unobservable window.

The GNSS update path supports two observation models. The baseline mode uses position and velocity fixes. The tighter mode uses pseudorange and Doppler directly. In that tighter mode the code explicitly checks Doppler geometry degeneracy and removes the weakest velocity direction when the observation geometry becomes ill conditioned. This choice prevents noise injection during reacquisition. Robustness is further shaped by chi square gating, Student t style covariance inflation, optional innovation adaptive estimation, and a dedicated pre blackout accelerometer bias calibration path.

## Acceptance and Verification

The acceptance logic follows the mission window defined in `tasks.md`. Accuracy is evaluated only in the two second recovery interval immediately after blackout exit. The thresholds are forty metres for position at three sigma, zero point five metres per second for velocity at three sigma, zero point fifteen degrees for attitude at three sigma, zero point zero zero five degrees per second for gyroscope bias at three sigma, and one hundred and twenty micro g for accelerometer bias at three sigma. Accuracy alone is not sufficient. The repository also tracks NEES and NIS so that statistical consistency is assessed together with recovery performance.

The current verification record is the one documented on March 11, 2026 in `verification.md`. In release mode, the 500 run Monte Carlo campaign passes all 500 runs in tight mode and 485 of 500 runs in PVT mode. An additional 100 run regression on the same date passes 100 of 100 runs in tight mode and 96 of 100 runs in PVT mode. The remaining PVT misses are dominated by accelerometer bias tail behaviour and occasional NEES excursions, which is consistent with the lower information content of position and velocity only updates during the recovery segment.

## Build and Reproduction

The repository targets stable Rust. A standard release build is sufficient for the benchmark harness and the diagnostics workflow.

```bash
cargo build --release
```

The default feature set includes simulation and diagnostics. A lean flight style build is available when only the core library is required.

```bash
cargo test --no-default-features --features flight
```

The main benchmark entry point is the Monte Carlo harness. Running without environment switches uses the PVT path. Setting `TIGHT=1` enables raw pseudorange and Doppler updates.

```bash
cargo run --release --bin harness -- 500
TIGHT=1 cargo run --release --bin harness -- 500
```

Single seed reproduction is provided for failure analysis and regression tracking.

```bash
cargo run --release --bin repro -- 231
TIGHT=1 cargo run --release --bin repro -- 231
```

When diagnosis is necessary, trace modes can be enabled for chronological inspection of propagation and measurement handling.

```bash
TRACE=1 cargo run --release --bin repro -- 206
TRACE_GNSS=1 TIGHT=1 cargo run --release --bin repro -- 231
```

The harness writes summary csv files, phasewise consistency diagnostics, blockwise consistency diagnostics, post blackout window statistics, and svg figures into the configured output directory.

## Repository Layout

The estimator implementation is organised under `src` with separate subtrees for propagation, measurement updates, history management, calibration, robustness support, and smoothing. The simulation and evaluation code is also under `src` and contains truth generation, sensor synthesis, Monte Carlo execution, and metric aggregation. `physics.rs` contains the shared nominal dynamics helper. `data.rs` loads the vehicle property, aerodynamic coefficient, and RCS preview tables from `data`. `types.rs` defines the public state and noise structures. `plot.rs` emits csv and svg artifacts when diagnostics are enabled.

Three supplementary texts support the repository. `design.md` explains the estimator decisions and numerical safeguards. `interface.md` specifies timestamp conventions and flight style calling patterns. `verification.md` records the reproducible validation procedure and the current benchmark status. The typeset technical report in the report subdirectory under `docs` provides the longer mathematical and experimental narrative.
