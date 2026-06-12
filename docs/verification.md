# Verification Procedure

This document defines the reproducible verification procedure for the current repository state. It is deliberately procedural. A validation claim is only meaningful if another reader can rerun the same commands, recover the same artifact set, and interpret the results against the same acceptance basis in `tasks.md`.

## Scope

Verification covers four layers. Unit and smoke tests confirm local invariants. Flight build checks confirm that the core library builds without diagnostics features. Monte Carlo campaigns evaluate the acceptance window and consistency metrics. Microbenchmark timing measures whether the current implementation remains comfortably below the per cycle compute budget. The report under `docs/report` should be rebuilt after the final benchmark run so that figures and narrative remain synchronized with the measured results.

## Required Baseline Checks

Run the default test suite first.

```bash
cargo test
```

Then confirm that the estimator also builds and tests in the lean flight configuration.

```bash
cargo test --no-default-features --features flight
```

These two commands establish that the codebase is internally consistent before expensive campaign runs begin.

The diagnostics smoke regressions replay the full mission timeline, including the blackout interval and the recovery window. They are covered by the default test suite, but the release smoke target is the practical command when the goal is to validate the blackout specific regressions introduced in the current repository state.

```bash
cargo test --release --test smoke
```

## Monte Carlo Campaign

The main acceptance evidence is produced by the benchmark harness. Running the harness without environment switches exercises the GNSS PVT path. Setting `TIGHT=1` enables the raw pseudorange and Doppler path. The default harness configuration uses a ten second fixed lag history horizon so that the validation path matches the task level memory and smoothing requirement.

```bash
cargo run --release --bin harness -- 500
TIGHT=1 cargo run --release --bin harness -- 500
```

For regression work that must not overwrite the canonical benchmark directory, direct the outputs to dedicated folders.

```bash
OUT_DIR=output_20260311_pvt target/release/harness 100
OUT_DIR=output_20260311_tight TIGHT=1 target/release/harness 100
```

The first pair of commands is the full campaign. The second pair is the quick regression pattern used in the March 11, 2026 validation record.

## Historical Record

The repository carries a historical release record dated March 11, 2026. In the 500 run campaign used at that time, tight mode passes 500 of 500 runs. PVT mode passes 485 of 500 runs. In the additional 100 run regression on the same date, tight mode passes 100 of 100 runs and PVT mode passes 96 of 100 runs.

Those numbers remain useful for comparison, but they predate the current acceptance alignment and summary schema. The harness now records recovery window RMS diagnostics, recovery window norm equivalent three sigma maxima, smoothed recovery window RMS diagnostics, and peak `Estimator::health()` history memory counters. A fresh full campaign is required before claiming a current validated record under the aligned acceptance logic.

On the repository state dated June 13, 2026, a fresh 10 run release regression passes in both modes. Running

```bash
OUT_DIR=output_current_pvt target/release/harness 10
OUT_DIR=output_current_tight TIGHT=1 target/release/harness 10
```

produces 10 of 10 passes in each mode. In that regression, PVT keeps `sigma3_vel_max_mps` between 0.415554 and 0.415679, with `sigma3_pos_max_m` near 9.8723 meters. Tight coupling keeps `sigma3_vel_max_mps` between 0.103921 and 0.104390, with `sigma3_pos_max_m` again near 9.8723 meters. These runs establish that the present code line clears the earlier PVT velocity acceptance miss and the earlier tight coupled forward filter divergence. They remain a quick regression rather than a replacement for the historical long campaign.

## Single Seed Reproduction

When a campaign identifies a failure or a suspicious tail case, reproduce the seed directly.

```bash
cargo run --release --bin repro -- 231
TIGHT=1 cargo run --release --bin repro -- 231
```

If chronological diagnosis is required, enable tracing.

```bash
TRACE=1 cargo run --release --bin repro -- 206
TRACE_GNSS=1 TIGHT=1 cargo run --release --bin repro -- 231
```

These traces are intended for targeted investigation. They should not be part of the routine acceptance workflow.

When a summary metric must be localized in time, the blackout trace recorder and the release smoke regressions are the preferred path. They carry position error, velocity error, position three sigma, velocity three sigma, and RV NEES on a shared time grid, using the same norm equivalent three sigma convention that drives pass and fail in the summary csv.

## Timing Check

The timing envelope of the critical kernels can be sampled with the microbenchmark binary.

```bash
cargo run --release --bin microbench
```

This reports release build timings for `propagate`, `time_update_sqrt`, `update_gnss`, and `update_gps_raw`. The benchmark is a development machine measurement. Target processor certification still requires hardware specific validation. It remains useful for detecting regressions that would threaten the five millisecond single cycle budget.

## Artifact Set

A successful harness run writes a stable artifact family into the chosen output directory. The key files are `mc_summary_{pvt|tight}.csv`, `mc_nees_blocks_{pvt|tight}.csv`, `mc_nees_phase_{pvt|tight}.csv`, `mc_post_window_{pvt|tight}.csv`, and the generated `mc_*.svg` figures. The summary file captures per run pass and metric values, including recovery window three sigma maxima and peak history memory counters. The phase and block diagnostics localize consistency problems by mission segment and state group.

## Report Rebuild

After the final benchmark run, rebuild the technical report from the report directory.

```bash
cd docs/report
xelatex main.tex
```

A second compilation pass may be needed if references or figure placements change. The report should only be considered current when its tables, figures, and narrative match the latest validated artifact set.

Verification is complete only when the commands above succeed, the generated artifacts exist, and the resulting pass rates and consistency summaries are judged against the quantitative requirements in `tasks.md`.
