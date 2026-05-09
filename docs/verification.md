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

## Monte Carlo Campaign

The main acceptance evidence is produced by the benchmark harness. Running the harness without environment switches exercises the GNSS PVT path. Setting `TIGHT=1` enables the raw pseudorange and Doppler path.

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

## Current Verified Record

The repository currently carries a validated release record dated March 11, 2026. In the 500 run campaign, tight mode passes 500 of 500 runs. PVT mode passes 485 of 500 runs. In the additional 100 run regression on the same date, tight mode passes 100 of 100 runs and PVT mode passes 96 of 100 runs.

These numbers indicate that tight mode provides stronger tail control during the post blackout recovery window. The residual PVT failures are dominated by accelerometer bias RMS tails and occasional NEES excursions, which is consistent with the lower observability provided by position and velocity only updates.

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

## Timing Check

The timing envelope of the critical kernels can be sampled with the microbenchmark binary.

```bash
cargo run --release --bin microbench
```

This reports release build timings for `propagate`, `time_update_sqrt`, `update_gnss`, and `update_gps_raw`. The benchmark is a development machine measurement. Target processor certification still requires hardware specific validation. It remains useful for detecting regressions that would threaten the five millisecond single cycle budget.

## Artifact Set

A successful harness run writes a stable artifact family into the chosen output directory. The key files are `mc_summary_{pvt|tight}.csv`, `mc_nees_blocks_{pvt|tight}.csv`, `mc_nees_phase_{pvt|tight}.csv`, `mc_post_window_{pvt|tight}.csv`, and the generated `mc_*.svg` figures. The summary file captures per run pass and metric values. The phase and block diagnostics localize consistency problems by mission segment and state group.

## Report Rebuild

After the final benchmark run, rebuild the technical report from the report directory.

```bash
cd docs/report
xelatex main.tex
```

A second compilation pass may be needed if references or figure placements change. The report should only be considered current when its tables, figures, and narrative match the latest validated artifact set.

Verification is complete only when the commands above succeed, the generated artifacts exist, and the resulting pass rates and consistency summaries are judged against the quantitative requirements in `tasks.md`.
