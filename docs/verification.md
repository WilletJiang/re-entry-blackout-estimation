# Reproducible verification

Use `Cargo.lock` and Rust 1.90.0 or a compatible stable release.
The authoritative task limits remain in [tasks.md](../tasks.md).

```bash
cargo test --release --locked
cargo test --release --locked --no-default-features --features flight
cargo build --release --locked --bin harness --bin repro --bin resourcebench
ASSESSMENT=deadline OUT_DIR=output_pvt RAYON_NUM_THREADS=4 target/release/harness 1000 --check
ASSESSMENT=causal OUT_DIR=output_tight RAYON_NUM_THREADS=4 TIGHT=1 target/release/harness 1000 --check
target/release/resourcebench
target/release/resourcebench --tight
```

`--check` returns a nonzero exit status when the selected statistical screen fails.
`ASSESSMENT=causal` (default) selects instantaneous 400 Hz outputs;
`ASSESSMENT=deadline` explicitly selects task §7's permitted retrospective result.
Without `--check`, the harness writes results without enforcing the statistical screen.
`SEED_OFFSET=100` starts a batch at seed 100. A single failure can be reproduced
with `target/release/repro <seed>` and the same `TIGHT`/`SCENARIO` selection.
The harness defaults to `SCENARIO=reference-entry`. `SCENARIO=legacy` retains the
original ascending trajectory for diagnosis; it is not entry-validation evidence.

## Statistical screen

Causal runs cover the complete 801-sample recovery grid at 400 Hz. The deadline
assessment covers every 100 Hz checkpoint from 160 through 162 s (201 samples),
using only data already received at 162 s. The full ten-second history remains
available. Smoothing at 162 s does not mean those values were available at their
historical timestamps; it also does not prove a bound between checkpoints.
Numerical validity is sticky over the mission. All six norm-equivalent bounds,
`3 sqrt(trace(P_block))`, must stay below their limits in every run. No component,
including angular rate, is omitted.

At each timestamp, the campaign also computes `3 sqrt(mean(||error||²))` across
independent seeds, including any systematic error. The maximum of this empirical
quantity across the window must meet each of the six limits. This is stronger
than checking a time-averaged RMS or inspecting covariance alone.

Both assessments use the **joint 18-dimensional** covariance of position,
velocity, left attitude error, angular rate, and both biases. All samples of the selected recovery grid
are included, irrespective of GNSS acceptance. Its central 95% chi-square
interval coverage must reach 90%; `mean(NEES)/18` must not exceed 1.30. Residual
force and torque are auxiliary model states and are not included in this minimum
navigation-state statistic.

GNSS, star, and radar NIS are recorded **before** robustification/gating, including
rejected observations, with exact precomputed central 95% chi-square intervals
for degrees of freedom 1–64. Each sensor class must reach 90% coverage. Outlier
samples remain in the denominator. The reference test requires at least 100
unique seeds. These time samples are correlated; their count is not an effective
number of independent trials and should not be used to invent a binomial
confidence interval.

## Output interpretation

- `campaign.txt`: causal joint statistical screen, empirical maxima, sensor coverages,
  and marginal NEES diagnostics.
- `recovery_ensemble.csv`: causal recovery timestamps and each state-block bound.
- `deadline/{campaign.txt,recovery_ensemble.csv,per_run.csv}`: the separate
  retrospective screen, 201-point ensemble curve, and each seed's covariance maxima.
  Per-run files also export unrounded NEES sums/hits and all NIS counts/hits, so
  untouched seed subsets can be audited without selection or rounding artifacts.
- `per_run.csv`: the corresponding per-seed causal covariance screen.
- `configuration.txt`: scenario, source revision, seed range, rates and stress
  probabilities. Run from a clean checkout to associate results with that revision.
- `mc_summary_{mode}.csv`: per-run maxima and older diagnostics. Its `passed`
  column is only the six-block covariance screen.
- `mc_trajectory_{mode}.csv`: truth altitude/Mach/domain coverage. Paired modes
  use identical physical truth and IMU random streams.
- Older RV-only NEES columns/plots retain their historical accepted-GNSS sampling
  convention. They do not replace the joint recovery NEES in `campaign.txt`.
- Smoothed RMS diagnostics are optional (`enable_smoothing_diagnostics` in the
  library). They are `NaN` when not requested; that is not a successful smoothing
  result. They are separate from the always-captured deadline checkpoint assessment.

## Resources

`resourcebench` runs 2,000 single-thread 100 Hz cycles with four IMU propagations
per cycle, GNSS updates, delayed stars/radar, a full reference aerodynamic table,
and ten seconds of history. The allocator measures live/peak allocation payload,
including nested raw-satellite buffers and replay workspaces. The result adds the
worker's entire reserved 512 KiB stack and enforces **20,000,000 bytes**. It reports
RTS time and checks that RTS performs no heap allocations. This excludes allocator
bookkeeping, executable/shared-library pages, and unrelated application memory;
the target integrator must account for their platform overhead.

Wall-clock cycle p99 and maximum include host scheduling interference. Linux
probes additionally report per-thread CPU p99/maximum and counts exceeding 5 ms.
CPU time cannot replace the wall-clock deadline. `--check-time` makes a measured
wall-clock cycle over 5 ms fail the probe; use it on the specified target. Full
10-second RTS and the two-second recovery pass are timed separately; neither is
hidden inside the causal cycle figure. They are synchronous and require their own
target scheduling budget. The retrospective statistical analysis alone is not
evidence that its batch computation fits a five-millisecond flight cycle.
The five-millisecond requirement is specifically for ARM A78AE-class hardware:
an x86 host measurement cannot certify it. Run the same probe and an integration
trace on the target, including realistic scheduling and I/O, before claiming
hardware compliance.
