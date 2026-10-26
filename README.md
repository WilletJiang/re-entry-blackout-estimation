# Entry Blackout Estimation

[中文](README_CN.md) · [Task specification](tasks.md)

A Rust navigation-estimation study using a 24-state square-root error-state EKF,
a causal angular-rate observer, timestamped measurement replay, and fixed-lag RTS
smoothing. It supports GNSS PVT and raw pseudorange/Doppler updates.

## Reference scenario

The command-line tools default to `SCENARIO=reference-entry`, an explicitly
synthetic descending trajectory. The preview aerodynamic and RCS tables are not a
calibrated full-flight model. `SCENARIO=legacy` retains the original trajectory;
`SimConfig::default()` preserves that initial condition for library callers.

IMU/GNSS/star/radar rates are 400/10/2/20 Hz. Stars and radar have 120/70 ms delays.
GNSS is absent during `[120,160) s`; stars are unavailable at measurement times
`[130,150] s`. Radar is valid below 25 km. The default stress case includes 5%
heavy-tailed contamination and 20% dropout. Truth noise and sensor random streams
are independent of estimator tuning; paired modes share truth and IMU streams.

## Build and run

```bash
cargo test --release --locked
cargo test --release --locked --no-default-features --features flight
cargo build --release --locked --bin harness --bin repro --bin resourcebench
ASSESSMENT=causal TIGHT=1 OUT_DIR=output_tight target/release/harness 1000 --check
ASSESSMENT=deadline OUT_DIR=output_pvt target/release/harness 1000 --check
target/release/repro 0
TIGHT=1 target/release/repro 0
target/release/resourcebench --tight
```

`--check` returns a nonzero status if the selected numerical assessment fails.
`ASSESSMENT=causal` (default) checks 801 instantaneous outputs in `[160,162] s`.
`ASSESSMENT=deadline` checks 201 retrospective checkpoints using only observations
received by 162 s. Both results are always written separately; retrospective
estimates must not be presented as instantaneous outputs.

`SEED_OFFSET` selects the first batch seed. `scripts/validate-reference.sh`
runs the test suites, both assessments, and the resource probes. Generated CSVs,
plots, and logs belong in the ignored `output*` directories.

## Integration

[Design](docs/design.md), [interfaces](docs/interface.md), and
[verification](docs/verification.md) describe state/error conventions, noise
units, recovery gates, replay requirements, and assessment definitions.

`resourcebench` measures estimator allocations plus a reserved thread stack and
reports wall-clock and Linux thread-CPU timing. RTS batch time is reported
separately. Use `--check-time` on the ARM reference target to check the 5 ms wall
budget; host measurements and this synthetic scenario are not flight certification.
