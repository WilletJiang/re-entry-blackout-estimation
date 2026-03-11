## Verification (Reproducible runs)

### 1) Unit / smoke tests

```bash
cargo test
```

Flight-style build check:

```bash
cargo test --no-default-features --features flight
```

### 2) Monte Carlo harness

PVT mode:

```bash
cargo run --release --bin harness -- 500
```

Tight coupling:

```bash
TIGHT=1 cargo run --release --bin harness -- 500
```

Current 500-run release campaign (verified on 2026-03-11):

- Tight: `500/500` pass
- PVT: `~97%` pass (tail dominated by accelerometer-bias RMS and rare NEES spikes)

Current release-build quick regression (verified on 2026-03-11):

- Tight: `100/100` pass
- PVT: `96/100` pass

Commands used for the 2026-03-11 quick regression:

```bash
OUT_DIR=output_20260311_pvt target/release/harness 100
OUT_DIR=output_20260311_tight TIGHT=1 target/release/harness 100
```

These runs write fresh artifacts under `output_20260311_pvt/` and `output_20260311_tight/` so they do not overwrite the final 500-run benchmark outputs under `output/`.

Outputs are written under `output/`:

- `mc_summary_{pvt|tight}.csv`
- `mc_nees_blocks_{pvt|tight}.csv`
- `mc_nees_phase_{pvt|tight}.csv`
- `mc_post_window_{pvt|tight}.csv`
- `mc_*.svg`

### 3) Single-seed reproduction

```bash
cargo run --release --bin repro -- 231
TIGHT=1 cargo run --release --bin repro -- 231
```

Traces (only when needed):

```bash
TRACE=1 cargo run --release --bin repro -- 206
TRACE_GNSS=1 TIGHT=1 cargo run --release --bin repro -- 231
```

### 4) Hot-path timing

```bash
cargo run --release --bin microbench
```

This prints release-build timings for:

- `propagate()`
- `time_update_sqrt()`
- `update_gnss()`
- `update_gps_raw()`

### 5) Report rebuild

```bash
cd docs/report
xelatex main.tex
```

Rebuild the PDF after the final Monte Carlo verification so the report text and figures match the current code and measured results.
