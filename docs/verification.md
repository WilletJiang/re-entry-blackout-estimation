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

Expected (verified on 2025-12-28):

- Tight: `500/500` pass
- PVT: `~97%` pass (tail dominated by accelerometer-bias RMS and rare NEES spikes)

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

