#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
runs="${1:-1000}"
output_root="${2:-output_validation}"
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-4}"
cargo test --release --locked
cargo test --release --locked --no-default-features --features flight
cargo build --release --locked --bin harness --bin resourcebench
ASSESSMENT=deadline OUT_DIR="$output_root/pvt" target/release/harness "$runs" --check
ASSESSMENT=causal OUT_DIR="$output_root/tight" TIGHT=1 target/release/harness "$runs" --check
target/release/resourcebench > "$output_root/resources_pvt.txt"
target/release/resourcebench --tight > "$output_root/resources_tight.txt"
printf '%s\n' "Numerical checks finished. ARM target timing and calibrated vehicle evidence remain separate."
