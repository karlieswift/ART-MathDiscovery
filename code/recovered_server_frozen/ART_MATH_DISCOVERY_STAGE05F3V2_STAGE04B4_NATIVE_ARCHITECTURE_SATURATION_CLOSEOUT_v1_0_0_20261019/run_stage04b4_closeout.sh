#!/usr/bin/env bash
set -euo pipefail
P=${1:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
THREADS=${2:-32}
mkdir -p logs
./run_selftest.sh | tee logs/selftest.log
STAGE="outputs.stage.$$"
rm -rf "$STAGE"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then echo "[HOLD] primary closeout stopped; inspect $STAGE and logs/run_stage04b4_console.log"; fi; exit $rc' EXIT
python scripts/stage04b4_closeout.py \
  --project-root "$P" \
  --threads "$THREADS" \
  --outdir "$STAGE" \
  2>&1 | tee logs/run_stage04b4_console.log
rm -rf outputs
mv "$STAGE" outputs
trap - EXIT
echo "PASS_STAGE04B4_ATOMIC_OUTPUT_COMMIT"
