#!/usr/bin/env bash
set -euo pipefail
P=${1:?project root required}
THREADS=${2:-32}
WORKERS=${3:-4}
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
mkdir -p logs work
./run_selftest.sh
python scripts/stage04b5.py preflight --project-root "$P" --home "$HERE" \
  2>&1 | tee logs/preflight.log
python scripts/stage04b5.py scan --home "$HERE" --workers "$WORKERS" \
  2>&1 | tee logs/scan.log
STAGE="outputs.stage.$$"
rm -rf "$STAGE"; mkdir -p "$STAGE"
python scripts/stage04b5.py close --home "$HERE" --outdir "$HERE/$STAGE" --threads "$THREADS" \
  2>&1 | tee logs/close.log
rm -rf outputs.previous
if [ -d outputs ]; then mv outputs outputs.previous; fi
mv "$STAGE" outputs
echo PASS_STAGE04B5_ATOMIC_OUTPUT_COMMIT
