#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <Stage04C2_R2C_root> <historical_clean_two_arm_root>" >&2
  exit 2
fi
ROOT="$(cd "$(dirname "$0")" && pwd)"
R2C="$(readlink -f "$1")"
HIST="$(readlink -f "$2")"
mkdir -p "$ROOT/outputs" "$ROOT/logs"
rm -rf "$ROOT/outputs/contigs"
rm -f "$ROOT/outputs"/*.tsv "$ROOT/outputs"/*.json

echo '[0/3] preflight: frozen R2C scores + historical exact clean rule'
python "$ROOT/scripts/00_preflight.py" --root "$ROOT" --r2c "$R2C" --historical-clean "$HIST" 2>&1 | tee "$ROOT/logs/00_preflight.log"

echo '[1/3] build 15 clean ACNV+recurrence boards from frozen score tables; truth not read'
python "$ROOT/scripts/01_build_clean_boards.py" --root "$ROOT" --r2c "$R2C" --historical-clean "$HIST" 2>&1 | tee "$ROOT/logs/01_build_clean_boards.log"

echo '[2/3] freeze clean board + score-only comparison; still no truth read'
python "$ROOT/scripts/02_finalize.py" --root "$ROOT" --r2c "$R2C" 2>&1 | tee "$ROOT/logs/02_finalize.log"

echo '[PASS] Stage04C2-R2D pre-unblinding clean board frozen'
