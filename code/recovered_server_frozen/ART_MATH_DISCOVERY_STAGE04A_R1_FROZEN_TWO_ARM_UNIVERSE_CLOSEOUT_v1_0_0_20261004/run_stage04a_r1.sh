#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
STAGE04A_ROOT="${1:?usage: ./run_stage04a_r1.sh /path/to/STAGE04A_ROOT}"
SRC="$STAGE04A_ROOT/outputs/STAGE04A_ALL_FROZEN_UNION_REGIONS.tsv"
UNI="$STAGE04A_ROOT/outputs/STAGE04A_UNIVERSE_MANIFEST.tsv"
python "$ROOT/scripts/build_stage04a_r1.py" --source "$SRC" --universe "$UNI" --outdir "$ROOT/outputs"
echo '[PASS] Stage04A-R1 frozen two-arm universe closeout complete'
