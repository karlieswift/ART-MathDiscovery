#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 1 ]; then echo "Usage: $0 <R2D_ROOT>" >&2; exit 2; fi
R2D=$1
HERE=$(cd "$(dirname "$0")" && pwd)
R4A="$HERE/outputs/r4a_truth"
OUT="$HERE/outputs/r4b_eval"
mkdir -p "$OUT" "$HERE/logs"
python "$HERE/scripts/03_evaluate_frozen_board.py" --r4a-out "$R4A" --r2d-root "$R2D" --out "$OUT" | tee "$HERE/logs/03_exact_array_eval.log"
python "$HERE/scripts/04_finalize.py" --outputs "$HERE/outputs" | tee "$HERE/logs/04_finalize_r4b.log"
