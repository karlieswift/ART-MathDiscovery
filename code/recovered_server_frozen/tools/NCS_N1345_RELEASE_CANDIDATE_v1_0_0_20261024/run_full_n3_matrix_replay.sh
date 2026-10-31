#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
OUT="${1:-$ROOT/full_matrix_replay}"
mkdir -p "$OUT"
python "$ROOT/scripts_optional/replay_sensitivity_from_full_matrices.py" \
 --matrix-dir "$ROOT/data/N3_original_matrices" --out "$OUT"
cmp "$OUT/N3_SENSITIVITY_6_NEW_MATRICES.tsv" "$ROOT/data/N3_SENSITIVITY_6_NEW_MATRICES.tsv"
cmp "$OUT/N3_REGION_RANKS_6_NEW_MATRICES.tsv.gz" "$ROOT/data/N3_REGION_RANKS_6_NEW_MATRICES.tsv.gz"
echo 'PASS_FULL_SIX_MATRIX_RECOMBINATION_TWO_OUTPUTS_BITWISE_MATCH'
