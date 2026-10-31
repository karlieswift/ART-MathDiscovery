#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
python - <<'PY'
import json
g=json.load(open("outputs/STAGE05F2E_GATE.json"))
assert g["gate_pass"] is True, g
assert g["query_n"]==30, g
assert g["coordinate_replay_query_n"]==30, g
assert g["direct_wholegenome_hash_recovery_n"]==30, g
assert g["mmseqs_final_exact_recovery_n"]==30, g
assert g["exact_query_recovery_n"]==30, g
assert g["genome_n"]==919, g
print("[PASS] Stage05F2E v1.0.1 gate")
PY
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05F2E_WHOLEGENOME_ANNOTATION_INDEPENDENT_FAMILY_EXPANSION_v1_0_1_RESULTS_${STAMP}.tar.gz"
FILES=(
  outputs
  raw/STAGE05F2E_QUERY30_VS_WHOLE919_SIXFRAME.m8
  raw/STAGE05F2E_QUERY30_VS_WHOLE919_SIXFRAME_PRIMARY_MASK0.m8
  raw/STAGE05F2E_RETURNED_TARGET_RVT_ALL.domtblout
  raw/STAGE05F2E_RETURNED_TARGET_RVT_ALL.tblout
  logs
  config
  docs
  README.md
)
if [ -s raw/STAGE05F2E_QUERY_MISSING_PRIMARY_EXHAUSTIVE.m8 ]; then
  FILES+=(raw/STAGE05F2E_QUERY_MISSING_PRIMARY_EXHAUSTIVE.m8)
fi
tar -czf "$OUT" "${FILES[@]}"
sha256sum "$OUT" | tee "$OUT.sha256"
echo "$ROOT/$OUT"
