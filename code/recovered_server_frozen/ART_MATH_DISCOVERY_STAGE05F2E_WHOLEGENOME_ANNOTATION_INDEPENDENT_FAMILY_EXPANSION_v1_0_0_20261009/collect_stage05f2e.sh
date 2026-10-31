#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
python - <<'PY'
import json
g=json.load(open("outputs/STAGE05F2E_GATE.json"))
assert g["gate_pass"] is True, g
assert g["query_n"]==30, g
assert g["exact_query_recovery_n"]==30, g
assert g["genome_n"]==919, g
print("[PASS] Stage05F2E gate")
PY
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05F2E_WHOLEGENOME_ANNOTATION_INDEPENDENT_FAMILY_EXPANSION_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" \
  outputs \
  raw/STAGE05F2E_QUERY30_VS_WHOLE919_SIXFRAME.m8 \
  raw/STAGE05F2E_RETURNED_TARGET_RVT_ALL.domtblout \
  raw/STAGE05F2E_RETURNED_TARGET_RVT_ALL.tblout \
  logs \
  config \
  docs \
  README.md
sha256sum "$OUT" | tee "$OUT.sha256"
echo "$ROOT/$OUT"
