#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
test -s outputs/STAGE05F2D_GATE.json
python - <<'PY'
import json
x=json.load(open('outputs/STAGE05F2D_GATE.json'))
assert x['gate_pass'] and x['gate']=='PASS_STAGE05F2D_ANNOTATION_INDEPENDENT_FAMILY_EXPANSION',x
print('[PASS] Stage05F2D gate')
PY
sha256sum raw/STAGE05F2D_QUERY30_VS_F1_UNIQUE.m8 outputs/STAGE05F2D_* > outputs/STAGE05F2D_SHA256.tsv
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05F2D_ANNOTATION_INDEPENDENT_FAMILY_EXPANSION_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" outputs raw logs config docs README.md
sha256sum "$OUT" | tee "$OUT.sha256"
echo "$ROOT/$OUT"
