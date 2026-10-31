#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd);cd "$ROOT"
python - <<'PY'
import json
for p in ['outputs/STAGE05F3A0_SOURCE_FREEZE_GATE.json','outputs/STAGE05F3A_GATE.json']:
 x=json.load(open(p));assert x['gate_pass'],x
print('[PASS] Stage05F3A gates')
PY
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05F3A_LARGE_UNIVERSE_FROZEN_QUERY_EXPANSION_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" outputs raw logs config docs README.md
sha256sum "$OUT" | tee "$OUT.sha256"
echo "$ROOT/$OUT"
