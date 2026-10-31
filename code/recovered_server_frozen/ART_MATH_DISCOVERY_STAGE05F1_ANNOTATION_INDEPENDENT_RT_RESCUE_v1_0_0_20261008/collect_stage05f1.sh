#!/usr/bin/env bash
set -euo pipefail
python - <<'PY'
import json, pathlib, sys
p=pathlib.Path('outputs/STAGE05F1_GATE.json')
if not p.exists(): raise SystemExit('[HOLD] missing gate')
d=json.loads(p.read_text())
if not d.get('gate_pass'): raise SystemExit('[HOLD] gate false')
print('[PASS] Stage05F1 gate')
PY
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05F1_ANNOTATION_INDEPENDENT_RT_RESCUE_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" outputs logs config docs README.md
sha256sum "$OUT" > "$OUT.sha256"
printf '%s\n' "$PWD/$OUT"
cat "$OUT.sha256"
