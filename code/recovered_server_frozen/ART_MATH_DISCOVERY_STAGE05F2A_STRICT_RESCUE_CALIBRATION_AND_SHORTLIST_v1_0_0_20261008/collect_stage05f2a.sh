#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python - <<'PY'
import json
from pathlib import Path
p=Path('outputs/STAGE05F2A_GATE.json')
if not p.exists(): raise SystemExit('[HOLD] missing gate')
g=json.loads(p.read_text())
if not g.get('gate_pass'): raise SystemExit('[HOLD] gate not pass')
print('[PASS] Stage05F2A gate')
PY
stamp=$(date +%Y%m%d_%H%M%S)
out="ART_STAGE05F2A_STRICT_RESCUE_CALIBRATION_AND_SHORTLIST_RESULTS_${stamp}.tar.gz"
tar -czf "$out" outputs logs config docs README.md
sha256sum "$out" | tee "$out.sha256"
printf '%s\n' "$PWD/$out"
