#!/usr/bin/env bash
set -euo pipefail
python - <<'PY'
import json
from pathlib import Path
p=Path('outputs/STAGE05F2C_BOLTZ_GATE.json')
if not p.exists(): raise SystemExit('missing Boltz gate')
g=json.loads(p.read_text())
if not g.get('gate_pass'): raise SystemExit('Boltz gate not PASS')
print('[PASS] Stage05F2C Boltz gate')
PY
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05F2C_NARROW_BOLTZ2_STRUCTURE_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" outputs logs inputs/boltz2 config/FROZEN_TARGET.json README.md docs/SCIENTIFIC_BOUNDARY.md
sha256sum "$OUT" | tee "$OUT.sha256"
echo "$(pwd)/$OUT"
