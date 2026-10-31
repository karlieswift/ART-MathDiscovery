#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
python - <<'PY'
import json
from pathlib import Path
p=Path('outputs/STAGE05E3E_GATE.json')
if not p.exists(): raise SystemExit('missing outputs/STAGE05E3E_GATE.json')
g=json.loads(p.read_text())
if not g.get('gate_pass'): raise SystemExit(f'gate not PASS: {g.get("gate")}')
print('[PASS] Stage05E3E gate')
PY
TS=$(date +%Y%m%d_%H%M%S)
PKG="ART_STAGE05E3E_PRIORITY_FAMILY_BIOLOGICAL_CLOSURE_RESULTS_${TS}.tar.gz"
tar -czf "$PKG" outputs logs config docs scripts tests README.md run_stage05e3e.sh collect_stage05e3e.sh
sha256sum "$PKG" > "${PKG}.sha256"
echo "$PWD/$PKG"
cat "${PKG}.sha256"
