#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
python - <<'PY'
import json
from pathlib import Path
p=Path('outputs/STAGE04C2_R3_GATE.json')
if not p.is_file(): raise SystemExit('HOLD: missing STAGE04C2_R3_GATE.json')
g=json.loads(p.read_text())
if not g.get('gate_pass'): raise SystemExit('HOLD: R3 gate_pass is not true')
print(g['gate'])
PY
TS=$(date +%Y%m%d_%H%M%S)
PKG="ART_MATH_DISCOVERY_STAGE04C2_R3_POSTFREEZE_CASSETTE_UNBLINDING_RESULTS_${TS}.tar.gz"
tar -czf "$PKG" outputs logs README.md PACKAGE_MANIFEST.sha256 run_selftest.sh run_stage04c2_r3.sh scripts tests
sha256sum "$PKG" > "${PKG}.sha256"
sha256sum -c "${PKG}.sha256"
ls -lh "$PKG" "${PKG}.sha256"
