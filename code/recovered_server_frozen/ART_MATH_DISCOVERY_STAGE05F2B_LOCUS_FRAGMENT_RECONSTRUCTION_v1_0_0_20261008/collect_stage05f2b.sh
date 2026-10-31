#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
python - <<'PY'
import json
from pathlib import Path
p=Path('outputs/STAGE05F2B_GATE.json')
if not p.exists(): raise SystemExit('missing gate')
g=json.loads(p.read_text())
if not g.get('gate_pass'): raise SystemExit('gate not PASS')
print('[PASS] Stage05F2B gate')
PY
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT="ART_STAGE05F2B_LOCUS_FRAGMENT_RECONSTRUCTION_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" README.md config docs outputs logs
sha256sum "$OUT" | tee "$OUT.sha256"
readlink -f "$OUT"
