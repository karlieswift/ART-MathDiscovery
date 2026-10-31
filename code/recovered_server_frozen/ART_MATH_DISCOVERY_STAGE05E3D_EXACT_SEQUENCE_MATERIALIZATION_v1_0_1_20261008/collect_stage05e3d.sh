#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
test -f outputs/STAGE05E3D_GATE.json || { echo "missing outputs/STAGE05E3D_GATE.json" >&2; exit 2; }
python - <<'PY'
import json
x=json.load(open('outputs/STAGE05E3D_GATE.json'))
if not x.get('gate_pass'): raise SystemExit('E3D gate is not PASS')
PY
TS=$(date +%Y%m%d_%H%M%S)
PKG="ART_STAGE05E3D_EXACT_SEQUENCE_MATERIALIZATION_RESULTS_${TS}.tar.gz"
tar -czf "$PKG" outputs logs config docs scripts README.md run_stage05e3d.sh collect_stage05e3d.sh tests
sha256sum "$PKG" > "$PKG.sha256"
echo "$HERE/$PKG"
cat "$PKG.sha256"
