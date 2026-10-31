#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
test -f outputs/STAGE05E3D_GATE.json || { echo "missing outputs/STAGE05E3D_GATE.json" >&2; exit 2; }
python - <<'PY'
import json
x=json.load(open('outputs/STAGE05E3D_GATE.json')); assert x.get('gate_pass') is True,x
PY
TS=$(date +%Y%m%d_%H%M%S)
PKG="ART_STAGE05E3D_EXACT_WITNESS_SEQUENCE_RESULTS_${TS}.tar.gz"
tar -czf "$PKG" outputs logs config docs scripts tests README.md run_stage05e3d.sh collect_stage05e3d.sh
sha256sum "$PKG" > "$PKG.sha256"
echo "$HERE/$PKG"; cat "$PKG.sha256"
