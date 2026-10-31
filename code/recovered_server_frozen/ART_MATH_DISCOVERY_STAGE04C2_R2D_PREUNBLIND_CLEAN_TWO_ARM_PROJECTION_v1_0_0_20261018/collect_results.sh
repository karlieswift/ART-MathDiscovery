#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
test -f outputs/STAGE04C2_R2D_GATE.json || { echo 'HOLD: no final gate'; exit 2; }
python - <<'PY'
import json
x=json.load(open('outputs/STAGE04C2_R2D_GATE.json'))
assert x.get('gate_pass') is True
assert x.get('truth_loaded') is False
assert x.get('cassette_overlap_computed') is False
print(x['gate'])
PY
TS=$(date +%Y%m%d_%H%M%S)
PKG="ART_MATH_DISCOVERY_STAGE04C2_R2D_PREUNBLIND_CLEAN_TWO_ARM_RESULTS_${TS}.tar.gz"
tar -czf "$PKG" outputs logs README.md scripts tests docs run_selftest.sh run_stage04c2_r2d.sh collect_results.sh PACKAGE_MANIFEST.sha256
sha256sum "$PKG" > "$PKG.sha256"
sha256sum -c "$PKG.sha256"
ls -lh "$PKG" "$PKG.sha256"
