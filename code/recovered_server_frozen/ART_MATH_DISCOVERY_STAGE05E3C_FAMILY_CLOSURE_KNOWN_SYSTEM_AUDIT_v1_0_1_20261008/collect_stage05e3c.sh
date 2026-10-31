#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
[ -f outputs/STAGE05E3C_GATE.json ] || { echo "missing outputs/STAGE05E3C_GATE.json" >&2; exit 2; }
python - <<'PY'
import json
x=json.load(open('outputs/STAGE05E3C_GATE.json'))
assert x.get('gate_pass') is True
PY
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT="ART_STAGE05E3C_FAMILY_CLOSURE_RESULTS_${STAMP}.tar.gz"
tar -czf "$OUT" POLICY.json README_RUN.txt docs refs outputs scripts run_stage05e3c.sh collect_stage05e3c.sh tests/run_selftest.sh
sha256sum "$OUT" | tee "$OUT.sha256"
echo "$HERE/$OUT"
