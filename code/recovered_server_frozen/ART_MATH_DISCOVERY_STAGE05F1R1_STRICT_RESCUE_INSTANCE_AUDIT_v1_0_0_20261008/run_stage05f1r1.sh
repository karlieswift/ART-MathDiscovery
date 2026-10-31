#!/usr/bin/env bash
set -euo pipefail
F1=${1:?usage: run_stage05f1r1.sh <F1 result directory or result tar> [project_root]}
P=${2:-.}
mkdir -p outputs work logs
python scripts/stage05f1r1.py "$F1" "$P" --outdir outputs 2>&1 | tee logs/stage05f1r1.log
python - <<'PYCHECK'
import json
x=json.load(open('outputs/STAGE05F1R1_GATE.json'))
assert x['gate_pass'] is True
print('[PASS] Stage05F1R1 strict instance deconvolution')
PYCHECK
