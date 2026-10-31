#!/usr/bin/env bash
set -euo pipefail
D1="${1:?usage: ./run_stage05d4.sh /path/to/Stage05D1}"
mkdir -p outputs logs
python scripts/build_stage05d4_universe.py "$D1" --out outputs | tee logs/run_stage05d4.log

echo
echo '=== GATE ==='
cat outputs/GATE_STAGE05D4.json

echo
echo '=== OUTPUTS ==='
ls -lh outputs
