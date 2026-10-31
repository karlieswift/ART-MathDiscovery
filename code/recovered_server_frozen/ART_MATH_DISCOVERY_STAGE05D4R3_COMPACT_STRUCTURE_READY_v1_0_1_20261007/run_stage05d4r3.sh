#!/usr/bin/env bash
set -euo pipefail
D1="${1:?D1 path required}"
D4="${2:?D4 path required}"
R1="${3:?R1 path required}"
R2="${4:?R2 path required}"
mkdir -p outputs logs
python scripts/build_structure_ready.py --d1 "$D1" --d4 "$D4" --r1 "$R1" --r2 "$R2" --out outputs 2>&1 | tee logs/run_stage05d4r3.log
