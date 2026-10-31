#!/usr/bin/env bash
set -euo pipefail
R3=${1:?R3_RESULTS_ROOT}; R5=${2:?R5_ROOT}
ROOT=$(cd "$(dirname "$0")" && pwd); cd "$ROOT"
rm -rf frozen inputs
mkdir -p frozen inputs
python scripts/prepare_p2.py --r3 "$R3" --r5 "$R5" --out frozen
cp -a frozen/p2_yaml inputs/p2_yaml
cat frozen/STAGE05D4R5D_FREEZE.json
