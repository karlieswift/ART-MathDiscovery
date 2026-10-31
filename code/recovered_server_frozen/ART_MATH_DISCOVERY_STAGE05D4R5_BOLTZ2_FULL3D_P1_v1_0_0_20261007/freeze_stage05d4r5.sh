#!/usr/bin/env bash
set -euo pipefail
R3=${1:?R3_DIR}; R4=${2:?R4_DIR}
rm -rf frozen inputs
mkdir -p frozen inputs
python scripts/freeze_r5.py --r3 "$R3" --r4 "$R4" --out frozen | tee stage05d4r5_freeze.log
cp -a frozen/p1_yaml inputs/
