#!/usr/bin/env bash
set -euo pipefail
R3=${1:?usage: run_selftest.sh STAGE05D4R3_ROOT}
python tests/selftest.py "$R3"
