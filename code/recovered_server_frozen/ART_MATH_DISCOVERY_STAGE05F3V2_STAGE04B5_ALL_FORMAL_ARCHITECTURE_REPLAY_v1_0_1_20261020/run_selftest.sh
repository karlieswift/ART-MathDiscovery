#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python scripts/selftest.py
python -m py_compile scripts/stage04b5.py scripts/selftest.py
echo PASS_STAGE04B5_PACKAGE_COMPILE
