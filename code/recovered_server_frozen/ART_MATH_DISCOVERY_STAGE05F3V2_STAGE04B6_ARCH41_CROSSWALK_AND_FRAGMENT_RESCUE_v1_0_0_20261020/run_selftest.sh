#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python scripts/selftest.py
python -m py_compile scripts/stage04b6.py scripts/selftest.py
echo PASS_STAGE04B6_PACKAGE_COMPILE
