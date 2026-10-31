#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
python "$HERE/scripts/stage04b8.py" selftest
python -m py_compile "$HERE/scripts/stage04b8.py"
echo PASS_STAGE04B8_PACKAGE_COMPILE
