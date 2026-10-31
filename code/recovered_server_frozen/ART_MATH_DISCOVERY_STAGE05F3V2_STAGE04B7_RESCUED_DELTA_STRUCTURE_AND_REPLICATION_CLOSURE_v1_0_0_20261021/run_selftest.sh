#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
python "$HERE/scripts/stage04b7.py" selftest
python -m py_compile "$HERE/scripts/stage04b7.py"
echo PASS_STAGE04B7_PACKAGE_COMPILE
