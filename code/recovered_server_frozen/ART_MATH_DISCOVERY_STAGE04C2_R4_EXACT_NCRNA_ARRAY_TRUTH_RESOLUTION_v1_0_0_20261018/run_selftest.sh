#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
python "$HERE/tests/test_logic.py"
echo "[PASS] Stage04C2-R4 selftest complete"
