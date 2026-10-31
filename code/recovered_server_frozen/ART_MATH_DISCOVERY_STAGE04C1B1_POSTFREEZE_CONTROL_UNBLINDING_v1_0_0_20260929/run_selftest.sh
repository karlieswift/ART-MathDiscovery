#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
python "$ROOT/tests/test_policy.py"
python -m py_compile "$ROOT/scripts/unblind.py"
echo '[PASS] Stage04C1B1 package selftest'
