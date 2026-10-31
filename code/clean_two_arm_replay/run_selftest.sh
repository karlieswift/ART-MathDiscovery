#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
python "$ROOT/tests/test_two_arm.py"
python -m py_compile "$ROOT"/scripts/*.py "$ROOT"/tests/*.py
echo '[PASS] clean two-arm replay package selftest'
