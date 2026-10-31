#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
python "$ROOT/tests/test_resolver.py"
echo '[PASS] Stage04C2 resolver selftest'
