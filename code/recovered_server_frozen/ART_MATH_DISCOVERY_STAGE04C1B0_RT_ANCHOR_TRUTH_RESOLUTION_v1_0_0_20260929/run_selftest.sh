#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
python "$ROOT/tests/test_policy.py"
python -m py_compile "$ROOT/scripts/parse_ddbj_flatfile.py" "$ROOT/scripts/resolve_truth.py"
echo '[PASS] Stage04C1B0 package selftest'
