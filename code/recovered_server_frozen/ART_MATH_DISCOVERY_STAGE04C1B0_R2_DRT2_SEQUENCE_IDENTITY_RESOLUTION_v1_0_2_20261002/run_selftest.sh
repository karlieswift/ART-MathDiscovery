#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
python -m py_compile "$ROOT/scripts/parse_ddbj_flatfile.py" "$ROOT/scripts/resolve_truth.py"
python "$ROOT/tests/test_policy.py"
python "$ROOT/tests/test_ec78_resolver_logic.py"
python "$ROOT/tests/test_drt2_resolver_logic.py"
echo '[PASS] Stage04C1B0-R2 package selftest'
