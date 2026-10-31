#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
python tests/test_policy.py
python scripts/validate_registry.py >/tmp/stage04c0_validate.$$ 
grep -q PASS_STAGE04C0_REGISTRY_SCHEMA_AND_HOLD_POLICY /tmp/stage04c0_validate.$$
rm -f /tmp/stage04c0_validate.$$
echo '[PASS] Stage04C0 registry package selftest'
