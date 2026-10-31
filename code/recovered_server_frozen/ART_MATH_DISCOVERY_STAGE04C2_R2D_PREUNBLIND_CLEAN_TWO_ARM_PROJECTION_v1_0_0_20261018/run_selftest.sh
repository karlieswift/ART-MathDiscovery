#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
COMMON="$ROOT/docs/HISTORICAL_common.py"
sha256sum "$COMMON" | grep -q '^a59b0e8439747487238b25249116f476baf54e55a5a6681411cb2272763cf255 '
TEST_HIST_COMMON="$COMMON" python "$ROOT/tests/test_adapter.py"
# Run original historical synthetic test against the exact copied common module by staging names.
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/scripts" "$TMP/tests"
cp "$COMMON" "$TMP/scripts/common.py"
cp "$ROOT/docs/HISTORICAL_test_two_arm.py" "$TMP/tests/test_two_arm.py"
python "$TMP/tests/test_two_arm.py"
echo '[PASS] Stage04C2-R2D selftest complete'
