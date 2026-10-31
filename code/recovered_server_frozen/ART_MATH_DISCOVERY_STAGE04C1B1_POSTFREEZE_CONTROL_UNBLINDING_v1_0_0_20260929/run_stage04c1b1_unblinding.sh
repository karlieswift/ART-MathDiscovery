#!/usr/bin/env bash
set -euo pipefail
if [ $# -ne 2 ];then echo "usage: $0 C1_SCORE_FREEZE_ROOT C1B0_TRUTH_ROOT";exit 2;fi
ROOT="$(cd "$(dirname "$0")" && pwd)"
"$ROOT/run_selftest.sh"
python "$ROOT/scripts/unblind.py" --root "$ROOT" --c1 "$(readlink -f "$1")" --c1b0 "$(readlink -f "$2")" 2>&1 | tee "$ROOT/logs/01_unblind.log"
sha256sum "$ROOT"/outputs/* > "$ROOT/outputs/C1B1_SHA256SUMS.txt"
echo '[PASS] Stage04C1B1 post-freeze unblinding complete'
