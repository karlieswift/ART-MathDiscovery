#!/usr/bin/env bash
set -euo pipefail
if [ $# -ne 3 ]; then echo "usage: $0 FAILED_C1B0_V100_ROOT C0_ROOT R3_ROOT"; exit 2; fi
ROOT="$(cd "$(dirname "$0")" && pwd)"
FAILED="$(readlink -f "$1")"; C0="$(readlink -f "$2")"; R3="$(readlink -f "$3")"
"$ROOT/run_selftest.sh"
"$ROOT/scripts/reuse_and_parse.sh" "$FAILED" 2>&1 | tee "$ROOT/logs/01_reuse_parse.log"
python "$ROOT/scripts/resolve_truth.py" --root "$ROOT" --c0 "$C0" --r3 "$R3" 2>&1 | tee "$ROOT/logs/02_resolve.log"
sha256sum "$ROOT"/inputs/genbank/*.gb "$ROOT"/outputs/* > "$ROOT/outputs/C1B0_R1_SHA256SUMS.txt"
echo '[PASS] Stage04C1B0-R1 truth resolution complete'
