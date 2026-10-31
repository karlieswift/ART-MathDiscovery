#!/usr/bin/env bash
set -euo pipefail
if [ $# -ne 2 ]; then echo "usage: $0 C0_ROOT R3_ROOT"; exit 2; fi
ROOT="$(cd "$(dirname "$0")" && pwd)"
C0=$(readlink -f "$1"); R3=$(readlink -f "$2")
"$ROOT/run_selftest.sh"
"$ROOT/scripts/fetch_and_parse.sh" 2>&1 | tee "$ROOT/logs/01_fetch_parse.log"
python "$ROOT/scripts/resolve_truth.py" --root "$ROOT" --c0 "$C0" --r3 "$R3" 2>&1 | tee "$ROOT/logs/02_resolve.log"
sha256sum "$ROOT"/inputs/genbank/*.gb "$ROOT"/outputs/* > "$ROOT/outputs/C1B0_SHA256SUMS.txt"
echo '[PASS] Stage04C1B0 truth resolution complete'
