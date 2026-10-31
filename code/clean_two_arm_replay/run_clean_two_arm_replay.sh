#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <MarsHill_Stage02B_root> <Stage03A_results_tar.gz>" >&2
  exit 2
fi
ROOT="$(cd "$(dirname "$0")" && pwd)"
MARS="$(readlink -f "$1")"
TAR="$(readlink -f "$2")"
mkdir -p "$ROOT/outputs" "$ROOT/logs"
echo '[0/4] frozen-input preflight'
python "$ROOT/scripts/00_preflight.py" --root "$ROOT" --marshill-root "$MARS" --stage03a-tar "$TAR" 2>&1 | tee "$ROOT/logs/00_preflight.log"
echo '[1/4] construct ACNV+recurrence-only boards from frozen pre-unblinding scores'
python "$ROOT/scripts/01_build_boards.py" --root "$ROOT" --marshill-root "$MARS" --stage03a-tar "$TAR" 2>&1 | tee "$ROOT/logs/01_build_boards.log"
echo '[2/4] post-freeze truth replay only after four boards are frozen'
python "$ROOT/scripts/02_postfreeze_unblind.py" --root "$ROOT" --marshill-root "$MARS" --stage03a-tar "$TAR" 2>&1 | tee "$ROOT/logs/02_postfreeze_unblind.log"
echo '[3/4] manifest'
python "$ROOT/scripts/03_manifest.py" --root "$ROOT" 2>&1 | tee "$ROOT/logs/03_manifest.log"
echo '[PASS] clean two-arm replay complete'
