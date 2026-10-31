#!/usr/bin/env bash
set -euo pipefail
P=${1:?project root required}
MODE=${2:-head}
UNITIG_DIR=${3:-}
ROOT=$(cd "$(dirname "$0")" && pwd)
OUT="$ROOT/outputs"
rm -rf "$OUT"; mkdir -p "$OUT" "$ROOT/logs"
ARGS=(--project-root "$P" --outdir "$OUT")
if [[ "$MODE" == "head" || "$MODE" == "local" ]]; then ARGS+=(--network-head); fi
if [[ "$MODE" == "local" ]]; then
  [[ -n "$UNITIG_DIR" ]] || { echo "MODE=local requires unitig dir" >&2; exit 2; }
  ARGS+=(--unitig-dir "$UNITIG_DIR")
fi
python "$ROOT/scripts/stage04b4r1.py" "${ARGS[@]}" 2>&1 | tee "$ROOT/logs/run_stage04b4r1.log"
sha256sum "$OUT"/* > "$OUT/STAGE04B4R1_OUTPUT_SHA256.tsv.tmp"
sed "s#  $OUT/#  outputs/#" "$OUT/STAGE04B4R1_OUTPUT_SHA256.tsv.tmp" > "$OUT/STAGE04B4R1_OUTPUT_SHA256.tsv"
rm -f "$OUT/STAGE04B4R1_OUTPUT_SHA256.tsv.tmp"
echo PASS_STAGE04B4R1_ATOMIC_CLOSEOUT
