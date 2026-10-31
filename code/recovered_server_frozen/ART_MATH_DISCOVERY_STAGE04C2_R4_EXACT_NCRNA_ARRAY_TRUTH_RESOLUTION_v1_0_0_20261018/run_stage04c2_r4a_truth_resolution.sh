#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <UG_SOURCE_ROOT> <R3_ROOT> [resolver args...]" >&2
  echo "Resolver args: --truth-tsv FILE | --cm-main FILE --cm-alt FILE | --unit-fasta FILE" >&2
  exit 2
fi
UG=$1; R3=$2; shift 2
HERE=$(cd "$(dirname "$0")" && pwd)
OUT="$HERE/outputs/r4a_truth"
mkdir -p "$OUT" "$HERE/logs"
python "$HERE/scripts/00_preflight_source.py" --ug-root "$UG" --r3-root "$R3" --out "$OUT" | tee "$HERE/logs/00_preflight.log"
python "$HERE/scripts/01_inventory_author_evidence.py" --ug-root "$UG" --out "$OUT" | tee "$HERE/logs/01_inventory.log"
set +e
python "$HERE/scripts/02_resolve_exact_truth.py" --ug-root "$UG" --r3-root "$R3" --out "$OUT" "$@" | tee "$HERE/logs/02_truth_resolve.log"
RC=${PIPESTATUS[0]}
set -e
python "$HERE/scripts/04_finalize.py" --outputs "$HERE/outputs" | tee "$HERE/logs/04_finalize_r4a.log"
exit "$RC"
