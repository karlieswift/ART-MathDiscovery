#!/usr/bin/env bash
set -euo pipefail
if [ $# -ne 1 ]; then echo "usage: $0 FAILED_C1B0_V100_ROOT"; exit 2; fi
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$(readlink -f "$1")"
for spec in "SEN2_14028S CP001363.1" "EC78_ECO7 CP026202.1" "DRT2_KPN AP006725.1"; do
  set -- $spec; cid=$1; acc=$2
  src="$SRC/inputs/genbank/${cid}__${acc}.gb"
  dst="$ROOT/inputs/genbank/${cid}__${acc}.gb"
  test -s "$src" || { echo "[FAIL] missing failed-run flatfile: $src"; exit 2; }
  grep -q '^LOCUS' "$src" || { echo "[FAIL] source is not GenBank flatfile: $src"; exit 2; }
  grep -q "$acc" "$src" || { echo "[FAIL] accession missing from source: $acc"; exit 2; }
  cp -f "$src" "$dst"
  python "$ROOT/scripts/parse_ddbj_flatfile.py" --flatfile "$dst" --out "$ROOT/work_${cid}_features.json"
  echo "[REUSE] $cid $acc"
done
