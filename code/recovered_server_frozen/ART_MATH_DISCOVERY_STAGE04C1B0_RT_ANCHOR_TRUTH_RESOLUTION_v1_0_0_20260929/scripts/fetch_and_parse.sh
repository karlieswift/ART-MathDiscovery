#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
for spec in "SEN2_14028S CP001363.1" "EC78_ECO7 CP026202.1" "DRT2_KPN AP006725.1"; do
  set -- $spec; cid=$1; acc=$2
  url="https://getentry.ddbj.nig.ac.jp/getentry/na/${acc}?format=flatfile&filetype=text"
  out="$ROOT/inputs/genbank/${cid}__${acc}.gb"
  echo "[FETCH] $cid $acc"
  curl -L --fail --retry 5 --retry-all-errors --connect-timeout 30 "$url" -o "$out"
  grep -q '^LOCUS' "$out" || { echo "[FAIL] not flatfile: $out"; exit 2; }
  grep -q "$acc" "$out" || { echo "[FAIL] accession missing: $acc"; exit 2; }
  python "$ROOT/scripts/parse_ddbj_flatfile.py" --flatfile "$out" --out "$ROOT/work_${cid}_features.json"
done
