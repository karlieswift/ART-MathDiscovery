#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 1 ]; then echo "Usage: $0 /path/to/C0_ROOT" >&2; exit 2; fi
SRC=$(cd "$1" && pwd); ROOT=$(cd "$(dirname "$0")/.." && pwd)
REG="$SRC/registry/RT_SYSTEM_CONTROL_REGISTRY.tsv"; MAP="$ROOT/config/ACCESSION_MAP.tsv"
test -s "$REG" || { echo "[FAIL] registry missing: $REG" >&2; exit 3; }
mkdir -p "$ROOT/inputs/sequences" "$ROOT/inputs/flatfiles" "$ROOT/outputs/record_gates" "$ROOT/logs"
rm -f "$ROOT/outputs/record_gates"/*.json "$ROOT/inputs/sequences"/* "$ROOT/inputs/flatfiles"/* 2>/dev/null || true
: > "$ROOT/logs/fetch_ddbj_mapped.log"
python "$ROOT/scripts/validate_map.py" --registry "$REG" --map "$MAP" --out "$ROOT/outputs/ACCESSION_MAP_GATE.json" | tee "$ROOT/logs/00_accession_map.log"
tail -n +2 "$MAP" | while IFS=$'\t' read -r cid regacc fetchacc relation explen note; do
  [ -n "$cid" ] || continue
  fa="$ROOT/inputs/sequences/${cid}__${fetchacc}.fasta"; gb="$ROOT/inputs/flatfiles/${cid}__${fetchacc}.flatfile.txt"
  echo "[FETCH] $cid registry=$regacc source=$fetchacc relation=$relation" | tee -a "$ROOT/logs/fetch_ddbj_mapped.log"
  tmpfa="${fa}.part"; tmpgb="${gb}.part"; rm -f "$tmpfa" "$tmpgb"
  curl -fL --retry 8 --retry-all-errors --retry-delay 2 --connect-timeout 30 --max-time 1800 \
    "https://getentry.ddbj.nig.ac.jp/getentry/na/${fetchacc}?format=fasta&filetype=text" -o "$tmpfa"
  curl -fL --retry 8 --retry-all-errors --retry-delay 2 --connect-timeout 30 --max-time 1800 \
    "https://getentry.ddbj.nig.ac.jp/getentry/na/${fetchacc}?format=flatfile&filetype=text" -o "$tmpgb"
  mv "$tmpfa" "$fa"; mv "$tmpgb" "$gb"
  python "$ROOT/scripts/validate_record.py" --control-id "$cid" --registry-accession "$regacc" --fetch-accession "$fetchacc" \
    --relation "$relation" --expected-length-bp "${explen:-0}" --fasta "$fa" --flatfile "$gb" --min-bp 1000000 \
    --out "$ROOT/outputs/record_gates/${cid}.json" 2>&1 | tee -a "$ROOT/logs/fetch_ddbj_mapped.log"
done
python "$ROOT/scripts/finalize.py" "$SRC"
