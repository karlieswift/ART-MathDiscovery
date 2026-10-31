#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 1 ]; then echo "Usage: $0 /path/to/STAGE04C0_R2_ROOT" >&2; exit 2; fi
R2=$(cd "$1" && pwd); ROOT=$(cd "$(dirname "$0")/.." && pwd)
MAP="$ROOT/config/ACCESSION_MAP.tsv"
mkdir -p "$ROOT/inputs/sequences" "$ROOT/inputs/flatfiles" "$ROOT/outputs/record_gates" "$ROOT/logs"
rm -f "$ROOT/outputs/record_gates"/*.json "$ROOT/inputs/sequences"/* "$ROOT/inputs/flatfiles"/* 2>/dev/null || true
: > "$ROOT/logs/01_reuse_validate.log"
tail -n +2 "$MAP" | while IFS=$'\t' read -r cid regacc fetchacc relation explen note; do
  [ -n "$cid" ] || continue
  srcfa="$R2/inputs/sequences/${cid}__${fetchacc}.fasta"
  srcgb="$R2/inputs/flatfiles/${cid}__${fetchacc}.flatfile.txt"
  test -s "$srcfa" || { echo "[FAIL] missing R2 FASTA: $srcfa" >&2; exit 3; }
  test -s "$srcgb" || { echo "[FAIL] missing R2 flatfile: $srcgb" >&2; exit 4; }
  dstfa="$ROOT/inputs/sequences/${cid}__${fetchacc}.fasta"
  dstgb="$ROOT/inputs/flatfiles/${cid}__${fetchacc}.flatfile.txt"
  cp --reflink=auto --preserve=timestamps "$srcfa" "$dstfa"
  cp --reflink=auto --preserve=timestamps "$srcgb" "$dstgb"
  echo "[REUSE] $cid $fetchacc expected_length=$explen" | tee -a "$ROOT/logs/01_reuse_validate.log"
  python "$ROOT/scripts/validate_record.py" --control-id "$cid" --registry-accession "$regacc" --fetch-accession "$fetchacc" \
    --relation "$relation" --expected-length-bp "$explen" --fasta "$dstfa" --flatfile "$dstgb" \
    --out "$ROOT/outputs/record_gates/${cid}.json" 2>&1 | tee -a "$ROOT/logs/01_reuse_validate.log"
done
python "$ROOT/scripts/finalize.py" "$R2" | tee "$ROOT/logs/02_finalize.log"
