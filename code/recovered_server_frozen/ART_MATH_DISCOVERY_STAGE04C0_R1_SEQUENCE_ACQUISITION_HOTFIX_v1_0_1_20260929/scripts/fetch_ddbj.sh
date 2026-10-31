#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 /path/to/ART_MATH_DISCOVERY_STAGE04C0_RT_SYSTEM_CONTROL_REGISTRY_v1_0_0_20260929" >&2
  exit 2
fi
SRC=$(cd "$1" && pwd)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
REG="$SRC/registry/RT_SYSTEM_CONTROL_REGISTRY.tsv"
test -s "$REG" || { echo "[FAIL] registry missing: $REG" >&2; exit 3; }
mkdir -p "$ROOT/inputs/sequences" "$ROOT/inputs/flatfiles" "$ROOT/outputs/record_gates" "$ROOT/logs"
: > "$ROOT/logs/fetch_ddbj.log"
python - "$REG" <<'PY' > "$ROOT/outputs/accessions.tsv"
import csv,sys
for r in csv.DictReader(open(sys.argv[1]),delimiter='\t'):
    if r.get('whole_sequence_accession'):
        print(r['control_id']+'\t'+r['whole_sequence_accession'])
PY
while IFS=$'\t' read -r cid acc; do
  [ -n "$acc" ] || continue
  fa="$ROOT/inputs/sequences/${cid}__${acc}.fasta"
  gb="$ROOT/inputs/flatfiles/${cid}__${acc}.flatfile.txt"
  echo "[FETCH] $cid $acc" | tee -a "$ROOT/logs/fetch_ddbj.log"
  curl -fL --retry 5 --retry-all-errors --connect-timeout 30 \
    "https://getentry.ddbj.nig.ac.jp/getentry/na/${acc}?format=fasta&filetype=text" \
    -o "$fa"
  curl -fL --retry 5 --retry-all-errors --connect-timeout 30 \
    "https://getentry.ddbj.nig.ac.jp/getentry/na/${acc}?format=flatfile&filetype=text" \
    -o "$gb"
  python "$ROOT/scripts/validate_record.py" \
    --control-id "$cid" --accession "$acc" --fasta "$fa" --flatfile "$gb" \
    --min-bp 1000000 --out "$ROOT/outputs/record_gates/${cid}.json" \
    2>&1 | tee -a "$ROOT/logs/fetch_ddbj.log"
done < "$ROOT/outputs/accessions.tsv"
python "$ROOT/scripts/finalize.py" "$SRC"
