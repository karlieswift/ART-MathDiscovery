#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$ROOT/inputs/sequences" "$ROOT/inputs/genbank" "$ROOT/logs"
python - "$ROOT" <<'PY'
import csv,sys
from pathlib import Path
root=Path(sys.argv[1])
rows=list(csv.DictReader((root/'registry/RT_SYSTEM_CONTROL_REGISTRY.tsv').open(),delimiter='\t'))
for r in rows:
    if r['whole_sequence_accession']:
        print(r['control_id']+'\t'+r['whole_sequence_accession'])
PY
while IFS=$'\t' read -r cid acc; do
  [ -n "$acc" ] || continue
  fa="$ROOT/inputs/sequences/${cid}__${acc}.fasta"
  gb="$ROOT/inputs/genbank/${cid}__${acc}.gb"
  echo "[FETCH] $cid $acc"
  curl -fL --retry 3 --retry-delay 2 \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text" \
    -o "$fa"
  curl -fL --retry 3 --retry-delay 2 \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=gbwithparts&retmode=text" \
    -o "$gb"
  test -s "$fa" && test -s "$gb"
done < <(python - "$ROOT" <<'PY'
import csv,sys
from pathlib import Path
root=Path(sys.argv[1])
for r in csv.DictReader((root/'registry/RT_SYSTEM_CONTROL_REGISTRY.tsv').open(),delimiter='\t'):
    if r['whole_sequence_accession']: print(r['control_id']+'\t'+r['whole_sequence_accession'])
PY
)
python "$ROOT/scripts/freeze_sequences.py"
