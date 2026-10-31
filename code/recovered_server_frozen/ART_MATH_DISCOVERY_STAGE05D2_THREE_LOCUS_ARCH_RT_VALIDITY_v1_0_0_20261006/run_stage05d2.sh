#!/usr/bin/env bash
set -euo pipefail
D1="${1:?usage: ./run_stage05d2.sh /path/to/Stage05D1 [threads]}"
THREADS="${2:-8}"
mkdir -p outputs work/ncbi
python scripts/audit_three_loci.py "$D1" config/FROZEN_LOCI.tsv outputs

# Current official annotation retrieval. Failure is non-fatal but explicitly recorded.
: > outputs/OFFICIAL_FETCH_STATUS.tsv
printf 'accession\tstatus\tpath\n' >> outputs/OFFICIAL_FETCH_STATUS.tsv
for ACC in NC_031007.1 NC_048173.1 NC_062581.1; do
  OUT="work/ncbi/${ACC}.gb"
  if curl -fL --retry 1 --retry-delay 1 --connect-timeout 10 --max-time 20 \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&rettype=gbwithparts&retmode=text" \
    -o "$OUT" && test -s "$OUT"; then
      printf '%s\tPASS\t%s\n' "$ACC" "$OUT" >> outputs/OFFICIAL_FETCH_STATUS.tsv
  else
      rm -f "$OUT"
      printf '%s\tFAIL\tNA\n' "$ACC" >> outputs/OFFICIAL_FETCH_STATUS.tsv
  fi
done
python scripts/official_annotation_audit.py "$D1" config/FROZEN_LOCI.tsv work/ncbi outputs
python scripts/make_architecture_svg.py outputs/THREE_LOCUS_ARCHITECTURE_RT_AUDIT.tsv outputs/THREE_LOCUS_ARCHITECTURE.svg
python scripts/write_gate.py "$D1" outputs
