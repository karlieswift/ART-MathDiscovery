#!/usr/bin/env bash
set -euo pipefail
D1="${1:?usage: ./run_stage05d2.sh /path/to/Stage05D1 [threads]}"
THREADS="${2:-8}"
mkdir -p outputs work/ncbi
python scripts/audit_three_loci.py "$D1" config/FROZEN_LOCI.tsv outputs

valid_gb() {
  local f="$1"
  test -s "$f" && grep -q '^LOCUS[[:space:]]' "$f" && grep -q '^ORIGIN' "$f" && grep -q '^//$' "$f"
}

fetch_gb() {
  local acc="$1" out="$2" tmp="${out}.tmp"
  rm -f "$tmp" "$out"

  # Primary: EDirect if available in the environment.
  if command -v efetch >/dev/null 2>&1; then
    if efetch -db nuccore -id "$acc" -format gb > "$tmp" 2>/dev/null && valid_gb "$tmp"; then
      mv "$tmp" "$out"
      echo "EDIRECT"
      return 0
    fi
  fi
  rm -f "$tmp"

  # Fallback: documented E-utilities GenBank flat-file endpoint.
  if curl -fsSL --retry 2 --retry-delay 2 --connect-timeout 15 --max-time 60 \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=gb&retmode=text" \
      -o "$tmp" && valid_gb "$tmp"; then
    mv "$tmp" "$out"
    echo "EUTILS_GB"
    return 0
  fi
  rm -f "$tmp"
  return 1
}

: > outputs/OFFICIAL_FETCH_STATUS.tsv
printf 'accession\tstatus\tmethod\tpath\n' >> outputs/OFFICIAL_FETCH_STATUS.tsv
for ACC in NC_031007.1 NC_048173.1 NC_062581.1; do
  OUT="work/ncbi/${ACC}.gb"
  if METHOD="$(fetch_gb "$ACC" "$OUT")"; then
    printf '%s\tPASS\t%s\t%s\n' "$ACC" "$METHOD" "$OUT" >> outputs/OFFICIAL_FETCH_STATUS.tsv
  else
    printf '%s\tFAIL_INVALID_OR_UNAVAILABLE\tNA\tNA\n' "$ACC" >> outputs/OFFICIAL_FETCH_STATUS.tsv
  fi
done

# Parse failures are recorded per accession and are non-fatal.
python scripts/official_annotation_audit.py "$D1" config/FROZEN_LOCI.tsv work/ncbi outputs
python scripts/make_architecture_svg.py outputs/THREE_LOCUS_ARCHITECTURE_RT_AUDIT.tsv outputs/THREE_LOCUS_ARCHITECTURE.svg
python scripts/write_gate.py "$D1" outputs
