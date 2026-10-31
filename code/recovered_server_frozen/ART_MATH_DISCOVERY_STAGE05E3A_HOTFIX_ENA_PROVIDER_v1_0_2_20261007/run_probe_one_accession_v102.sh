#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
  echo "usage: $0 <E3A_ROOT> [accession]" >&2
  exit 2
fi
E3A="$(cd "$1" && pwd)"
ACC="${2:-JN849462.1}"
FIX="$(cd "$(dirname "$0")" && pwd)"
rm -rf /tmp/e3a_ena_probe
mkdir -p /tmp/e3a_ena_probe/out
printf '%s\n' "$ACC" > /tmp/e3a_ena_probe/ids.txt
python "$FIX/scripts/fetch_ena_to_genbank_v102.py" \
  --ids /tmp/e3a_ena_probe/ids.txt \
  --outdir /tmp/e3a_ena_probe/out \
  --manifest /tmp/e3a_ena_probe/manifest.tsv \
  --diagnostics /tmp/e3a_ena_probe/diagnostics.tsv \
  --sleep 0
column -t -s $'\t' /tmp/e3a_ena_probe/diagnostics.tsv || cat /tmp/e3a_ena_probe/diagnostics.tsv
cat /tmp/e3a_ena_probe/manifest.tsv
