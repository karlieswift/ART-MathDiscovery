#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
./run_selftest.sh
mkdir -p raw/genbank outputs logs
python scripts/fetch_genbank.py \
  --ids frozen/STAGE05E3A_SEQUENCE_IDS.txt \
  --outdir raw/genbank \
  --manifest raw/STAGE05E3A_GENBANK_RETRIEVAL_MANIFEST.tsv \
  2>&1 | tee logs/fetch_genbank.log
python scripts/analyze_official_context.py \
  --frozen frozen \
  --gbdir raw/genbank \
  --outdir outputs \
  2>&1 | tee logs/analyze_official_context.log
