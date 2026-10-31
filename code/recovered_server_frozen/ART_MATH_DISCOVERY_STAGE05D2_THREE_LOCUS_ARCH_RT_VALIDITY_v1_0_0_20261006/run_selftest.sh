#!/usr/bin/env bash
set -euo pipefail
D1="${1:?usage: ./run_selftest.sh /path/to/Stage05D1}"
req=(
  config/NOVELTY_FIRST_999_CANDIDATES.tsv
  config/NOVELTY_FIRST_CANDIDATE_ORF_LINKS.tsv
  config/NOVELTY_FIRST_NONCODING_REPEAT_GROUPS.tsv
  config/NOVELTY_FIRST_PROTEIN_FANOUT.tsv
  inputs/NOVELTY_FIRST_UNIQUE_PROTEINS.faa
  outputs/d1a_frozen/RVT_ALL_E1e4.tblout
  outputs/d1a_frozen/RVT_ALL_E1e4.domtblout
)
for f in "${req[@]}"; do test -s "$D1/$f" || { echo "[FAIL] missing $D1/$f"; exit 2; }; done
python scripts/selftest.py "$D1"
