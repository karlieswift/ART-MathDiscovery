#!/usr/bin/env bash
set -euo pipefail
D1="${1:-}"
python tests/selftest.py
if [[ -n "$D1" ]]; then
  for f in \
    config/NOVELTY_FIRST_999_CANDIDATES.tsv \
    config/NOVELTY_FIRST_CANDIDATE_ORF_LINKS.tsv \
    config/NOVELTY_FIRST_24020_ORF_REGISTRY.tsv \
    inputs/NOVELTY_FIRST_UNIQUE_PROTEINS.faa; do
      [[ -s "$D1/$f" ]] || { echo "[FAIL] missing $D1/$f"; exit 2; }
  done
  echo '[PASS] Stage05D1 source files present'
fi
