#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"; cd "$ROOT"
TS=$(date +%Y%m%d_%H%M%S)
OUT="ART_STAGE05D1_NOVELTY_FIRST_RESULTS_${TS}.tar.gz"
FILES=(config/PREFREEZE_SUMMARY.json config/NOVELTY_FIRST_999_CANDIDATES.tsv references/broad_rt/REFERENCE_SHA256.txt references/broad_rt/REFERENCE_COUNTS.txt outputs/GATE_STAGE05D1.json outputs/SUMMARY_STAGE05D1.json outputs/NOVELTY_FIRST_BOARD.tsv outputs/MYRT_PRIMARY_DOMAINS.tsv outputs/CDD_PFAM_REFINEMENT_DOMAINS.tsv outputs/YxDD_MOTIF_ONLY.tsv outputs/PROTEIN_TO_CANDIDATE_BROAD_RT_LINKS.tsv logs/myrt.log logs/cdd.log)
tar -czf "$OUT" "${FILES[@]}"
sha256sum "$OUT" > "$OUT.sha256"
echo "$OUT"
echo "$OUT.sha256"
