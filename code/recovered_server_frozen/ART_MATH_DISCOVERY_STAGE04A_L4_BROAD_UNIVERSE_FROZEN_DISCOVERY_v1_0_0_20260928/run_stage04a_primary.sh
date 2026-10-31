#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
COHORT="${1:?usage: run_stage04a_primary.sh COHORT_TSV [GENOME_WORKERS=2] [INNER_THREADS=16] [KEEP_RAW=0]}"
GW="${2:-2}"
IT="${3:-16}"
KR="${4:-0}"
cd "$ROOT"
mkdir -p outputs logs inputs work
printf '%s\n' '[0/5] byte-exact Stage03A frozen-engine preflight'
python scripts/00_preflight.py --root "$ROOT" | tee logs/00_preflight.log
printf '%s\n' '[1/5] freeze Stage04 universe; no ART truth/annotation'
python scripts/02_freeze_universe.py --root "$ROOT" --cohort-tsv "$COHORT" | tee logs/01_universe.log
printf '%s\n' '[2/5] run frozen primary engine across all sequences'
python scripts/03_run_primary.py --root "$ROOT" --genome-workers "$GW" --inner-threads "$IT" --keep-raw "$KR" | tee logs/02_primary.log
printf '%s\n' '[3/5] summarize candidate burden and retain all union regions'
python scripts/04_summarize_primary.py --root "$ROOT" | tee logs/03_summary.log
printf '%s\n' '[4/5] manifest'
python scripts/05_manifest.py --root "$ROOT" | tee logs/04_manifest.log
printf '%s\n' '[PASS] Stage04A primary scoring complete. Do not annotate candidates until outputs/STAGE04A_FINAL_GATE.json is archived and hashed.'
