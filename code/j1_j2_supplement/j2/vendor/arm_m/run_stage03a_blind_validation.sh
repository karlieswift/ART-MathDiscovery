#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
THREADS="${1:-32}"
KEEP_RAW="${2:-0}"
FASTA_SOURCE_DIR="${3:-${ART_FASTA_SOURCE_DIR:-}}"
cd "$ROOT"
mkdir -p outputs logs inputs work

echo '[0/5] freeze development engine before any holdout scoring'
python scripts/00_freeze_engine.py --root "$ROOT" | tee logs/00_engine_freeze.log

echo '[1/5] fetch/freeze raw validation FASTAs only (no annotations, no ART truth)'
FETCH_ARGS=(--root "$ROOT")
if [[ -n "$FASTA_SOURCE_DIR" ]]; then
  echo "[INFO] raw FASTA local source dir: $FASTA_SOURCE_DIR"
  FETCH_ARGS+=(--source-dir "$FASTA_SOURCE_DIR")
fi
python scripts/01_fetch_validation_fastas.py "${FETCH_ARGS[@]}" | tee logs/01_fetch_fastas.log

echo '[2/5] L3 genome-blind primary scans'
tail -n +2 config/VALIDATION_COHORT.tsv | cut -f1 | while read -r ACC; do
  echo "=== PRIMARY ${ACC} ==="
  python scripts/02_primary_scan.py --root "$ROOT" --accession "$ACC" --threads "$THREADS" --keep-raw "$KEEP_RAW" | tee "logs/${ACC}_02_primary.log"
done

echo '[3/5] all-region orthogonal confirmation on each frozen blind board'
tail -n +2 config/VALIDATION_COHORT.tsv | cut -f1 | while read -r ACC; do
  echo "=== CONFIRM ${ACC} ==="
  python scripts/03_confirmation_scan.py --root "$ROOT" --accession "$ACC" --workers "$THREADS" --keep-raw "$KEEP_RAW" | tee "logs/${ACC}_03_confirmation.log"
done

echo '[4/5] freeze cohort blind scores; still no ART truth'
python scripts/04_finalize_cohort.py --root "$ROOT" | tee logs/04_finalize.log

echo '[PASS] Stage03A blind validation complete'
echo '[SEALED] Do not inspect/import ART truth until STAGE03A_FINAL_GATE.json is PASS and returned for audit.'
