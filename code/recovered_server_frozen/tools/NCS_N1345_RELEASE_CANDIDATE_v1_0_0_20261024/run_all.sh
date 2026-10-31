#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
OUT="${1:-$ROOT/results}"
mkdir -p "$OUT/N1" "$OUT/N3" "$OUT/N4"
python "$ROOT/src/n1_benchmark.py" admit --intake "$ROOT/data/N1_ART_NEW_INTAKE.tsv" --exposure "$ROOT/data/EXPOSURE_919_AND_HISTORICAL_ACCESSIONS.tsv" --out "$OUT/N1/N1_ADMISSION.tsv" | tee "$OUT/N1/N1_STATUS.json"
python "$ROOT/src/n3_analysis.py" --region-ranks "$ROOT/data/N3_REGION_RANKS_6_NEW_MATRICES.tsv.gz" --out "$OUT/N3" | tee "$OUT/N3/N3_EXEC_LOG.json"
python "$ROOT/src/n4_analysis.py" --review "$ROOT/data/N4_Review69_SANITIZED.tsv" --pairs "$ROOT/data/N4_Pairs2144_SANITIZED.tsv" --out "$OUT/N4" | tee "$OUT/N4/N4_EXEC_LOG.json"
python "$ROOT/src/release_check.py" --root "$ROOT" --out "$OUT/N5_VERIFY.json"
python "$ROOT/src/verify_expected.py" "$ROOT" "$OUT"
printf 'RESULTS_DIRECTORY=%s\n' "$OUT"
