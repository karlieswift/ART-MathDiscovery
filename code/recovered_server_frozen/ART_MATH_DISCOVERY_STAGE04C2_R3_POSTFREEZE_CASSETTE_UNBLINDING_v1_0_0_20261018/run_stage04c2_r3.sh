#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <Stage04C2_R2D_root> <Stage04C2_R2C_root> <Stage04C2_UG27_source_root> <R2D_frozen_result_archive.tar.gz>" >&2
  exit 2
fi
ROOT="$(cd "$(dirname "$0")" && pwd)"
R2D="$(readlink -f "$1")"
R2C="$(readlink -f "$2")"
UG="$(readlink -f "$3")"
ARC="$(readlink -f "$4")"
mkdir -p "$ROOT/outputs" "$ROOT/logs"
rm -f "$ROOT/outputs"/*.tsv "$ROOT/outputs"/*.json

echo '[0/4] verify archived + live R2D clean board; truth not read'
python "$ROOT/scripts/00_score_archive_preflight.py" --root "$ROOT" --r2d "$R2D" --r2c "$R2C" --r2d-archive "$ARC" 2>&1 | tee "$ROOT/logs/00_score_archive_preflight.log"

echo '[1/4] open previously frozen Stage04C2B-R3 bridge and lock five cassette coordinates'
python "$ROOT/scripts/01_truth_bridge_lock.py" --root "$ROOT" --r2c "$R2C" --ug "$UG" 2>&1 | tee "$ROOT/logs/01_truth_bridge_lock.log"

echo '[2/4] post-freeze cassette-level unblinding against frozen R2D boards'
python "$ROOT/scripts/02_cassette_unblind.py" --root "$ROOT" --r2d "$R2D" 2>&1 | tee "$ROOT/logs/02_cassette_unblind.log"

echo '[3/4] finalize no-retune cassette-level gate + hashes'
python "$ROOT/scripts/03_finalize.py" --root "$ROOT" 2>&1 | tee "$ROOT/logs/03_finalize.log"

echo '[PASS] Stage04C2-R3 post-freeze cassette-level unblinding complete'
