#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <Stage05D1_root> <Stage05D4R1_results_dir_or_tar.gz>" >&2
  exit 2
fi
D1=$(realpath "$1"); R1IN=$(realpath "$2")
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
rm -rf work outputs logs
mkdir -p work outputs logs
if [ -d "$R1IN" ]; then
  R1="$R1IN"
else
  mkdir -p work/r1
  tar -xzf "$R1IN" -C work/r1
  R1="$HERE/work/r1"
fi
FAM="$R1/outputs/STAGE05D4R1_FAM50_CLUSTERS.tsv"
[ -f "$FAM" ] || { echo "[FAIL] missing $FAM"; exit 1; }
python scripts/build_architecture_families.py \
  --d1-root "$D1" \
  --fam50-tsv "$FAM" \
  --outdir outputs \
  2>&1 | tee logs/run_stage05e1.log
cat outputs/GATE_STAGE05E1.json
