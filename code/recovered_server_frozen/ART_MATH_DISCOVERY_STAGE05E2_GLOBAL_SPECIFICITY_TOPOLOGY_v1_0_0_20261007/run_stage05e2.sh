#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <D1_ROOT> <E1_ROOT> <R1_RESULTS_TAR_OR_DIR>" >&2
  exit 2
fi
D1=$(readlink -f "$1")
E1=$(readlink -f "$2")
R1=$(readlink -f "$3")
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
rm -rf outputs work/r1
mkdir -p outputs work
if [ -d "$R1" ]; then
  FAM=$(find "$R1" -type f -name STAGE05D4R1_FAM50_CLUSTERS.tsv -print -quit)
else
  mkdir -p work/r1
  tar -xzf "$R1" -C work/r1
  FAM=$(find work/r1 -type f -name STAGE05D4R1_FAM50_CLUSTERS.tsv -print -quit)
fi
[ -n "${FAM:-}" ] && [ -f "$FAM" ] || { echo '[FAIL] FAM50 table not found' >&2; exit 3; }
python scripts/build_stage05e2.py --d1-root "$D1" --e1-root "$E1" --fam50-tsv "$FAM" --outdir outputs
