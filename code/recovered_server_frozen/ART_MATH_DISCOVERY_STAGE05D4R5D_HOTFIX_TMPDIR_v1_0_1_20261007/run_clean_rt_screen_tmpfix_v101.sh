#!/usr/bin/env bash
set -euo pipefail
R5D=${1:?R5D_ROOT}
R5C=${2:?R5C_ROOT}
THREADS=${3:-32}
cd "$R5D"
# Implementation-only hotfix: Foldseek easy-search needs the parent tmp/ directory.
# The frozen R5D scientific settings, P2 universe, Boltz structures, clean-RT panel,
# Foldseek parameters and thresholds are unchanged.
if [ -e tmp ] && [ ! -d tmp ]; then
  stamp=$(date +%Y%m%d_%H%M%S)
  mv tmp "tmp_non_dir_${stamp}"
  echo "[PRESERVE] non-directory tmp -> tmp_non_dir_${stamp}"
fi
mkdir -p tmp
rm -rf tmp/p2_clean_rt
exec ./run_clean_rt_screen.sh "$R5C" "$THREADS"
