#!/usr/bin/env bash
set -euo pipefail
P=${1:?project root}; THREADS=${2:-32}; HERE=$(cd "$(dirname "$0")" && pwd)
command -v foldseek >/dev/null || { echo 'HOLD foldseek not found; activate art_novelty'; exit 2; }
RTDB="$P/ART_MATH_DISCOVERY_STAGE05D4R5C_CLEAN_RT_CALIBRATION_v1_0_0_20261007/references/clean_rt_afsp"
ls "${RTDB}"* >/dev/null 2>&1 || { echo "HOLD clean RT DB missing: $RTDB"; exit 2; }
mkdir -p "$HERE/work/foldseek" "$HERE/work/foldseek_tmp" "$HERE/logs"
while IFS=$'\t' read -r q rest; do
  [[ "$q" == query_id ]] && continue
  cif="$HERE/outputs/structures/${q}.cif"
  [[ -s "$cif" ]] || { echo "HOLD missing CIF $q"; exit 2; }
  out="$HERE/work/foldseek/${q}.m8"; tmp="$HERE/work/foldseek_tmp/${q}"
  rm -rf "$tmp"
  foldseek easy-search "$cif" "$RTDB" "$out" "$tmp" -e 1.0 --max-seqs 100 --threads "$THREADS" --exact-tmscore 1 --format-output query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore > "$HERE/logs/foldseek_${q}.log" 2>&1
  echo "[PASS FOLDSEEK] $q"
done < "$HERE/outputs/DELTA13_STRUCTURAL_INPUT_MANIFEST.tsv"
python "$HERE/scripts/stage04b7.py" finalize --home "$HERE"
echo PASS_STAGE04B7_ATOMIC_ANALYSIS_CLOSEOUT
