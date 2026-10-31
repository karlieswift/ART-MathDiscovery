#!/usr/bin/env bash
set -euo pipefail
B7=${1:?B7 root}
P=${2:?project root}
GPUCSV=${3:-0,1,2,3}
B7=$(cd "$B7" && pwd)
command -v boltz >/dev/null || { echo 'HOLD boltz not found; activate art_structure'; exit 2; }
MAN="$B7/outputs/DELTA13_STRUCTURAL_INPUT_MANIFEST.tsv"
[[ -s "$MAN" ]] || { echo "HOLD missing $MAN"; exit 2; }
CACHE="$P/ART_MATH_DISCOVERY_STAGE05D4R5_BOLTZ2_FULL3D_P1_v1_0_0_20261007/boltz_cache"
[[ -d "$CACHE" ]] && export BOLTZ_CACHE="$CACHE"
mkdir -p "$B7/outputs/structures" "$B7/work/boltz_retry" "$B7/logs/b7r1_retry"
# Recover any CIFs already present in raw shard outputs.
while IFS=$'\t' read -r q rest; do
  [[ "$q" == query_id ]] && continue
  [[ -s "$B7/outputs/structures/${q}.cif" ]] && continue
  cif=$(find "$B7/work/boltz_shards" -type f -name "${q}_model_0.cif" 2>/dev/null | head -1 || true)
  [[ -n "$cif" ]] && cp -f "$cif" "$B7/outputs/structures/${q}.cif"
done < "$MAN"

echo '===== ORIGINAL SHARD DIAGNOSTIC ====='
for si in 0 1 2 3; do
  n=$(find "$B7/work/boltz_shards/shard${si}" -type f -name 'B7Q*_model_0.cif' 2>/dev/null | wc -l)
  echo "shard${si}_cif_n=$n"
  if [[ -s "$B7/logs/boltz_shard${si}.log" ]]; then
    grep -Ei 'traceback|error|exception|cuda|out of memory|oom|msa|429|503|failed|killed' "$B7/logs/boltz_shard${si}.log" | tail -n 12 || true
  fi
done
mapfile -t MISS < <(tail -n +2 "$MAN" | cut -f1 | while read -r q; do [[ -s "$B7/outputs/structures/${q}.cif" ]] || echo "$q"; done)
echo "MISSING_N=${#MISS[@]}"
printf '%s\n' "${MISS[@]:-}"
if [[ ${#MISS[@]} -eq 0 ]]; then
  echo PASS_STAGE04B7R1_NO_MISSING_BOLTZ
  exit 0
fi
IFS=',' read -r -a GPUS <<< "$GPUCSV"
[[ ${#GPUS[@]} -ge 1 ]] || { echo HOLD no GPUs; exit 2; }
# Retry each missing query independently. This is a technical retry only; protocol is byte-equivalent scientifically.
i=0
for q in "${MISS[@]}"; do
  y="$B7/inputs/boltz_yaml/${q}.yaml"
  [[ -s "$y" ]] || { echo "HOLD missing YAML $y"; exit 2; }
  gpu=${GPUS[$((i % ${#GPUS[@]}))]}; i=$((i+1))
  idir="$B7/work/boltz_retry/${q}_input"; odir="$B7/work/boltz_retry/${q}_out"
  rm -rf "$idir" "$odir"; mkdir -p "$idir" "$odir"
  ln -s "$y" "$idir/${q}.yaml"
  echo "[RETRY] q=$q gpu=$gpu"
  rc=0
  (
    export CUDA_VISIBLE_DEVICES="$gpu"
    boltz predict "$idir" --out_dir "$odir" --model boltz2 --no_kernels --use_msa_server --seed 20261007 --recycling_steps 3 --sampling_steps 200 --diffusion_samples 1 --devices 1 --accelerator gpu --output_format mmcif --preprocessing-threads 8
  ) > "$B7/logs/b7r1_retry/${q}.log" 2>&1 || rc=$?
  cif=$(find "$odir" -type f -name "${q}_model_0.cif" | head -1 || true)
  if [[ -n "$cif" ]]; then
    cp -f "$cif" "$B7/outputs/structures/${q}.cif"
    echo "[PASS RETRY] $q gpu=$gpu rc=$rc"
  else
    echo "[HOLD RETRY] $q gpu=$gpu rc=$rc"
    tail -n 40 "$B7/logs/b7r1_retry/${q}.log" || true
  fi
done
n=$(find "$B7/outputs/structures" -maxdepth 1 -name 'B7Q*.cif' | wc -l)
echo "BOLTZ_PREDICTED_AFTER_RETRY=$n/13"
if [[ "$n" -ne 13 ]]; then
  echo HOLD_STAGE04B7R1_BOLTZ_STILL_INCOMPLETE
  exit 2
fi
echo PASS_STAGE04B7R1_BOLTZ13_TECHNICAL_RETRY_COMPLETE
