#!/usr/bin/env bash
set -euo pipefail
P=${1:?project root}; GPUCSV=${2:-0,1,2,3}; HERE=$(cd "$(dirname "$0")" && pwd)
command -v boltz >/dev/null || { echo 'HOLD boltz not found; activate art_structure'; exit 2; }
CACHE="$P/ART_MATH_DISCOVERY_STAGE05D4R5_BOLTZ2_FULL3D_P1_v1_0_0_20261007/boltz_cache"
[[ -d "$CACHE" ]] && export BOLTZ_CACHE="$CACHE"
IFS=',' read -r -a GPUS <<< "$GPUCSV"
if [[ ${#GPUS[@]} -lt 1 ]]; then echo HOLD no GPUs; exit 2; fi
mkdir -p "$HERE/work/boltz_shards" "$HERE/outputs/structures" "$HERE/logs"
pids=()
for si in 0 1 2 3; do
  shdir="$HERE/inputs/shard$si"
  [[ -d "$shdir" ]] || continue
  n=$(find "$shdir" -maxdepth 1 -name '*.yaml' | wc -l)
  [[ "$n" -gt 0 ]] || continue
  gpu=${GPUS[$((si % ${#GPUS[@]}))]}
  out="$HERE/work/boltz_shards/shard$si"
  mkdir -p "$out"
  echo "[BOLTZ SHARD $si] gpu=$gpu n=$n"
  (
    export CUDA_VISIBLE_DEVICES="$gpu"
    boltz predict "$shdir" --out_dir "$out" --model boltz2 --no_kernels --use_msa_server --seed 20261007 --recycling_steps 3 --sampling_steps 200 --diffusion_samples 1 --devices 1 --accelerator gpu --output_format mmcif --preprocessing-threads 8
  ) > "$HERE/logs/boltz_shard${si}.log" 2>&1 &
  pids+=("$!")
done
rc=0
for p in "${pids[@]}"; do wait "$p" || rc=1; done
# gather CIFs independent of minor boltz output-dir version nesting
for q in $(tail -n +2 "$HERE/outputs/DELTA13_STRUCTURAL_INPUT_MANIFEST.tsv" | cut -f1); do
  cif=$(find "$HERE/work/boltz_shards" -type f -name "${q}_model_0.cif" | head -1 || true)
  if [[ -n "$cif" ]]; then cp -f "$cif" "$HERE/outputs/structures/${q}.cif"; fi
done
n=$(find "$HERE/outputs/structures" -maxdepth 1 -name 'B7Q*.cif' | wc -l)
echo "BOLTZ_PREDICTED=$n/13"
if [[ "$n" -ne 13 ]]; then echo HOLD_STAGE04B7_BOLTZ_INCOMPLETE; exit 2; fi
[[ "$rc" -eq 0 ]] || echo 'NOTE one shard command returned nonzero but all 13 CIFs were recovered'
echo PASS_STAGE04B7_BOLTZ13_COMPLETE
