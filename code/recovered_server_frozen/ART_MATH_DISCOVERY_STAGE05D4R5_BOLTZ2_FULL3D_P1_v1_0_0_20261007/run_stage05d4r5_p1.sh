#!/usr/bin/env bash
set -euo pipefail
THREADS=${1:-8}
: "${BOLTZ_CACHE:=$PWD/boltz_cache}"
export BOLTZ_CACHE
[ -d inputs/p1_yaml ] || { echo '[FAIL] run freeze_stage05d4r5.sh first'; exit 2; }
mkdir -p boltz_out logs
CUDA_VISIBLE_DEVICES=0,1,2,3 boltz predict inputs/p1_yaml --out_dir boltz_out --cache "$BOLTZ_CACHE" --model boltz2 --no_kernels --use_msa_server --devices 4 --accelerator gpu --recycling_steps 3 --sampling_steps 200 --diffusion_samples 1 --output_format mmcif --seed 20261007 --preprocessing-threads "$THREADS" 2>&1 | tee logs/boltz_p1.log
n=$(find boltz_out/predictions -name '*_model_0.cif' | wc -l)
[ "$n" -eq 43 ] || { echo "[FAIL] expected 43 P1 CIFs, got $n"; exit 3; }
echo '[PASS] 43/43 P1 full-3D predictions present'
