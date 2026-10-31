#!/usr/bin/env bash
set -euo pipefail
: "${BOLTZ_CACHE:=$PWD/boltz_cache}"
export BOLTZ_CACHE
mkdir -p smoke_input smoke_out logs
q=$(find inputs/p1_yaml -maxdepth 1 -name '*.yaml' | sort | head -1)
cp -f "$q" smoke_input/
CUDA_VISIBLE_DEVICES=0 boltz predict smoke_input --out_dir smoke_out --cache "$BOLTZ_CACHE" --model boltz2 --no_kernels --use_msa_server --devices 1 --accelerator gpu --recycling_steps 3 --sampling_steps 200 --diffusion_samples 1 --output_format mmcif --seed 20261007 --preprocessing-threads 1 2>&1 | tee logs/smoke_one.log
n=$(find smoke_out/predictions -name '*_model_0.cif' 2>/dev/null | wc -l)
[ "$n" -eq 1 ] || { echo "[FAIL] expected 1 CIF, got $n"; exit 2; }
echo '[PASS] Boltz-2 one-query smoke test'
