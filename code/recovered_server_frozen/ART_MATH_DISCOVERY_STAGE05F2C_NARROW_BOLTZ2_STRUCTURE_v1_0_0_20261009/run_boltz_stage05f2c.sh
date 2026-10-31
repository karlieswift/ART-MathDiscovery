#!/usr/bin/env bash
set -euo pipefail
YAML="inputs/boltz2/SF_00172800.yaml"
OUT="boltz_out"
[ -s "$YAML" ] || { echo "missing $YAML; run prepare first" >&2; exit 2; }
command -v boltz >/dev/null 2>&1 || { echo "boltz not found in current environment" >&2; exit 2; }
python - <<'PY'
import torch
print('torch=', torch.__version__)
print('cuda_available=', torch.cuda.is_available())
print('gpu_n=', torch.cuda.device_count())
if not torch.cuda.is_available(): raise SystemExit(2)
PY
mkdir -p "$OUT" logs
set -o pipefail
boltz predict "$YAML" \
  --out_dir "$OUT" \
  --model boltz2 \
  --accelerator gpu \
  --devices 1 \
  --recycling_steps 3 \
  --sampling_steps 200 \
  --diffusion_samples 1 \
  --seed 20261007 \
  --use_msa_server \
  --no_kernels \
  --output_format mmcif \
  2>&1 | tee logs/boltz.log
python scripts/postprocess_boltz.py "$OUT" | tee logs/postprocess.log
