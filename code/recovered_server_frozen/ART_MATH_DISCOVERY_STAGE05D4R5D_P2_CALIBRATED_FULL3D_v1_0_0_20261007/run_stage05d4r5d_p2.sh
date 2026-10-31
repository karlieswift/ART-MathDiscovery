#!/usr/bin/env bash
set -euo pipefail
R5=${1:?R5_ROOT}; THREADS=${2:-2}
ROOT=$(cd "$(dirname "$0")" && pwd); cd "$ROOT"
[ -d inputs/p2_yaml ] || { echo '[FAIL] run prepare_stage05d4r5d.sh first'; exit 2; }
: "${BOLTZ_CACHE:=$R5/boltz_cache}"
export BOLTZ_CACHE
rm -rf inputs/p2_shards boltz_p2_shards boltz_p2_out
mkdir -p inputs/p2_shards/{gpu0,gpu1,gpu2,gpu3} boltz_p2_shards logs outputs
python - <<'PY'
from pathlib import Path
files=sorted(Path('inputs/p2_yaml').glob('*.yaml'))
assert len(files)==40,len(files)
for i,p in enumerate(files):
    g=i%4
    (Path('inputs/p2_shards')/f'gpu{g}'/p.name).write_bytes(p.read_bytes())
for g in range(4): print(f'[SHARD] gpu{g} inputs={len(list((Path("inputs/p2_shards")/f"gpu{g}").glob("*.yaml")))}')
PY
pids=()
for g in 0 1 2 3; do
  echo "[LAUNCH] GPU $g"
  (
    CUDA_VISIBLE_DEVICES="$g" boltz predict "inputs/p2_shards/gpu$g" \
      --out_dir "boltz_p2_shards/gpu$g" \
      --cache "$BOLTZ_CACHE" \
      --model boltz2 \
      --no_kernels \
      --use_msa_server \
      --devices 1 \
      --accelerator gpu \
      --recycling_steps 3 \
      --sampling_steps 200 \
      --diffusion_samples 1 \
      --output_format mmcif \
      --seed 20261007 \
      --preprocessing-threads "$THREADS" \
      > "logs/boltz_p2_gpu${g}.log" 2>&1
  ) & pids+=("$!")
done
fail=0
for i in 0 1 2 3; do
  if wait "${pids[$i]}"; then echo "[PASS] gpu$i worker completed"; else echo "[FAIL] gpu$i worker failed"; fail=1; fi
done
[ "$fail" -eq 0 ] || exit 3
python scripts/collect_boltz_p2.py --shards boltz_p2_shards --out boltz_p2_out
cat outputs/GATE_STAGE05D4R5D_PREDICTION.json
