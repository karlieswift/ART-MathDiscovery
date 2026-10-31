#!/usr/bin/env bash
set -euo pipefail
THREADS_PER_WORKER=${1:-2}
: "${BOLTZ_CACHE:=$PWD/boltz_cache}"
export BOLTZ_CACHE

[ -d inputs/p1_yaml ] || { echo '[FAIL] inputs/p1_yaml missing; run freeze_stage05d4r5.sh first'; exit 2; }
command -v boltz >/dev/null || { echo '[FAIL] boltz not found'; exit 2; }

n_yaml=$(find inputs/p1_yaml -maxdepth 1 -type f -name '*.yaml' | wc -l)
[ "$n_yaml" -eq 43 ] || { echo "[FAIL] expected 43 P1 YAMLs, got $n_yaml"; exit 2; }

mkdir -p inputs/p1_shards logs boltz_out_shards outputs
rm -rf inputs/p1_shards/gpu0 inputs/p1_shards/gpu1 inputs/p1_shards/gpu2 inputs/p1_shards/gpu3
mkdir -p inputs/p1_shards/gpu0 inputs/p1_shards/gpu1 inputs/p1_shards/gpu2 inputs/p1_shards/gpu3

# Deterministic sorted round-robin partition: 11/11/11/10.
i=0
while IFS= read -r f; do
  g=$((i % 4))
  ln -sfn "$(readlink -f "$f")" "inputs/p1_shards/gpu${g}/$(basename "$f")"
  i=$((i + 1))
done < <(find inputs/p1_yaml -maxdepth 1 -type f -name '*.yaml' | LC_ALL=C sort)

for g in 0 1 2 3; do
  c=$(find "inputs/p1_shards/gpu${g}" -maxdepth 1 -type l -name '*.yaml' | wc -l)
  echo "[SHARD] gpu${g} inputs=${c}"
done

# Preserve any failed DDP attempt; never delete it.
if [ -e boltz_out ] && [ ! -e boltz_out/.single_gpu_shards_complete ]; then
  ts=$(date +%Y%m%d_%H%M%S)
  mv boltz_out "boltz_out_ddp_failed_${ts}"
  echo "[PRESERVE] prior boltz_out -> boltz_out_ddp_failed_${ts}"
fi
rm -rf boltz_out
mkdir -p boltz_out/predictions

pids=()
for g in 0 1 2 3; do
  mkdir -p "boltz_out_shards/gpu${g}"
  echo "[LAUNCH] GPU ${g}"
  (
    export CUDA_VISIBLE_DEVICES="$g"
    boltz predict "inputs/p1_shards/gpu${g}" \
      --out_dir "boltz_out_shards/gpu${g}" \
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
      --preprocessing-threads "$THREADS_PER_WORKER"
  ) > "logs/boltz_p1_gpu${g}.log" 2>&1 &
  pids[$g]=$!
done

rc=0
for g in 0 1 2 3; do
  if wait "${pids[$g]}"; then
    echo "[PASS] gpu${g} worker completed"
  else
    echo "[FAIL] gpu${g} worker failed; inspect logs/boltz_p1_gpu${g}.log" >&2
    rc=1
  fi
done
[ "$rc" -eq 0 ] || exit 3

# Assemble canonical boltz_out/predictions as symlinks to complete query folders.
for g in 0 1 2 3; do
  pdir="boltz_out_shards/gpu${g}/predictions"
  [ -d "$pdir" ] || { echo "[FAIL] missing $pdir"; exit 4; }
  while IFS= read -r d; do
    name=$(basename "$d")
    [ ! -e "boltz_out/predictions/$name" ] || { echo "[FAIL] duplicate prediction dir $name"; exit 4; }
    ln -s "$(readlink -f "$d")" "boltz_out/predictions/$name"
  done < <(find "$pdir" -mindepth 1 -maxdepth 1 -type d | LC_ALL=C sort)
done

touch boltz_out/.single_gpu_shards_complete
n_cif=$(find -L boltz_out/predictions -type f -name '*_model_0.cif' | wc -l)
n_conf=$(find -L boltz_out/predictions -type f -name 'confidence_*_model_0.json' | wc -l)

echo "[COUNT] CIF=${n_cif} confidence=${n_conf}"
[ "$n_cif" -eq 43 ] || { echo "[FAIL] expected 43 P1 CIFs, got $n_cif"; exit 5; }
[ "$n_conf" -eq 43 ] || { echo "[FAIL] expected 43 confidence JSONs, got $n_conf"; exit 5; }

python - <<PY
import json
from pathlib import Path
obj={
  'stage':'Stage05D4-R5-prediction',
  'gate':'PASS_BOLTZ2_P1_43_FULL3D_SINGLEGPU_SHARDS',
  'gate_pass':True,
  'n_P1':43,
  'n_cif':int('$n_cif'),
  'n_confidence_json':int('$n_conf'),
  'implementation_hotfix':'v1.0.1: replace Boltz/PyTorch-Lightning 4-device DDP with four deterministic independent single-GPU workers to avoid NCCL/TCPStore failure',
  'scientific_settings_changed':False,
  'candidate_filtering_performed':False,
  'dna_score_recomputed':False,
  'dna_rank_recomputed':False,
  'boltz_model':'boltz2',
  'seed':20261007,
  'recycling_steps':3,
  'sampling_steps':200,
  'diffusion_samples':1,
  'msa_server':True,
  'no_kernels':True,
}
Path('outputs').mkdir(exist_ok=True)
Path('outputs/GATE_STAGE05D4R5_PREDICTION.json').write_text(json.dumps(obj,indent=2)+'\n')
PY
cat outputs/GATE_STAGE05D4R5_PREDICTION.json
echo '[PASS] 43/43 P1 full-3D predictions present via four independent single-GPU shards'
