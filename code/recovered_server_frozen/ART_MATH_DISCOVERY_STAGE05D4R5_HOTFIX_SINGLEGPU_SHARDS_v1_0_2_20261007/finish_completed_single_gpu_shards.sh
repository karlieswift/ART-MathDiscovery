#!/usr/bin/env bash
set -euo pipefail

# Collector-only hotfix for Boltz 2.2.x output nesting.
# Reuses already completed gpu0..gpu3 workers; does NOT rerun prediction.

for g in 0 1 2 3; do
  root="boltz_out_shards/gpu${g}"
  [ -d "$root" ] || { echo "[FAIL] missing shard root $root"; exit 2; }

  mapfile -t pdirs < <(find "$root" -type d -name predictions | LC_ALL=C sort)
  if [ "${#pdirs[@]}" -ne 1 ]; then
    echo "[FAIL] expected exactly one predictions directory under $root; found ${#pdirs[@]}" >&2
    printf '  %s\n' "${pdirs[@]:-<none>}" >&2
    echo "[INFO] directory tree (depth<=4):" >&2
    find "$root" -maxdepth 4 -type d | LC_ALL=C sort >&2
    exit 3
  fi
  echo "[FOUND] gpu${g} predictions=${pdirs[0]}"
done

# Preserve any existing canonical assembly; prediction shard outputs remain untouched.
if [ -e boltz_out ]; then
  ts=$(date +%Y%m%d_%H%M%S)
  mv boltz_out "boltz_out_precollector_${ts}"
  echo "[PRESERVE] prior boltz_out -> boltz_out_precollector_${ts}"
fi
mkdir -p boltz_out/predictions outputs

for g in 0 1 2 3; do
  root="boltz_out_shards/gpu${g}"
  pdir=$(find "$root" -type d -name predictions | LC_ALL=C sort | head -n 1)
  while IFS= read -r d; do
    [ -n "$d" ] || continue
    name=$(basename "$d")
    [ ! -e "boltz_out/predictions/$name" ] || { echo "[FAIL] duplicate prediction dir $name"; exit 4; }
    ln -s "$(readlink -f "$d")" "boltz_out/predictions/$name"
  done < <(find "$pdir" -mindepth 1 -maxdepth 1 -type d | LC_ALL=C sort)
done

touch boltz_out/.single_gpu_shards_complete
n_cif=$(find -L boltz_out/predictions -type f -name '*_model_0.cif' | wc -l)
n_conf=$(find -L boltz_out/predictions -type f -name 'confidence_*_model_0.json' | wc -l)
n_pred_dirs=$(find -L boltz_out/predictions -mindepth 1 -maxdepth 1 -type d | wc -l)

echo "[COUNT] prediction_dirs=${n_pred_dirs} CIF=${n_cif} confidence=${n_conf}"
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
  'n_prediction_dirs':int('$n_pred_dirs'),
  'implementation_hotfix':'v1.0.2 collector-only: discover Boltz 2.2.x nested predictions directories recursively under each completed single-GPU shard; no predictions rerun',
  'prediction_workers_rerun':False,
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
Path('outputs/GATE_STAGE05D4R5_PREDICTION.json').write_text(json.dumps(obj,indent=2)+'\n')
PY
cat outputs/GATE_STAGE05D4R5_PREDICTION.json
echo '[PASS] 43/43 P1 full-3D predictions assembled from completed single-GPU shards'
