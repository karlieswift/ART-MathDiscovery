#!/usr/bin/env bash
set -euo pipefail

# Stage05D4-R5 v1.0.3 collector/layout hotfix.
# Converts the v1.0.2 canonical assembly from symlinked prediction directories
# into real canonical directories containing file-level symlinks.
# No Boltz prediction is rerun and no scientific setting is changed.

[ -d boltz_out/predictions ] || { echo '[FAIL] missing boltz_out/predictions'; exit 2; }

n_cif_before=$(find -L boltz_out/predictions -type f -name '*_model_0.cif' | wc -l)
n_conf_before=$(find -L boltz_out/predictions -type f -name 'confidence_*_model_0.json' | wc -l)
[ "$n_cif_before" -eq 43 ] || { echo "[FAIL] expected 43 existing CIFs before layout fix, got $n_cif_before"; exit 3; }
[ "$n_conf_before" -eq 43 ] || { echo "[FAIL] expected 43 existing confidence JSONs before layout fix, got $n_conf_before"; exit 3; }

tmp="boltz_out_layoutfix_tmp_$$"
rm -rf "$tmp"
mkdir -p "$tmp/predictions"

n_dirs=0
for d in boltz_out/predictions/*; do
  [ -d "$d" ] || continue
  q=$(basename "$d")
  mkdir -p "$tmp/predictions/$q"
  while IFS= read -r f; do
    [ -n "$f" ] || continue
    b=$(basename "$f")
    ln -sfn "$(readlink -f "$f")" "$tmp/predictions/$q/$b"
  done < <(find -L "$d" -mindepth 1 -maxdepth 1 -type f | LC_ALL=C sort)
  n_dirs=$((n_dirs+1))
done

n_cif=$(find -L "$tmp/predictions" -type f -name '*_model_0.cif' | wc -l)
n_conf=$(find -L "$tmp/predictions" -type f -name 'confidence_*_model_0.json' | wc -l)
[ "$n_dirs" -eq 43 ] || { echo "[FAIL] expected 43 prediction dirs, got $n_dirs"; rm -rf "$tmp"; exit 4; }
[ "$n_cif" -eq 43 ] || { echo "[FAIL] expected 43 CIFs after layout fix, got $n_cif"; rm -rf "$tmp"; exit 4; }
[ "$n_conf" -eq 43 ] || { echo "[FAIL] expected 43 confidence JSONs after layout fix, got $n_conf"; rm -rf "$tmp"; exit 4; }

# Preserve the v1.0.2 symlink-directory assembly.
ts=$(date +%Y%m%d_%H%M%S)
mv boltz_out "boltz_out_pre_v103_layout_${ts}"
mv "$tmp" boltz_out
# preserve the completion marker and prediction gate from the prior assembly if present
cp -f "boltz_out_pre_v103_layout_${ts}/.single_gpu_shards_complete" boltz_out/.single_gpu_shards_complete 2>/dev/null || true

mkdir -p outputs
python - <<PY
import json
from pathlib import Path
obj={
  'stage':'Stage05D4-R5-prediction-layout',
  'gate':'PASS_BOLTZ2_P1_43_CANONICAL_REALDIR_LAYOUT',
  'gate_pass':True,
  'n_P1':43,
  'n_prediction_dirs':43,
  'n_cif':43,
  'n_confidence_json':43,
  'implementation_hotfix':'v1.0.3 canonical-layout-only: replace symlinked prediction directories with real canonical directories containing file-level symlinks so downstream find/pathlib consumers can traverse them without -L',
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
  'no_kernels':True
}
Path('outputs/GATE_STAGE05D4R5_LAYOUT_v1_0_3.json').write_text(json.dumps(obj,indent=2)+'\n')
PY

cat outputs/GATE_STAGE05D4R5_LAYOUT_v1_0_3.json
echo '[PASS] canonical Boltz prediction layout repaired without rerunning predictions'
