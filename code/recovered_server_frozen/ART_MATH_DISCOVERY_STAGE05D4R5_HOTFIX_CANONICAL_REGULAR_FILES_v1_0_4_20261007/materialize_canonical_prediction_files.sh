#!/usr/bin/env bash
set -euo pipefail

R5_DIR="${1:-$PWD}"
cd "$R5_DIR"

SRC_ROOT="boltz_out/predictions"
if [[ ! -d "$SRC_ROOT" ]]; then
  echo "[FAIL] missing $SRC_ROOT" >&2
  exit 2
fi

# Require the already-repaired 43 canonical query directories.
mapfile -t QDIRS < <(find "$SRC_ROOT" -mindepth 1 -maxdepth 1 -type d | sort)
if [[ ${#QDIRS[@]} -ne 43 ]]; then
  echo "[FAIL] expected 43 canonical query directories, got ${#QDIRS[@]}" >&2
  exit 3
fi

STAMP=$(date +%Y%m%d_%H%M%S)
NEW_ROOT="boltz_out_v104_real/predictions"
rm -rf boltz_out_v104_real
mkdir -p "$NEW_ROOT"

n_cif=0
n_conf=0
n_dirs=0
for qdir in "${QDIRS[@]}"; do
  q=$(basename "$qdir")
  out="$NEW_ROOT/$q"
  mkdir -p "$out"
  n_dirs=$((n_dirs+1))

  shopt -s nullglob
  files=("$qdir"/*)
  shopt -u nullglob
  if [[ ${#files[@]} -eq 0 ]]; then
    echo "[FAIL] empty canonical query directory: $qdir" >&2
    exit 4
  fi

  for f in "${files[@]}"; do
    # Resolve file-level symlinks produced by v1.0.3. Also accepts regular files.
    if [[ -L "$f" ]]; then
      src=$(readlink -f "$f")
    elif [[ -f "$f" ]]; then
      src="$f"
    else
      continue
    fi
    if [[ ! -f "$src" ]]; then
      echo "[FAIL] unresolved/non-regular source for $f -> $src" >&2
      exit 5
    fi
    dest="$out/$(basename "$f")"
    # Prefer a hard link so the canonical path is a true regular file without duplicating bytes.
    # Fall back to a real copy if hard-linking is unavailable.
    if ! ln "$src" "$dest" 2>/dev/null; then
      cp -a --reflink=auto "$src" "$dest"
    fi
  done

  cifs=("$out"/*_model_0.cif)
  confs=("$out"/confidence_*_model_0.json)
  [[ ${#cifs[@]} -eq 1 && -f "${cifs[0]}" && ! -L "${cifs[0]}" ]] || {
    echo "[FAIL] expected one regular CIF in $out" >&2; exit 6; }
  [[ ${#confs[@]} -eq 1 && -f "${confs[0]}" && ! -L "${confs[0]}" ]] || {
    echo "[FAIL] expected one regular confidence JSON in $out" >&2; exit 7; }
  n_cif=$((n_cif+1))
  n_conf=$((n_conf+1))
done

# Verify using the exact semantics used by downstream consumers: no -L, regular files only.
plain_cif=$(find boltz_out_v104_real/predictions -type f -name '*_model_0.cif' | wc -l)
plain_conf=$(find boltz_out_v104_real/predictions -type f -name 'confidence_*_model_0.json' | wc -l)
plain_dirs=$(find boltz_out_v104_real/predictions -mindepth 1 -maxdepth 1 -type d | wc -l)

if [[ "$plain_dirs" -ne 43 || "$plain_cif" -ne 43 || "$plain_conf" -ne 43 ]]; then
  echo "[FAIL] regular-file verification dirs=$plain_dirs CIF=$plain_cif confidence=$plain_conf" >&2
  exit 8
fi

if [[ -d boltz_out ]]; then
  mv boltz_out "boltz_out_pre_v104_regular_${STAMP}"
fi
mv boltz_out_v104_real boltz_out

mkdir -p outputs
cat > outputs/GATE_STAGE05D4R5_PREDICTION_LAYOUT_V104.json <<JSON
{
  "stage": "Stage05D4-R5-prediction-layout",
  "gate": "PASS_BOLTZ2_P1_43_CANONICAL_REGULAR_FILE_LAYOUT",
  "gate_pass": true,
  "n_P1": 43,
  "n_prediction_dirs": 43,
  "n_cif_regular_files": 43,
  "n_confidence_json_regular_files": 43,
  "implementation_hotfix": "v1.0.4 layout-only: materialize v1.0.3 file-level symlinks as true regular files (hard links when possible; copy fallback) so downstream find -type f/pathlib consumers see all 43 outputs without -L",
  "prediction_workers_rerun": false,
  "scientific_settings_changed": false,
  "candidate_filtering_performed": false,
  "dna_score_recomputed": false,
  "dna_rank_recomputed": false,
  "boltz_model": "boltz2",
  "seed": 20261007,
  "recycling_steps": 3,
  "sampling_steps": 200,
  "diffusion_samples": 1,
  "msa_server": true,
  "no_kernels": true
}
JSON

cat outputs/GATE_STAGE05D4R5_PREDICTION_LAYOUT_V104.json
printf '[PASS] canonical regular-file layout: dirs=%s CIF=%s confidence=%s\n' "$plain_dirs" "$plain_cif" "$plain_conf"
