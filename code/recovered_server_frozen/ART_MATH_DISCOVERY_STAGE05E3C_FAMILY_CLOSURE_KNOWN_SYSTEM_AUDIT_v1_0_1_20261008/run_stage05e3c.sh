#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -lt 3 ]; then
  echo "usage: $0 <E3B_RESULTS.tar.gz|dir> <E3A_RESULTS.tar.gz|dir> <PROJECT_ROOT>" >&2
  exit 2
fi
E3B="$1"; E3A="$2"; P="$3"
HERE="$(cd "$(dirname "$0")" && pwd)"

resolve_input() {
  local requested="$1"
  local project_root="$2"
  local label="$3"
  if [ -e "$requested" ]; then
    printf '%s\n' "$requested"
    return 0
  fi
  local base
  base="$(basename "$requested")"
  mapfile -d '' hits < <(find "$project_root" -type f -name "$base" -print0 2>/dev/null || true)
  if [ "${#hits[@]}" -eq 1 ]; then
    echo "[INFO] $label input not at requested path; uniquely resolved by basename:" >&2
    echo "       ${hits[0]}" >&2
    printf '%s\n' "${hits[0]}"
    return 0
  fi
  echo "[FAIL] cannot resolve $label input: $requested" >&2
  if [ "${#hits[@]}" -eq 0 ]; then
    echo "       no file named $base found below $project_root" >&2
  else
    echo "       multiple files named $base found; pass the exact path:" >&2
    printf '       %s\n' "${hits[@]}" >&2
  fi
  return 3
}

E3B="$(resolve_input "$E3B" "$P" E3B)"
E3A="$(resolve_input "$E3A" "$P" E3A)"

cd "$HERE"
python scripts/stage05e3c.py \
  --e3b "$E3B" \
  --e3a "$E3A" \
  --project-root "$P" \
  --outdir outputs \
  --workdir work \
  --expected-families 12 \
  --art95-manifest "${E3C_ART95_MANIFEST:-}" \
  --hmm-dir "${E3C_KNOWN_HMM_DIR:-}"
python - <<'PY'
import json
x=json.load(open('outputs/STAGE05E3C_GATE.json'))
assert x['gate_pass'] is True
assert x['observed_family_n']==12
assert x['family_membership_changed'] is False
assert x['dna_rank_recomputed'] is False
print('[PASS] Stage05E3C gate')
PY
