#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -lt 3 ]; then
  echo "usage: $0 <E3B_RESULTS.tar.gz|dir> <E3A_RESULTS.tar.gz|dir> <PROJECT_ROOT>" >&2
  exit 2
fi
E3B="$1"; E3A="$2"; P="$3"
HERE="$(cd "$(dirname "$0")" && pwd)"
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
