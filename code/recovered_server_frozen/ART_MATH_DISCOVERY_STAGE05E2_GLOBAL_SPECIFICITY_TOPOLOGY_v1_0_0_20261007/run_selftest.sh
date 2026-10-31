#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
python -m py_compile scripts/build_stage05e2.py
python - <<'PY'
import json
p=json.load(open('config/POLICY.json'))
assert p['compact_orf_radius_bp']==250
assert p['exact_primary_assay']=='EXACT_CORE_PRIMARY_DESCRIPTION'
assert p['candidate_filtering_performed'] is False
assert p['uses_structure_results'] is False
print('[PASS] Stage05E2 policy/script selftest')
PY
