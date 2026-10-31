#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
python -m py_compile "$HERE/scripts/analyze_stage05e3b.py"
python - <<'PY'
import json
p=json.load(open('config/POLICY.json'))
assert p['frozen_universe']['families']==12
assert p['frozen_universe']['candidates']==119
assert p['frozen_universe']['sequence_ids']==61
assert p['frozen_universe']['exact_primary_joint_witness_candidates']==80
assert 'candidate filtering' in p['prohibited']
print('[PASS] Stage05E3B policy/script selftest')
PY
