#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
python -m py_compile scripts/build_calibration_panel.py scripts/summarize_r5c.py
python - <<'PY'
import json
p=json.load(open('config/POLICY.json'))
assert p['query_universe_expected']==43
assert p['negative_controls_per_rt']==3
assert p['reference_selection_uses_candidate_results'] is False
assert p['strong_rt_fold_rule'].startswith('alntmscore >= 0.50')
assert p['suggestive_rt_fold_rule'].startswith('alntmscore >= 0.40')
print('[PASS] Stage05D4-R5C policy selftest')
PY
