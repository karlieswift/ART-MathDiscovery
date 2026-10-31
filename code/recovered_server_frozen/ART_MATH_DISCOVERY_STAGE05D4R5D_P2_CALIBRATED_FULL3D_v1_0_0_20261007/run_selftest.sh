#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd); cd "$ROOT"
python - <<'PY'
import json,py_compile
p=json.load(open('config/POLICY.json'))
assert p['source_r5_p2_expected']==40
assert p['boltz_model']=='boltz2' and p['no_kernels'] is True and p['msa_server'] is True
assert p['strong_rt_fold_rule'].startswith('alntmscore >= 0.50')
assert p['suggestive_rt_fold_rule'].startswith('alntmscore >= 0.40')
for x in ['scripts/prepare_p2.py','scripts/collect_boltz_p2.py','scripts/summarize_p2.py']:py_compile.compile(x,doraise=True)
print('[PASS] Stage05D4-R5D policy/script selftest')
PY
