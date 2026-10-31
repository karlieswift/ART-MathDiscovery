#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
python -m py_compile scripts/stage05f2d.py
python - <<'PY'
import importlib.util
from pathlib import Path
p=Path('scripts/stage05f2d.py')
spec=importlib.util.spec_from_file_location('m',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
r={'evalue':'1e-8','qcov':'0.80','fident':'0.31'}
assert m.hit_class(r,False)=='CORE_SEQUENCE_HOMOLOG'
r={'evalue':'1e-4','qcov':'0.60','fident':'0.12'}
assert m.hit_class(r,False)=='BROAD_SEQUENCE_HOMOLOG_SENSITIVITY'
assert m.hit_class(r,True)=='EXACT_AA'
assert m.longer(100,130)
assert m.major_ext(100,180)
assert not m.major_ext(205,300)
print('[PASS] Stage05F2D classification/protocol selftest')
PY
python - <<'PY'
import json
x=json.load(open('config/FROZEN_POLICY.json'))
assert x['expected_query_n']==30
assert x['expected_target_unique_aa_n']==380969
assert x['expected_target_instance_n']==498879
assert x['candidate_filtering_performed'] is False
print('[PASS] Stage05F2D frozen policy selftest')
PY
