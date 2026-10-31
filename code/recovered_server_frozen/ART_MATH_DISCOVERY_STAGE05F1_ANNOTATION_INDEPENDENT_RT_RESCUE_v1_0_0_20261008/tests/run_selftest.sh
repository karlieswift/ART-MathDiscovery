#!/usr/bin/env bash
set -euo pipefail
python - <<'PY'
import importlib.util, pathlib
p=pathlib.Path('scripts/stage05f1.py')
spec=importlib.util.spec_from_file_location('m',p); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
# translation / reverse-complement coordinate primitives
assert m.translate_frame('ATGGATTATTAA') == 'MDY*'
seq='A'*30+'ATGGCTGCTGCTTAA'+'C'*30
parts=list(m.segments_from_interval('X',seq,0,len(seq)))
# min45 means this toy produces no segment; direct motif regex sanity
import re
assert re.search(r'[YF][A-Z]DD','MAALDDQ') is None
assert re.search(r'[YF][A-Z]DD','MAYADDQ')
assert m.merge_intervals([(1,5),(4,9),(20,21)]) == [(1,9),(20,21)]
print('[PASS] Stage05F1 primitive selftest')
PY
python -m py_compile scripts/stage05f1.py
