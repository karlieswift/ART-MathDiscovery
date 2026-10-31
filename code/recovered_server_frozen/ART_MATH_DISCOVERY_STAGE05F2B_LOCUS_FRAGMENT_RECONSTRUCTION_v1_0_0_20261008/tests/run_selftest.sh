#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
export STAGE05F2B_SELFTEST_HERE="$HERE"
python - <<'PY2'
import importlib.util, pathlib, os
p=pathlib.Path(os.environ['STAGE05F2B_SELFTEST_HERE'])/'scripts/stage05f2b.py'
spec=importlib.util.spec_from_file_location('m',p); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
assert m.translate('ATGAAATAA')=='MK*'
assert m.union_len([(1,10),(8,12),(20,25)])==18
seq='A'*30+'ATG'*40+'TAA'+'GCT'*30+'TGA'+'C'*30
x=m.segments('x',seq,0,len(seq),25)
assert any(r['aa_length']>=30 for r in x)
print('[PASS] Stage05F2B primitive reconstruction selftest')
PY2
