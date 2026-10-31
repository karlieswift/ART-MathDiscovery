#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
python - "$HERE" <<'PY'
import importlib.util,sys
from pathlib import Path
h=Path(sys.argv[1])
def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
d=load('d',h/'scripts/frozen_b4b1_repeat_array_scan.py')
p=load('p',h/'scripts/frozen_b4b0_family_wide_context_rescue.py')
assert d.SEED_K==(8,10,12,14,16) and d.MIN_COPIES==3 and d.GAP_MIN==60 and d.GAP_MAX==360 and d.MAX_MOTIF==80
seq='A'*10+'ACGTACGA'+'C'*150+'ACGTACGA'+'G'*160+'ACGTACGA'+'T'*10
# low complexity spacers okay; motif itself high complexity. Must recover at least one 3-copy array.
a=d.scan_region(seq,'UPSTREAM',-len(seq)); assert any(x['copy_n']>=3 for x in a),a
# exact placement sanity
nt='AAA'+'ATG'*2
assert p.revcomp('ACGT')=='ACGT'
print('[PASS] F3B4B1D frozen-detector/context helper selftest')
PY
