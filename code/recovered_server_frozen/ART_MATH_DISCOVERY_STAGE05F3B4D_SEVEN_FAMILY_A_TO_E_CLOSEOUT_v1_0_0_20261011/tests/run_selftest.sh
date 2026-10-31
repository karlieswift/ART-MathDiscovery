#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")/.." && pwd)
python -m py_compile "$S/scripts"/*.py
python - "$S/scripts/analyze_abc.py" <<'PY'
import importlib.util,sys
p=sys.argv[1]
s=importlib.util.spec_from_file_location('m',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
assert m.lev('TTGAAGAA','ATTGAAGAA')==1
assert m.lev('ATGGAGTA','ATGGAGTA')==0
assert m.jacc(m.kmers('A'*40+'C'*40,31),m.kmers('A'*40+'C'*40,31))==1.0
cs=m.components(['a','b','c'],[('a','b')]);assert sorted(map(tuple,cs))==[('a','b'),('c',)]
print('[PASS] F3B4D A/B/C helper selftest')
PY
bash -n "$S/run_stage05f3b4d_prepare.sh"
bash -n "$S/run_stage05f3b4d_ace.sh"
bash -n "$S/run_stage05f3b4d_d_boltz.sh"
bash -n "$S/run_stage05f3b4d_d_foldseek.sh"
bash -n "$S/finalize_stage05f3b4d.sh"
bash -n "$S/collect_stage05f3b4d.sh"
echo '[PASS] F3B4D package syntax/selftest'
