#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
PYTHONPATH="$HERE/scripts" python - <<'PY'
from rule_core import select_families
rows=[]
for i in range(40):
    rows.append({
      'family_rep_id':f'G{i:08d}',
      'accession_n':40-i,
      'candidate_instance_n':100-i,
      'max_aa_length':100+i*10,
      'max_formal_model_n':1+(i%11),
      'max_qualifying_model_n':2+(i%13),
      'motif_c_like_family':(i%2==0),
      'motif_c_like_candidate_n':10-(i%10) if i%2==0 else 0,
      'priority_candidate_best_domain_i_evalue':10.0**(-(2+(i%8))),
      'priority_candidate_aa_length':100+i*10,
    })
sel,reason,rank,_=select_families(rows)
assert len(sel)==32 and len(set(sel))==32
assert rank[sel[0]]==1 and rank[sel[-1]]==32
assert any('REPLICATION_TOP8' in reason[x] for x in sel)
assert any('LENGTH_TOP8' in reason[x] for x in sel)
assert any('MODEL_COVERAGE_TOP8' in reason[x] for x in sel)
assert any('MOTIF_TOPOLOGY_TOP8' in reason[x] for x in sel)
print('PASS_STAGE04B2_RULE_CORE_SELFTEST')
PY
python -m py_compile "$HERE/scripts/"*.py "$HERE/frozen/historical_stage05f3b3a_exact.py"
echo PASS_STAGE04B2_PYTHON_COMPILE_SELFTEST
