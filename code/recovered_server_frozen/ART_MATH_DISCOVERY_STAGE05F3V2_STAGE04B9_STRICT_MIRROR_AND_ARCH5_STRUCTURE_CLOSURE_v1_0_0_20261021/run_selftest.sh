#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
python -m py_compile "$HERE"/scripts/*.py "$HERE"/frozen/historical_stage05f3b3a_exact.py
python - "$HERE" <<'PY'
import sys
from pathlib import Path
sys.path.insert(0,str(Path(sys.argv[1])/'scripts'))
from rule_core import select_families
from classify_foldseek import tier
rows=[]
for i in range(40):
 rows.append({'family_rep_id':f'F{i:03d}','accession_n':40-i,'candidate_instance_n':80-i,'max_aa_length':100+i,'max_formal_model_n':1+(i%7),'max_qualifying_model_n':2+(i%9),'motif_c_like_candidate_n':1 if i<12 else 0,'motif_c_like_family':i<12,'priority_candidate_best_domain_i_evalue':1e-8*(i+1),'priority_candidate_aa_length':100+i})
s,reason,rank,_=select_families(rows)
assert len(s)==32 and sorted(rank.values())==list(range(1,33))
assert tier(.51,.51,.51)=='STRONG' and tier(.41,.31,.31)=='SUGGESTIVE' and tier(.7,.8,.1)=='WEAK'
print('PASS_STAGE04B9_SELFTEST')
PY
sha=$(sha256sum "$HERE/frozen/historical_stage05f3b3a_exact.py"|awk '{print $1}')
[[ "$sha" == d222f638c627566909765012874c7aaa531a1a79b5e67a0d73f726f184a30ed8 ]] || { echo "HOLD bundled historical selector SHA=$sha"; exit 2; }
echo PASS_STAGE04B9_HISTORICAL_SELECTOR_BYTE_EXACT
