#!/usr/bin/env python3
import sys,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import pandas as pd
from common import build_two_arm_board
# synthetic: two overlapping selected families merge; isolated weak ACNV creates second region
D=pd.DataFrame([
 dict(sequence_id='x',window_len=128,start0=0,end0=128,q_acnv_family=.001,q_recurrence=.5,acnv_family_best_member='A',recurrence_best_k=8,recurrence_best_lag_bp=50),
 dict(sequence_id='x',window_len=128,start0=64,end0=192,q_acnv_family=.5,q_recurrence=.002,acnv_family_best_member='A',recurrence_best_k=8,recurrence_best_lag_bp=60),
 dict(sequence_id='x',window_len=128,start0=500,end0=628,q_acnv_family=.004,q_recurrence=.8,acnv_family_best_member='B',recurrence_best_k=6,recurrence_best_lag_bp=70),
])
sel,r=build_two_arm_board(D,.005,'T')
assert len(sel)==3 and len(r)==2
assert int(r.iloc[0].support_family_n)==2 and int(r.iloc[0].start0)==0 and int(r.iloc[0].end0)==192
assert list(r.region_rank)==[1,2]
print('[PASS] clean two-arm board selftest')
