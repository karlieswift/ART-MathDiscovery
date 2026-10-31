#!/usr/bin/env python3
import csv, tempfile, sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent))
import stage04b6 as s
H=Path(__file__).resolve().parents[1]
sel=s.tsv_read(H/'frozen/ARCHITECTURE_SELECTED_FAMILIES.tsv')
assert len(sel)==41
assert sum(r['arm']=='MIXED_REPLICATION_ARM' for r in sel)==23
assert sum(r['arm']=='DELTA_ONLY_NOVELTY_ARM' for r in sel)==18
m=s.motif_hits('AAAYMDDTTTFVDD')
assert [x[1] for x in m]==['YMDD','FVDD']
# synthetic exact unitig rescue: candidate sits inside stop-bounded same-frame segment
cand='ATGAAAACC' # MKT
unit='CCCTAAGGG'+ 'ATGAAAACC' + 'GGGCCC' + 'TAA' + 'CCC'
# ensure a candidate exact lookup helper can find only one hit via temp fasta
with tempfile.TemporaryDirectory() as td:
 p=Path(td)/'u.fa'; p.write_text('>u1\n'+unit+'\n')
 rec,n=s.recover_maximal(p,cand)
 assert n==1 and rec is not None
 assert rec['right_stop_bounded'] is True
print('PASS_STAGE04B6_SELFTEST')
