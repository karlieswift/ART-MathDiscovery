#!/usr/bin/env python3
import csv, json, tempfile, subprocess, sys
from pathlib import Path
HERE=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader(open(HERE/'frozen/DELTA32_STRUCTURAL_COMPATIBILITY.tsv'),delimiter='\t'))
pos=[r for r in rows if r['tier'] in ('STRONG','SUGGESTIVE')]
assert len(rows)==32, len(rows)
assert len(pos)==12, len(pos)
assert sum(r['tier']=='STRONG' for r in pos)==3
assert sum(r['tier']=='SUGGESTIVE' for r in pos)==9
contract=json.load(open(HERE/'config/CONTRACT.json'))
assert contract['expected_structural_positive_n']==12
# import and synthetic repeat check
sys.path.insert(0,str(HERE/'scripts'))
import stage04b4_closeout as m
seq='A'*200 + 'ACGTACGT' + 'C'*150 + 'ACGTACGT' + 'G'*150 + 'ACGTACGT' + 'T'*300
arr=m.repeat_arrays(seq,900,1000)
assert any(x['copy_n']>=3 for x in arr), arr
print('PASS_STAGE04B4_SELFTEST')
