#!/usr/bin/env python3
import csv,json,sys,os,tempfile,subprocess
r3=sys.argv[1]
g=json.load(open(os.path.join(r3,'outputs','GATE_STAGE05D4R3.json')))
assert g['gate_pass'] and g['n_tier_A']==350
p=os.path.join(r3,'outputs','STAGE05D4R3_NR90_STRUCTURE_PRIORITY.tsv')
with open(p) as f:
    r=csv.DictReader(f,delimiter='\t'); n=sum(x['structure_priority_tier']=='TIER_A_COMPACT_STRUCTURE_READY' for x in r)
assert n==350,n
fa=os.path.join(r3,'outputs','STAGE05D4R3_TIER_A_COMPACT_STRUCTURE_READY.faa')
assert sum(1 for x in open(fa) if x.startswith('>'))==350
print('[PASS] Stage05D4-R4 selftest: frozen Tier-A=350; no DNA rescoring')
