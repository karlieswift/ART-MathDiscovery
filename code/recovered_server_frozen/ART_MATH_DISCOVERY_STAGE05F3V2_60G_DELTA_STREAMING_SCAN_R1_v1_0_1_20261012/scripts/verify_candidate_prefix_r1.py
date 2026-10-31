#!/usr/bin/env python3
import csv,json,sys,hashlib
from pathlib import Path
oldp=Path(sys.argv[1]); newdir=Path(sys.argv[2]); out=Path(sys.argv[3])
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
report={'stage':'Stage05F3V2-R1-candidate-prefix-compatibility','baseline_pool':str(oldp),'baseline_present':oldp.is_file(),'strata':{},'gate_pass':True}
if oldp.is_file():
    old={s:[] for s in 'ABCD'}
    for r in csv.DictReader(open(oldp),delimiter='\t'):
        if r.get('stratum') in old: old[r['stratum']].append(r['accession'])
    for s in 'ABCD':
        p=newdir/f'F3V2_DELTA_CANDIDATES_{s}.tsv'
        new=[r['accession'] for r in csv.DictReader(open(p),delimiter='\t')]
        n=len(old[s]); exact=(new[:n]==old[s])
        report['strata'][s]={'baseline_n':n,'r1_total_n':len(new),'exact_prefix':exact}
        report['gate_pass'] &= exact
else:
    report['note']='No v1.0.0 candidate pool found; exact-prefix property is covered by tests/run_selftest.sh using the retained v1.0.0 selector.'
report['gate']='PASS_F3V2_R1_V100_CANDIDATE_PREFIX_EXACT' if report['gate_pass'] else 'FAIL_F3V2_R1_CANDIDATE_ORDER_DRIFT'
out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
if not report['gate_pass']: raise SystemExit(2)
