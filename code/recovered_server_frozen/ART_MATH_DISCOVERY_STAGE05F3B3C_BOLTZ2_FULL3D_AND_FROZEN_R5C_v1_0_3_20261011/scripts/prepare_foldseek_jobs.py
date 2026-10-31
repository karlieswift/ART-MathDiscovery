#!/usr/bin/env python3
import csv, json, sys
from pathlib import Path

status=Path(sys.argv[1])
out=Path(sys.argv[2])
with status.open(newline='') as f:
    rows=list(csv.DictReader(f, delimiter='\t'))
assert len(rows)==32, len(rows)
pass_rows=[]
tech=[]
for r in rows:
    cid=r['candidate_id'].strip()
    bst=r['boltz_status'].strip()
    cif=r['cif_path'].strip().rstrip('\r')
    if bst.startswith('PASS'):
        p=Path(cif)
        if not p.is_file() or p.stat().st_size==0:
            raise SystemExit(f'FAIL missing/nonempty CIF for Boltz PASS candidate {cid}: {cif!r}')
        rr=dict(r)
        rr['candidate_id']=cid
        rr['boltz_status']=bst
        rr['cif_path']=str(p)
        pass_rows.append(rr)
    else:
        tech.append(cid)
cols=['shortlist_rank','family_rep_id','candidate_id','accession','aa_length','boltz_status','cif_path']
with out.open('w', newline='\n') as f:
    w=csv.DictWriter(f, delimiter='\t', fieldnames=cols, lineterminator='\n')
    w.writeheader()
    for r in pass_rows:
        w.writerow({k:r[k] for k in cols})
print(json.dumps({
    'structure_handoff_n':len(rows),
    'boltz_pass_n':len(pass_rows),
    'boltz_technical_failure_n':len(tech),
    'lf_sanitized_job_file':str(out)
}, indent=2))
print('PASS_F3B3C_V103_FOLDSEEK_JOB_HANDOFF_PREPARED')
