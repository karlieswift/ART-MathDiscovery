#!/usr/bin/env python3
import csv, hashlib, json, re, sys
from pathlib import Path

status_path=Path(sys.argv[1])
logdir=Path(sys.argv[2])
primary_summary=Path(sys.argv[3])
outdir=Path(sys.argv[4])
outdir.mkdir(parents=True,exist_ok=True)

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

with status_path.open(newline='') as f:
    status=list(csv.DictReader(f,delimiter='\t'))
with primary_summary.open(newline='') as f:
    summ={r['candidate_id']:r for r in csv.DictReader(f,delimiter='\t')}

fails=[r for r in status if r['boltz_status'].strip()=='TECHNICAL_FAILURE']
assert len(status)==32, len(status)
assert len(fails)==10, len(fails)

rules=[
 ('CUDA_OOM', re.compile(r'out of memory|torch\.OutOfMemoryError|CUDA error:\s*out of memory|CUBLAS_STATUS_ALLOC_FAILED',re.I)),
 ('MSA_SERVER_OR_NETWORK', re.compile(r'429|too many requests|502|503|504|connectionerror|connecttimeout|readtimeout|remotedisconnected|name.?resolution|failed to establish a new connection|temporary failure|requests\.exceptions|httperror|server disconnected|msa.{0,40}server|mmseqs.{0,40}server|colabfold.{0,40}server',re.I|re.S)),
 ('CUDA_RUNTIME', re.compile(r'CUDA error|cuDNN|CUDNN|CUBLAS|device-side assert|illegal memory access',re.I)),
 ('INPUT_OR_VALIDATION', re.compile(r'validationerror|invalid yaml|yaml.*error|parse.*error|invalid sequence|unsupported residue',re.I)),
 ('MODEL_NUMERICAL_OR_RUNTIME', re.compile(r'nan|inf detected|runtimeerror|assertionerror',re.I)),
]

out=[]
for r in fails:
    cid=r['candidate_id'].strip()
    lp=logdir/f'boltz_{cid}.log'
    text=lp.read_text(errors='replace') if lp.exists() else ''
    cause='UNCLASSIFIED_RUNTIME_FAILURE'
    matched=''
    for name,rx in rules:
        m=rx.search(text)
        if m:
            cause=name
            matched=m.group(0).replace('\t',' ')[:160]
            break
    nonempty=[x.strip() for x in text.splitlines() if x.strip()]
    tail=' || '.join(nonempty[-12:])[:5000]
    s=summ.get(cid,{})
    out.append({
      'shortlist_rank':r['shortlist_rank'],
      'family_rep_id':r['family_rep_id'],
      'candidate_id':cid,
      'accession':r['accession'],
      'aa_length':r['aa_length'],
      'primary_return_code':r['return_code'],
      'primary_elapsed_seconds':r['elapsed_seconds'],
      'primary_best_tier':s.get('best_tier',''),
      'log_exists':str(lp.exists()),
      'log_sha256':sha(lp) if lp.exists() else '',
      'diagnostic_category':cause,
      'matched_pattern':matched,
      'log_tail':tail,
    })

cols=list(out[0])
with (outdir/'F3B3D_PRIMARY_FAILURE_DIAGNOSTICS.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,delimiter='\t',fieldnames=cols,lineterminator='\n'); w.writeheader(); w.writerows(out)

counts={}
for r in out: counts[r['diagnostic_category']]=counts.get(r['diagnostic_category'],0)+1
meta={
 'stage':'Stage05F3B3D-diagnostic',
 'primary_technical_failure_n':len(out),
 'diagnostic_category_counts':counts,
 'diagnostic_role':'descriptive only; categories do not control retry eligibility',
 'all_primary_technical_failures_retry_eligible_under_frozen_exact-protocol_rule':True,
}
with (outdir/'F3B3D_PRIMARY_FAILURE_DIAGNOSTICS.json').open('w') as f:
    json.dump(meta,f,indent=2); f.write('\n')
print(json.dumps(meta,indent=2))
print('PASS_F3B3D_PRIMARY_FAILURE_DIAGNOSIS_COMPLETE')
