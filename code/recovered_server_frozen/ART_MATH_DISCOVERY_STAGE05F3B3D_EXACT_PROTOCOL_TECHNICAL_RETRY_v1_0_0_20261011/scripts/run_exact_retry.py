#!/usr/bin/env python3
import csv, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

manifest_path=Path(sys.argv[1])
primary_status_path=Path(sys.argv[2])
outroot=Path(sys.argv[3])
logroot=Path(sys.argv[4])
status_path=Path(sys.argv[5])
ngpu=int(sys.argv[6])
outroot.mkdir(parents=True,exist_ok=True); logroot.mkdir(parents=True,exist_ok=True)

with manifest_path.open(newline='') as f:
    manifest={r['candidate_id']:r for r in csv.DictReader(f,delimiter='\t')}
with primary_status_path.open(newline='') as f:
    primary=list(csv.DictReader(f,delimiter='\t'))

targets=[r for r in primary if r['boltz_status'].strip()=='TECHNICAL_FAILURE']
assert len(primary)==32, len(primary)
assert len(targets)==10, len(targets)
for r in targets:
    assert r['candidate_id'] in manifest, r['candidate_id']

# Frozen rule: retry all and only the 10 primary technical failures exactly once.
# No length/cause/content-based filtering is allowed.
targets=sorted(targets,key=lambda r:int(r['shortlist_rank']))
shards=[[] for _ in range(ngpu)]
for i,r in enumerate(targets): shards[i%ngpu].append(r)

def find_cif(d):
    p=Path(d)
    xs=sorted(p.rglob('*_model_0.cif'))
    if not xs: xs=sorted(p.rglob('*.cif'))
    return str(xs[0]) if xs else ''

def worker(gpu,items):
    recs=[]
    for p in items:
        cid=p['candidate_id']; m=manifest[cid]
        od=outroot/cid; od.mkdir(parents=True,exist_ok=True)
        cif=find_cif(od)
        if cif:
            recs.append((p,'PASS_EXISTING_RETRY',0,cif,0.0)); continue
        log=logroot/f'retry_boltz_{cid}.log'
        env=os.environ.copy(); env['CUDA_VISIBLE_DEVICES']=str(gpu)
        cmd=[
          'boltz','predict',m['yaml_path'],
          '--out_dir',str(od),
          '--model','boltz2',
          '--no_kernels',
          '--use_msa_server',
          '--seed','20261007',
          '--recycling_steps','3',
          '--sampling_steps','200',
          '--diffusion_samples','1',
          '--devices','1',
          '--accelerator','gpu',
          '--output_format','mmcif'
        ]
        t=time.time()
        with log.open('w') as lf:
            lf.write('COMMAND: '+' '.join(cmd)+'\n'); lf.flush()
            proc=subprocess.run(cmd,stdout=lf,stderr=subprocess.STDOUT,env=env)
        elapsed=time.time()-t
        cif=find_cif(od)
        status='PASS_RETRY' if proc.returncode==0 and cif else 'TECHNICAL_FAILURE_RETRY'
        recs.append((p,status,proc.returncode,cif,elapsed))
    return recs

allrecs=[]
with ThreadPoolExecutor(max_workers=ngpu) as ex:
    futs=[ex.submit(worker,gpu,items) for gpu,items in enumerate(shards)]
    for fut in futs: allrecs.extend(fut.result())
allrecs.sort(key=lambda x:int(x[0]['shortlist_rank']))
cols=['shortlist_rank','family_rep_id','candidate_id','accession','aa_length','gpu','retry_boltz_status','return_code','elapsed_seconds','cif_path']
with status_path.open('w',newline='') as f:
    w=csv.DictWriter(f,delimiter='\t',fieldnames=cols,lineterminator='\n'); w.writeheader()
    for r,status,rc,cif,elapsed in allrecs:
        w.writerow({
          'shortlist_rank':r['shortlist_rank'],'family_rep_id':r['family_rep_id'],'candidate_id':r['candidate_id'],
          'accession':r['accession'],'aa_length':r['aa_length'],'gpu':(int(r['shortlist_rank'])-1)%ngpu,
          'retry_boltz_status':status,'return_code':rc,'elapsed_seconds':f'{elapsed:.3f}','cif_path':cif
        })
print('retry_target_n =',len(allrecs))
print('retry_boltz_pass_n =',sum(x[1].startswith('PASS') for x in allrecs))
print('retry_boltz_failure_n =',sum(x[1]=='TECHNICAL_FAILURE_RETRY' for x in allrecs))
print('PASS_F3B3D_EXACT_PROTOCOL_BOLTZ_RETRY_COMPLETE')
