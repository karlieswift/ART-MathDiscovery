#!/usr/bin/env python3
import csv, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

manifest=Path(sys.argv[1]); outroot=Path(sys.argv[2]); logroot=Path(sys.argv[3]); status_path=Path(sys.argv[4]); ngpu=int(sys.argv[5])
outroot.mkdir(parents=True, exist_ok=True); logroot.mkdir(parents=True, exist_ok=True)
with manifest.open() as f: rows=list(csv.DictReader(f, delimiter='\t'))
if len(rows)!=32: raise SystemExit('FAIL manifest not 32')

# deterministic round-robin by frozen shortlist rank
shards=[[] for _ in range(ngpu)]
for i,r in enumerate(rows): shards[i%ngpu].append(r)

def find_cif(d):
    p=Path(d)
    xs=sorted(p.rglob('*_model_0.cif'))
    if not xs: xs=sorted(p.rglob('*.cif'))
    return str(xs[0]) if xs else ''

def worker(gpu, items):
    recs=[]
    for r in items:
        cid=r['candidate_id']; od=outroot/cid; od.mkdir(parents=True,exist_ok=True)
        cif=find_cif(od)
        if cif:
            recs.append((r,'PASS_EXISTING',0,cif,0.0)); continue
        log=logroot/f'boltz_{cid}.log'
        env=os.environ.copy(); env['CUDA_VISIBLE_DEVICES']=str(gpu)
        cmd=[
            'boltz','predict',r['yaml_path'],
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
            lf.write('COMMAND: '+' '.join(cmd)+'\n')
            lf.flush()
            p=subprocess.run(cmd,stdout=lf,stderr=subprocess.STDOUT,env=env)
        elapsed=time.time()-t
        cif=find_cif(od)
        status='PASS' if p.returncode==0 and cif else 'TECHNICAL_FAILURE'
        recs.append((r,status,p.returncode,cif,elapsed))
    return recs

allrecs=[]
with ThreadPoolExecutor(max_workers=ngpu) as ex:
    futs=[ex.submit(worker,gpu,items) for gpu,items in enumerate(shards)]
    for fut in futs: allrecs.extend(fut.result())
allrecs.sort(key=lambda x:int(x[0]['shortlist_rank']))
cols=['shortlist_rank','family_rep_id','candidate_id','accession','aa_length','gpu','boltz_status','return_code','elapsed_seconds','cif_path']
# GPU is deterministic rank-1 mod ngpu
with status_path.open('w',newline='') as f:
    w=csv.DictWriter(f,delimiter='\t',fieldnames=cols); w.writeheader()
    for r,status,rc,cif,elapsed in allrecs:
        w.writerow({
            'shortlist_rank':r['shortlist_rank'],'family_rep_id':r['family_rep_id'],'candidate_id':r['candidate_id'],
            'accession':r['accession'],'aa_length':r['aa_length'],'gpu':(int(r['shortlist_rank'])-1)%ngpu,
            'boltz_status':status,'return_code':rc,'elapsed_seconds':f'{elapsed:.3f}','cif_path':cif
        })
print('query_n =',len(allrecs))
print('boltz_pass_n =',sum(x[1].startswith('PASS') for x in allrecs))
print('technical_failure_n =',sum(x[1]=='TECHNICAL_FAILURE' for x in allrecs))
print('PASS_F3B3C_BOLTZ_PRIMARY_ATTEMPT')
