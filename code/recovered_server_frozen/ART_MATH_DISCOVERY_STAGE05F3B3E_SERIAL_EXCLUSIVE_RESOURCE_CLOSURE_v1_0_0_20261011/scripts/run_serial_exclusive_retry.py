#!/usr/bin/env python3
import csv, json, os, shutil, subprocess, sys, time
from pathlib import Path

manifest_path=Path(sys.argv[1])
target_path=Path(sys.argv[2])
outroot=Path(sys.argv[3])
logroot=Path(sys.argv[4])
state_root=Path(sys.argv[5])
status_path=Path(sys.argv[6])
gpu_id=int(sys.argv[7])

outroot.mkdir(parents=True,exist_ok=True)
logroot.mkdir(parents=True,exist_ok=True)
state_root.mkdir(parents=True,exist_ok=True)

with manifest_path.open(newline='') as f:
    manifest={r['candidate_id']:r for r in csv.DictReader(f,delimiter='\t')}
with target_path.open(newline='') as f:
    targets=list(csv.DictReader(f,delimiter='\t'))

assert len(targets)==7, len(targets)
assert len({r['candidate_id'] for r in targets})==7
for r in targets:
    assert r['candidate_id'] in manifest, r['candidate_id']
targets.sort(key=lambda r:int(r['shortlist_rank']))

def find_cif(d):
    p=Path(d)
    xs=sorted(p.rglob('*_model_0.cif'))
    if not xs: xs=sorted(p.rglob('*.cif'))
    return str(xs[0]) if xs else ''

def classify_failure(log_path):
    try:
        s=Path(log_path).read_text(errors='replace').lower()
    except Exception:
        return 'UNCLASSIFIED_RUNTIME_FAILURE'
    oom_tokens=['out of memory','cuda oom','torch.cuda.outofmemoryerror','cuda out of memory']
    if any(x in s for x in oom_tokens): return 'CUDA_OOM'
    if 'msa' in s and any(x in s for x in ['timeout','server','connection','network','request failed']): return 'MSA_SERVER_OR_NETWORK'
    if 'cuda' in s: return 'CUDA_RUNTIME'
    if any(x in s for x in ['validation','invalid','yaml','input']): return 'INPUT_OR_VALIDATION'
    return 'UNCLASSIFIED_RUNTIME_FAILURE'

records=[]
for t in targets:
    cid=t['candidate_id']; m=manifest[cid]
    state=state_root/f'{cid}.json'
    od=outroot/cid
    log=logroot/f'serial_boltz_{cid}.log'
    if state.exists():
        rec=json.loads(state.read_text())
        records.append(rec)
        print(f"[RESUME] {cid} {rec['serial_boltz_status']}",flush=True)
        continue

    # New stage output only. No primary/retry output is modified or reused as a result.
    if od.exists(): shutil.rmtree(od)
    od.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy()
    env['CUDA_VISIBLE_DEVICES']=str(gpu_id)
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
    t0=time.time()
    with log.open('w') as lf:
        lf.write('RESOURCE_MODE: serial exclusive single-GPU technical closure\n')
        lf.write('PHYSICAL_GPU_ID: '+str(gpu_id)+'\n')
        lf.write('COMMAND: '+' '.join(cmd)+'\n')
        lf.flush()
        p=subprocess.run(cmd,stdout=lf,stderr=subprocess.STDOUT,env=env)
    elapsed=time.time()-t0
    cif=find_cif(od)
    status='PASS_SERIAL_EXCLUSIVE' if p.returncode==0 and cif else 'TECHNICAL_FAILURE_SERIAL_EXCLUSIVE'
    category='' if status.startswith('PASS') else classify_failure(log)
    rec={
      'shortlist_rank':t['shortlist_rank'],'family_rep_id':t['family_rep_id'],'candidate_id':cid,
      'accession':t['accession'],'aa_length':t['aa_length'],'physical_gpu_id':str(gpu_id),
      'serial_boltz_status':status,'return_code':str(p.returncode),'elapsed_seconds':f'{elapsed:.3f}',
      'failure_category':category,'cif_path':cif
    }
    tmp=state.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(rec,indent=2)+'\n')
    tmp.replace(state)
    records.append(rec)
    print(f"[{status}] {cid} aa={t['aa_length']} elapsed={elapsed:.1f}s category={category}",flush=True)

records.sort(key=lambda r:int(r['shortlist_rank']))
cols=['shortlist_rank','family_rep_id','candidate_id','accession','aa_length','physical_gpu_id','serial_boltz_status','return_code','elapsed_seconds','failure_category','cif_path']
with status_path.open('w',newline='') as f:
    w=csv.DictWriter(f,delimiter='\t',fieldnames=cols,lineterminator='\n'); w.writeheader(); w.writerows(records)
print('serial_target_n =',len(records))
print('serial_boltz_pass_n =',sum(r['serial_boltz_status'].startswith('PASS') for r in records))
print('serial_boltz_failure_n =',sum(r['serial_boltz_status'].startswith('TECHNICAL_FAILURE') for r in records))
print('PASS_F3B3E_SERIAL_EXCLUSIVE_BOLTZ_ATTEMPT_COMPLETE')
