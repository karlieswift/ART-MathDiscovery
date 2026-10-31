#!/usr/bin/env python3
import csv,json,os,subprocess,sys,time
from pathlib import Path

def find_cif(d):
    xs=sorted(Path(d).rglob('*_model_0.cif'))
    if not xs: xs=sorted(Path(d).rglob('*.cif'))
    return str(xs[0]) if xs else ''

def classify(log):
    t=Path(log).read_text(errors='ignore').lower() if Path(log).exists() else ''
    if 'out of memory' in t or 'cuda oom' in t: return 'CUDA_OOM'
    if 'msa' in t and ('error' in t or 'failed' in t): return 'MSA_OR_SERVER_FAILURE'
    if 'timed out' in t or 'timeout' in t: return 'TIMEOUT'
    return 'OTHER_TECHNICAL_FAILURE'

def main():
    if len(sys.argv)!=7: raise SystemExit('usage: run_d_boltz.py TARGETS OUTROOT LOGROOT STATEDIR STATUS GPU_ID')
    targets_p,outroot,logroot,statedir,status_p,gpu_id=sys.argv[1:]
    rows=list(csv.DictReader(open(targets_p,encoding='utf-8'),delimiter='\t'))
    if len(rows)!=6: raise SystemExit(f'FAIL structure target n={len(rows)}')
    outroot=Path(outroot);logroot=Path(logroot);statedir=Path(statedir);status_p=Path(status_p)
    for d in (outroot,logroot,statedir): d.mkdir(parents=True,exist_ok=True)
    recs=[]
    for r in rows:
        cid=r['candidate_id']; state=statedir/f'{cid}.json'; od=outroot/cid; log=logroot/f'd_boltz_{cid}.log'
        if state.exists():
            z=json.loads(state.read_text()); recs.append(z); print(f'[RESUME] {cid} {z["boltz_status"]}',flush=True); continue
        if od.exists():
            import shutil; shutil.rmtree(od)
        od.mkdir(parents=True,exist_ok=True)
        env=os.environ.copy(); env['CUDA_VISIBLE_DEVICES']=str(gpu_id)
        cmd=['boltz','predict',r['yaml_path'],'--out_dir',str(od),'--model','boltz2','--no_kernels','--use_msa_server','--seed','20261007','--recycling_steps','3','--sampling_steps','200','--diffusion_samples','1','--devices','1','--accelerator','gpu','--output_format','mmcif']
        t0=time.time()
        with log.open('w') as lf:
            lf.write('ROLE: post-discovery seven-family architecture-driven structure extension\n')
            lf.write('PHYSICAL_GPU_ID: '+str(gpu_id)+'\nCOMMAND: '+' '.join(cmd)+'\n');lf.flush()
            p=subprocess.run(cmd,stdout=lf,stderr=subprocess.STDOUT,env=env)
        elapsed=time.time()-t0; cif=find_cif(od)
        st='PASS' if p.returncode==0 and cif else 'TECHNICAL_FAILURE'
        z={'extension_rank':r['extension_rank'],'family_rep_id':r['family_rep_id'],'candidate_id':cid,'accession':r['accession'],'aa_length':r['aa_length'],'physical_gpu_id':str(gpu_id),'boltz_status':st,'return_code':str(p.returncode),'elapsed_seconds':f'{elapsed:.3f}','failure_category':'' if st=='PASS' else classify(log),'cif_path':cif}
        tmp=state.with_suffix('.tmp');tmp.write_text(json.dumps(z,indent=2)+'\n');tmp.replace(state);recs.append(z)
        print(f'[{st}] {cid} family={r["family_rep_id"]} aa={r["aa_length"]} elapsed={elapsed:.1f}s category={z["failure_category"]}',flush=True)
    recs.sort(key=lambda z:int(z['extension_rank']))
    cols=['extension_rank','family_rep_id','candidate_id','accession','aa_length','physical_gpu_id','boltz_status','return_code','elapsed_seconds','failure_category','cif_path']
    with status_p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(recs)
    print('target_n =',len(recs));print('boltz_pass_n =',sum(r['boltz_status']=='PASS' for r in recs));print('boltz_failure_n =',sum(r['boltz_status']!='PASS' for r in recs));print('PASS_F3B4D_D_BOLTZ_EXTENSION_ATTEMPT_COMPLETE')
if __name__=='__main__': main()
