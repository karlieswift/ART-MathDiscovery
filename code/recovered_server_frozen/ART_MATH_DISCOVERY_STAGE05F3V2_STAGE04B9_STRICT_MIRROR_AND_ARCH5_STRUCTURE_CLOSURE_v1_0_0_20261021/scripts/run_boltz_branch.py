#!/usr/bin/env python3
import argparse,csv,os,shutil,subprocess,time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

def rows(p):
 with open(p,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def find_cif(d,q):
 x=Path(d)/f'{q}.cif'
 if x.is_file():return str(x)
 xs=sorted(Path(d).rglob(f'{q}_model_0.cif'))
 if not xs:xs=sorted(Path(d).rglob('*_model_0.cif'))
 return str(xs[0]) if xs else ''
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--manifest',required=True);ap.add_argument('--structures-dir',required=True);ap.add_argument('--work-dir',required=True);ap.add_argument('--logs-dir',required=True);ap.add_argument('--gpus',required=True);ap.add_argument('--status',required=True);a=ap.parse_args()
 m=rows(a.manifest);gpus=[x.strip() for x in a.gpus.split(',') if x.strip()];
 if not gpus:raise SystemExit('HOLD no GPUs')
 SD=Path(a.structures_dir);WD=Path(a.work_dir);LD=Path(a.logs_dir);SD.mkdir(parents=True,exist_ok=True);WD.mkdir(parents=True,exist_ok=True);LD.mkdir(parents=True,exist_ok=True)
 todo=[];done=[]
 for r in m:
  q=r['query_id'];existing=SD/f'{q}.cif'
  if existing.is_file():done.append((r,'PASS_EXISTING',0,str(existing),0.0));continue
  todo.append(r)
 shards=[[] for _ in gpus]
 for i,r in enumerate(todo):shards[i%len(gpus)].append(r)
 def worker(gpu,items):
  out=[]
  for r in items:
   q=r['query_id'];od=WD/q;od.mkdir(parents=True,exist_ok=True);log=LD/f'boltz_{q}.log';env=os.environ.copy();env['CUDA_VISIBLE_DEVICES']=gpu
   cmd=['boltz','predict',r['yaml_path'],'--out_dir',str(od),'--model','boltz2','--no_kernels','--use_msa_server','--seed','20261007','--recycling_steps','3','--sampling_steps','200','--diffusion_samples','1','--devices','1','--accelerator','gpu','--output_format','mmcif']
   t=time.time()
   with open(log,'w') as lf:
    lf.write('COMMAND: '+' '.join(cmd)+'\n');lf.flush();p=subprocess.run(cmd,stdout=lf,stderr=subprocess.STDOUT,env=env)
   xs=sorted(od.rglob('*_model_0.cif'))
   cif=str(xs[0]) if xs else ''
   if p.returncode==0 and cif:
    shutil.copy2(cif,SD/f'{q}.cif');status='PASS'
   else:status='TECHNICAL_FAILURE'
   out.append((r,status,p.returncode,cif,time.time()-t))
  return out
 rec=done[:]
 with ThreadPoolExecutor(max_workers=len(gpus)) as ex:
  for z in ex.map(lambda t:worker(*t),[(gpus[i],shards[i]) for i in range(len(gpus))]):rec.extend(z)
 rec.sort(key=lambda x:x[0]['query_id']);fields=['query_id','candidate_id','aa_length','boltz_status','return_code','elapsed_seconds','cif_path']
 with open(a.status,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
  for r,s,rc,cif,el in rec:w.writerow({'query_id':r['query_id'],'candidate_id':r['candidate_id'],'aa_length':r['aa_length'],'boltz_status':s,'return_code':rc,'elapsed_seconds':f'{el:.3f}','cif_path':str(SD/f"{r['query_id']}.cif") if (SD/f"{r['query_id']}.cif").is_file() else cif})
 n=sum((SD/f"{r['query_id']}.cif").is_file() for r in m);print(f'BOLTZ_AVAILABLE={n}/{len(m)}')
 if n!=len(m):raise SystemExit('HOLD_BOLTZ_INCOMPLETE')
 print('PASS_BOLTZ_BRANCH_COMPLETE')
if __name__=='__main__':main()
