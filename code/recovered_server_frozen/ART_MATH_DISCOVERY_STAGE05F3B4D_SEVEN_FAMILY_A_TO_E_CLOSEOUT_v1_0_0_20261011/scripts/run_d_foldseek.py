#!/usr/bin/env python3
import csv,subprocess,sys,shutil
from pathlib import Path

def main():
    if len(sys.argv)!=7: raise SystemExit('usage: run_d_foldseek.py BOLTZ_STATUS RTDB OUTDIR TMPDIR LOGDIR THREADS')
    boltz_p,rtdb,outdir,tmpdir,logdir,threads=sys.argv[1:]
    outdir=Path(outdir);tmpdir=Path(tmpdir);logdir=Path(logdir)
    for d in (outdir,tmpdir,logdir):d.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader(open(boltz_p,encoding='utf-8'),delimiter='\t'))
    rec=[]
    for r in rows:
        cid=r['candidate_id']
        if r['boltz_status']!='PASS':
            rec.append({'candidate_id':cid,'foldseek_status':'NOT_RUN_BOLTZ_FAILURE','return_code':'','m8_path':''});continue
        cif=r['cif_path'];m8=outdir/f'{cid}.m8';tmp=tmpdir/cid;log=logdir/f'd_foldseek_{cid}.log'
        if tmp.exists():shutil.rmtree(tmp)
        if m8.exists():m8.unlink()
        cmd=['foldseek','easy-search',cif,rtdb,str(m8),str(tmp),'-e','1.0','--max-seqs','100','--threads',str(threads),'--exact-tmscore','1','--format-output','query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore']
        with log.open('w') as lf:p=subprocess.run(cmd,stdout=lf,stderr=subprocess.STDOUT)
        st='PASS' if p.returncode==0 and m8.exists() else 'TECHNICAL_FAILURE'
        rec.append({'candidate_id':cid,'foldseek_status':st,'return_code':p.returncode,'m8_path':str(m8)})
        print(f'[{st}] {cid} rc={p.returncode}',flush=True)
    sp=Path(outdir).parent.parent/'outputs/F3B4D_D_FOLDSEEK_STATUS.tsv';sp.parent.mkdir(parents=True,exist_ok=True)
    cols=['candidate_id','foldseek_status','return_code','m8_path']
    with sp.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rec)
    print('foldseek_pass_n =',sum(r['foldseek_status']=='PASS' for r in rec));print('PASS_F3B4D_D_FOLDSEEK_RUN_COMPLETE')
if __name__=='__main__':main()
