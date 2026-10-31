#!/usr/bin/env python3
import argparse,csv,hashlib,json,os,re,shutil,subprocess,sys,tempfile
from collections import defaultdict
from pathlib import Path

STRATA=("virome","marine","freshwater","soil","wastewater","human","livestock","other_metagenome")

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def read_fasta_stream(proc):
    name=None; seq=[]
    for raw in proc.stdout:
        line=raw.decode('ascii','ignore').strip()
        if not line: continue
        if line.startswith('>'):
            if name is not None: yield name,''.join(seq)
            name=line[1:].split()[0];seq=[]
        else: seq.append(line.upper())
    if name is not None: yield name,''.join(seq)

def get_uri(acc): return f's3://logan-pub/c/{acc}/{acc}.contigs.fa.zst'
def http_uri(acc): return f'https://s3.amazonaws.com/logan-pub/c/{acc}/{acc}.contigs.fa.zst'

def download(acc,dst):
    if dst.exists() and dst.stat().st_size>0:return True,'cache'
    dst.parent.mkdir(parents=True,exist_ok=True)
    if shutil.which('aws'):
        cmd=['aws','s3','cp',get_uri(acc),str(dst),'--no-sign-request','--only-show-errors']
    elif shutil.which('wget'):
        cmd=['wget','-q','-O',str(dst),http_uri(acc)]
    else: raise RuntimeError('need aws CLI or wget')
    rc=subprocess.run(cmd).returncode
    if rc!=0 or not dst.exists() or dst.stat().st_size==0:
        try: dst.unlink()
        except: pass
        return False,'download_failed'
    return True,'downloaded'

def count_eligible(zst,min_nt):
    p=subprocess.Popen(['zstdcat',str(zst)],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
    n=bp=0
    for _,seq in read_fasta_stream(p):
        if len(seq)>=min_nt:
            n+=1;bp+=len(seq)
    rc=p.wait()
    if rc!=0: raise RuntimeError(f'zstdcat failed {zst}')
    return n,bp

def append_eligible(zst,acc,min_nt,out):
    p=subprocess.Popen(['zstdcat',str(zst)],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
    n=bp=0
    for name,seq in read_fasta_stream(p):
        if len(seq)<min_nt: continue
        safe=re.sub(r'[^A-Za-z0-9_.:-]+','_',name)
        out.write(f'>{acc}|{safe}\n')
        for i in range(0,len(seq),80):out.write(seq[i:i+80]+'\n')
        n+=1;bp+=len(seq)
    rc=p.wait()
    if rc!=0: raise RuntimeError(f'zstdcat failed second pass {zst}')
    return n,bp

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--candidates',required=True);ap.add_argument('--external-root',required=True);ap.add_argument('--out-fasta',required=True);ap.add_argument('--out-manifest',required=True);ap.add_argument('--out-gate',required=True)
    ap.add_argument('--target-bp',type=int,default=6000000000);ap.add_argument('--min-bp',type=int,default=5000000000);ap.add_argument('--max-bp',type=int,default=10000000000);ap.add_argument('--min-contig-nt',type=int,default=300);ap.add_argument('--max-accession-bp',type=int,default=250000000)
    a=ap.parse_args(); ext=Path(a.external_root);cache=ext/'downloads';cache.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader(open(a.candidates),delimiter='\t')); by=defaultdict(list)
    for r in rows:by[r['stratum']].append(r)
    for s in by:by[s].sort(key=lambda r:int(r['stable_rank']))
    idx={s:0 for s in STRATA}; bp_by=defaultdict(int); acc_by=defaultdict(int); accepted=[]; failed=[]; total=0
    outp=Path(a.out_fasta);outp.parent.mkdir(parents=True,exist_ok=True)
    with outp.open('w') as out:
        exhausted=False
        while total<a.target_bp and not exhausted:
            progress=False
            for st in STRATA:
                if total>=a.target_bp: break
                q=by.get(st,[])
                while idx[st]<len(q):
                    r=q[idx[st]];idx[st]+=1;acc=r['accession'];z=cache/f'{acc}.contigs.fa.zst'
                    ok,route=download(acc,z)
                    if not ok:
                        failed.append({'accession':acc,'stratum':st,'reason':route});continue
                    zsha=sha256_file(z); zbytes=z.stat().st_size
                    try:n,bp=count_eligible(z,a.min_contig_nt)
                    except Exception as e:
                        failed.append({'accession':acc,'stratum':st,'reason':f'parse_failed:{e}'});continue
                    if bp==0:
                        failed.append({'accession':acc,'stratum':st,'reason':'no_eligible_contigs'});continue
                    if bp>a.max_accession_bp:
                        failed.append({'accession':acc,'stratum':st,'reason':f'accession_bp_over_cap:{bp}'});continue
                    if total+bp>a.max_bp:
                        failed.append({'accession':acc,'stratum':st,'reason':f'would_exceed_global_max:{bp}'});continue
                    n2,bp2=append_eligible(z,acc,a.min_contig_nt,out)
                    if (n2,bp2)!=(n,bp):raise RuntimeError('two-pass count mismatch')
                    rec=dict(r);rec.update({'source_uri':get_uri(acc),'compressed_bytes':zbytes,'compressed_sha256':zsha,'retained_contig_n':n,'retained_bp':bp})
                    accepted.append(rec);total+=bp;bp_by[st]+=bp;acc_by[st]+=1;progress=True
                    # raw compressed file can be re-fetched from immutable URI; retain only if env requests it
                    if os.environ.get('F3_KEEP_LOGAN_ZST','0')!='1':
                        try:z.unlink()
                        except:pass
                    break
            if not progress:exhausted=True
    strata100=sum(1 for s in STRATA if bp_by[s]>=100_000_000)
    gate_pass=(a.min_bp<=total<=a.max_bp and strata100>=5 and len(accepted)>=24)
    cols=['stratum','stable_rank','accession','assay_type','organism','scientific_name','stable_sha256','source_uri','compressed_bytes','compressed_sha256','retained_contig_n','retained_bp']
    with open(a.out_manifest,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(accepted)
    failp=Path(a.out_manifest).with_name('STAGE05F3A_SOURCE_SKIPS.tsv')
    with failp.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['accession','stratum','reason'],delimiter='\t');w.writeheader();w.writerows(failed)
    g={'stage':'Stage05F3A0','gate':'PASS_F3A0_EXTERNAL_SOURCE_UNIVERSE_FROZEN' if gate_pass else 'HOLD_F3A0_SOURCE_UNIVERSE_INCOMPLETE','gate_pass':gate_pass,'accepted_accession_n':len(accepted),'total_bp':total,'total_contig_n':sum(int(x['retained_contig_n']) for x in accepted),'strata_bp':dict(bp_by),'strata_accession_n':dict(acc_by),'strata_ge_100mb_n':strata100,'target_bp':a.target_bp,'min_bp':a.min_bp,'max_bp':a.max_bp,'min_contig_nt':a.min_contig_nt,'max_accession_retained_bp':a.max_accession_bp,'external_fasta_sha256':sha256_file(outp),'source_manifest_sha256':sha256_file(Path(a.out_manifest)),'source_selection_used_query_sequences':False,'source_selection_used_rt_results':False}
    Path(a.out_gate).write_text(json.dumps(g,indent=2)+'\n')
    print(json.dumps(g,indent=2))
    if not gate_pass:raise SystemExit(3)
if __name__=='__main__':main()
