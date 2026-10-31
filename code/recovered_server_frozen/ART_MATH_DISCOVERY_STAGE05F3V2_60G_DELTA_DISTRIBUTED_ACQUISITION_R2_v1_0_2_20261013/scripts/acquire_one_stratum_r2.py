#!/usr/bin/env python3
import argparse,csv,fcntl,hashlib,json,os,re,socket,subprocess,sys,time
from pathlib import Path

TARGET={'A':13_500_000_000,'B':18_000_000_000,'C':13_500_000_000,'D':9_000_000_000}
BASE_FIELDS=['stratum','accession','nucleotide_bp','contig_n','compressed_bytes','zst_sha256','s3_uri']
FIELDS=['deterministic_rank']+BASE_FIELDS
LEDGER_FIELDS=['event_id','stratum','deterministic_rank','accession','status','compressed_bytes','nucleotide_bp','contig_n','detail']
PERMANENT={'LOGAN_OBJECT_NOT_FOUND','COMPRESSED_SIZE_OUT_OF_RANGE','EMPTY_OR_INVALID_FASTA','REALIZED_BP_OUT_OF_RANGE'}

def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def read_tsv(p):
    p=Path(p)
    if not p.is_file() or p.stat().st_size==0: return []
    with open(p,newline='') as f: return list(csv.DictReader(f,delimiter='\t'))

def init_tsv(p,fields):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    if not p.exists() or p.stat().st_size==0:
        with open(p,'w',newline='') as f: csv.DictWriter(f,delimiter='\t',fieldnames=fields).writeheader()

def append_tsv(p,fields,row):
    with open(p,'a',newline='') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=fields,extrasaction='ignore'); w.writerow(row); f.flush(); os.fsync(f.fileno())

def run(cmd): return subprocess.run(cmd,text=True,capture_output=True)

def count_fasta(zst,count_script):
    p1=subprocess.Popen(['zstd','-dc',str(zst)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    p2=subprocess.Popen([sys.executable,str(count_script)],stdin=p1.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    p1.stdout.close(); out2,err2=p2.communicate(); err1=p1.stderr.read().decode(errors='replace'); rc1=p1.wait()
    if rc1!=0 or p2.returncode!=0: raise RuntimeError((err1+' '+err2).strip()[:500])
    j=json.loads(out2); return int(j['bp']),int(j['records'])

def link_or_keep(src,dst):
    src=Path(src).resolve(); dst=Path(dst); dst.parent.mkdir(parents=True,exist_ok=True)
    if dst.exists() or dst.is_symlink():
        try:
            if os.path.samefile(src,dst): return 'existing-samefile'
        except Exception: pass
        if dst.is_symlink() or dst.is_file(): dst.unlink()
        else: raise RuntimeError(f'cannot replace non-file {dst}')
    try:
        os.link(src,dst); return 'hardlink'
    except OSError:
        try:
            os.symlink(src,dst); return 'symlink'
        except OSError:
            import shutil; shutil.copy2(src,dst); return 'copy'

def find_missing_ranks(candidate_tsv, wanted):
    wanted=set(wanted); out={}
    if not wanted: return out
    with open(candidate_tsv,newline='') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            a=r['accession']
            if a in wanted:
                out[a]=r['deterministic_rank']; wanted.remove(a)
                if not wanted: break
    if wanted: raise SystemExit(f'FATAL: cannot recover deterministic rank for imported accessions: {sorted(wanted)[:10]}')
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--stage',required=True)
    ap.add_argument('--source-r1-stage',required=True)
    ap.add_argument('--stratum',choices=list('ABCD'),required=True)
    ap.add_argument('--max-accession-bp',type=int,default=int(os.environ.get('MAX_ACCESSION_BP','1500000000')))
    ap.add_argument('--min-zst-bytes',type=int,default=int(os.environ.get('MIN_ZST_BYTES','10000000')))
    ap.add_argument('--max-zst-bytes',type=int,default=int(os.environ.get('MAX_ZST_BYTES','2000000000')))
    ap.add_argument('--transport-fail-limit',type=int,default=int(os.environ.get('F3V2_TRANSPORT_FAIL_LIMIT','5')))
    ap.add_argument('--download-retries',type=int,default=int(os.environ.get('F3V2_DOWNLOAD_RETRIES','3')))
    a=ap.parse_args()
    S=Path(a.stage).resolve(); R=Path(a.source_r1_stage).resolve(); STR=a.stratum
    cp=R/'inputs'/f'F3V2_DELTA_CANDIDATES_{STR}.tsv'; oldp=R/'inputs/OLD_F3_ACCESSIONS.txt'
    if not cp.is_file() or not oldp.is_file(): raise SystemExit('FATAL: source R1 stage is missing frozen candidate inputs; run R1 selector first')
    old={x.strip() for x in open(oldp) if x.strip() and not x.startswith('#')}
    shard=S/'work/shards'/STR; datadir=S/'data/shards'/STR; shard.mkdir(parents=True,exist_ok=True); datadir.mkdir(parents=True,exist_ok=True)
    acceptedp=shard/'accepted.partial.tsv'; ledgerp=shard/'ledger.tsv'; init_tsv(acceptedp,FIELDS); init_tsv(ledgerp,LEDGER_FIELDS)
    lockdir=S/'work/locks'; lockdir.mkdir(parents=True,exist_ok=True); lockpath=lockdir/f'{STR}.lock'
    lockfh=open(lockpath,'a+');
    try: fcntl.flock(lockfh.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError: raise SystemExit(f'FATAL: stratum {STR} is already leased by another R2 process: {lockpath}')
    lockfh.seek(0); lockfh.truncate(); lockfh.write(f'host={socket.gethostname()} pid={os.getpid()} started={time.time()}\n'); lockfh.flush()

    event_id=max([int(r.get('event_id') or 0) for r in read_tsv(ledgerp)] or [0])
    def log(row,status,detail='',compressed_bytes='',bp='',nr=''):
        nonlocal event_id
        event_id+=1
        append_tsv(ledgerp,LEDGER_FIELDS,{'event_id':event_id,'stratum':STR,'deterministic_rank':row.get('deterministic_rank',''),'accession':row.get('accession',''),'status':status,'compressed_bytes':compressed_bytes,'nucleotide_bp':bp,'contig_n':nr,'detail':detail[:500]})

    # Import current R1 partial acceptances for this stratum; this is why the serial R1 process must be stopped first.
    accepted=read_tsv(acceptedp); aset={r['accession'] for r in accepted}; imported=[]
    srcacc=R/'work/download/accepted.partial.tsv'
    if not srcacc.is_file(): srcacc=R/'frozen/F3V2_DELTA_SELECTED_ACCESSIONS.tsv'
    srcrows=[r for r in read_tsv(srcacc) if r.get('stratum')==STR and r.get('accession') not in aset]
    # Recover ranks from R1 ledger where possible.
    rank_by_acc={}
    for lr in read_tsv(R/'work/download/F3V2_ACQUISITION_LEDGER.tsv'):
        if lr.get('stratum')==STR and lr.get('accession') and lr.get('deterministic_rank'):
            rank_by_acc[lr['accession']]=lr['deterministic_rank']
    missing=[r['accession'] for r in srcrows if r['accession'] not in rank_by_acc]
    rank_by_acc.update(find_missing_ranks(cp,missing))
    for r in srcrows:
        acc=r['accession']
        if acc in old: raise SystemExit(f'FATAL: source R1 partial contains old-universe accession {acc}')
        src=R/'data/contigs'/f'{acc}.contigs.fa.zst'
        if not src.is_file() or src.stat().st_size==0: raise SystemExit(f'FATAL: source R1 accepted contig missing: {src}')
        actual=sha256(src); expected=r.get('zst_sha256','')
        if expected and expected!=actual: raise SystemExit(f'FATAL: source R1 SHA mismatch for {acc}')
        dst=datadir/f'{acc}.contigs.fa.zst'; mode=link_or_keep(src,dst)
        rec={'deterministic_rank':rank_by_acc[acc],**{k:r.get(k,'') for k in BASE_FIELDS}}
        append_tsv(acceptedp,FIELDS,rec); accepted.append(rec); aset.add(acc); imported.append(acc)
        log(rec,'ACCEPTED_IMPORTED_R1',f'SHA-verified serial R1 partial; {mode}',rec['compressed_bytes'],rec['nucleotide_bp'],rec['contig_n'])
    if imported: print(f'[R2 IMPORT] stratum={STR} imported_r1_accepted={len(imported)}',flush=True)

    # Reuse permanent decisions from the stopped serial R1 ledger and this shard ledger.
    last={}
    for lr in read_tsv(R/'work/download/F3V2_ACQUISITION_LEDGER.tsv'):
        if lr.get('stratum')==STR and lr.get('accession'): last[lr['accession']]=lr.get('status','')
    for lr in read_tsv(ledgerp):
        if lr.get('accession'): last[lr['accession']]=lr.get('status','')

    accepted=read_tsv(acceptedp); aset={r['accession'] for r in accepted}
    got=sum(int(r['nucleotide_bp']) for r in accepted if r['stratum']==STR)
    total=max(sum(1 for _ in open(cp))-1,0)
    print(f'[R2 STRATUM {STR}] target_bp={TARGET[STR]} deterministic_candidates={total} resume_bp={got}',flush=True)
    count_script=S/'scripts/count_fasta_stream.py'; transport_fail=0; local_seen=set()
    with open(cp,newline='') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            if got>=TARGET[STR]: break
            acc=row['accession']; rank=row['deterministic_rank']; row['stratum']=STR
            if acc in local_seen:
                log(row,'DUPLICATE_ACCESSION','duplicate accession inside deterministic candidate table'); continue
            local_seen.add(acc)
            if acc in old:
                log(row,'ALREADY_IN_OLD_UNIVERSE','fatal invariant violation'); raise SystemExit(f'FATAL old accession leaked into delta: {acc}')
            if acc in aset or last.get(acc) in PERMANENT: continue
            log(row,'SELECTED','candidate reached by deterministic acquisition stream')
            key=f'c/{acc}/{acc}.contigs.fa.zst'; uri=f's3://logan-pub/{key}'; z=datadir/f'{acc}.contigs.fa.zst'; bytes_n=0
            if not z.is_file() or z.stat().st_size==0:
                q=run(['aws','s3api','head-object','--bucket','logan-pub','--key',key,'--no-sign-request','--query','ContentLength','--output','text'])
                if q.returncode!=0:
                    msg=(q.stderr or q.stdout or '').strip().replace('\n',' ')
                    if re.search(r'(^|\D)404(\D|$)|NoSuchKey|Not Found',msg,re.I):
                        log(row,'LOGAN_OBJECT_NOT_FOUND',msg); last[acc]='LOGAN_OBJECT_NOT_FOUND'; transport_fail=0; continue
                    log(row,'HEAD_FAILED_TRANSIENT',msg); transport_fail+=1
                    print(f'[TRANSIENT HEAD] {STR} {acc} consecutive={transport_fail}/{a.transport_fail_limit}: {msg[:180]}',file=sys.stderr,flush=True)
                    if transport_fail>=a.transport_fail_limit: raise SystemExit('FATAL: repeated S3 transport/head failures; stop fail-closed and rerun same stratum later')
                    continue
                transport_fail=0
                try: bytes_n=int(q.stdout.strip())
                except Exception:
                    log(row,'HEAD_FAILED_TRANSIENT',f'unparseable ContentLength: {q.stdout!r}'); transport_fail+=1; continue
                if bytes_n<a.min_zst_bytes or bytes_n>a.max_zst_bytes:
                    log(row,'COMPRESSED_SIZE_OUT_OF_RANGE',f'frozen range [{a.min_zst_bytes},{a.max_zst_bytes}]',bytes_n); last[acc]='COMPRESSED_SIZE_OUT_OF_RANGE'; continue
                print(f'[DOWNLOAD] {STR} rank={rank} {acc} compressed={bytes_n}',flush=True)
                ok=False; msg=''
                for attempt in range(1,a.download_retries+1):
                    part=Path(str(z)+'.part'); part.unlink(missing_ok=True)
                    d=run(['aws','s3','cp',uri,str(part),'--no-sign-request','--only-show-errors'])
                    if d.returncode==0 and part.is_file() and part.stat().st_size>0:
                        os.replace(part,z); ok=True; break
                    msg=(d.stderr or d.stdout or '').strip().replace('\n',' '); part.unlink(missing_ok=True); time.sleep(min(attempt*2,5))
                if not ok:
                    log(row,'DOWNLOAD_FAILED',msg,bytes_n); transport_fail+=1
                    if transport_fail>=a.transport_fail_limit: raise SystemExit('FATAL: repeated S3 download failures; stop fail-closed and rerun same stratum later')
                    continue
                transport_fail=0
            else: bytes_n=z.stat().st_size
            try: bp,nr=count_fasta(z,count_script)
            except Exception as ex:
                z.unlink(missing_ok=True); log(row,'DECOMPRESS_FAILED',str(ex),bytes_n); continue
            if bp<=0 or nr<=0:
                z.unlink(missing_ok=True); log(row,'EMPTY_OR_INVALID_FASTA','zero bp or zero records',bytes_n,bp,nr); last[acc]='EMPTY_OR_INVALID_FASTA'; continue
            if bp>a.max_accession_bp:
                z.unlink(missing_ok=True); log(row,'REALIZED_BP_OUT_OF_RANGE',f'bp>{a.max_accession_bp}',bytes_n,bp,nr); last[acc]='REALIZED_BP_OUT_OF_RANGE'; continue
            sh=sha256(z)
            rec={'deterministic_rank':rank,'stratum':STR,'accession':acc,'nucleotide_bp':bp,'contig_n':nr,'compressed_bytes':bytes_n,'zst_sha256':sh,'s3_uri':uri}
            append_tsv(acceptedp,FIELDS,rec); accepted.append(rec); aset.add(acc); got+=bp; last[acc]='ACCEPTED'
            log(rec,'ACCEPTED','accepted under unchanged R1 rules',bytes_n,bp,nr)
            print(f'[ACCEPT] {STR} rank={rank} {acc} bp={bp} cumulative={got}/{TARGET[STR]}',flush=True)
    if got<TARGET[STR]: raise SystemExit(f'FATAL: exhausted eligible deterministic candidates for {STR}: {got} < {TARGET[STR]}')
    rows=sorted(read_tsv(acceptedp),key=lambda r:int(r['deterministic_rank']))
    outf=S/'frozen/shards'/f'{STR}_SELECTED.tsv'; outf.parent.mkdir(parents=True,exist_ok=True)
    with open(outf,'w',newline='') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=FIELDS); w.writeheader(); w.writerows(rows)
    meta={'stage':'Stage05F3V2-R2-distributed-acquisition','stratum':STR,'gate_pass':True,'gate':f'PASS_F3V2_R2_STRATUM_{STR}_QUOTA_COMPLETE','target_bp':TARGET[STR],'accepted_bp':got,'accepted_n':len(rows),'candidate_n':total,'source_r1_stage':str(R),'candidate_table_sha256':sha256(cp),'selected_sha256':sha256(outf),'host':socket.gethostname(),'science_change':False,'claim_boundary':'Operational parallelization across independent frozen strata only; within-stratum deterministic order and all R1 acquisition rules are unchanged.'}
    (S/'frozen/shards'/f'{STR}_FREEZE.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps(meta,indent=2),flush=True); print(f'PASS_F3V2_R2_STRATUM_{STR}_COMPLETE',flush=True)

if __name__=='__main__': main()
