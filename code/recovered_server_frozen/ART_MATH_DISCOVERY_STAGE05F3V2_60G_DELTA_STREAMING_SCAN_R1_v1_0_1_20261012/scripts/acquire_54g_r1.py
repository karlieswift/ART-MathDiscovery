#!/usr/bin/env python3
import argparse,csv,hashlib,json,os,re,shutil,subprocess,sys,time
from pathlib import Path

TARGET={'A':13_500_000_000,'B':18_000_000_000,'C':13_500_000_000,'D':9_000_000_000}
FIELDS=['stratum','accession','nucleotide_bp','contig_n','compressed_bytes','zst_sha256','s3_uri']
LEDGER_FIELDS=['event_id','stratum','deterministic_rank','accession','status','compressed_bytes','nucleotide_bp','contig_n','detail']
PERMANENT={'LOGAN_OBJECT_NOT_FOUND','COMPRESSED_SIZE_OUT_OF_RANGE','EMPTY_OR_INVALID_FASTA','REALIZED_BP_OUT_OF_RANGE'}
RETRYABLE={'HEAD_FAILED_TRANSIENT','DOWNLOAD_FAILED','DECOMPRESS_FAILED'}

def eprint(*x): print(*x,file=sys.stderr,flush=True)
def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def read_tsv(p):
    if not Path(p).is_file() or Path(p).stat().st_size==0: return []
    with open(p,newline='') as f: return list(csv.DictReader(f,delimiter='\t'))

def init_tsv(p,fields):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    if not p.exists() or p.stat().st_size==0:
        with open(p,'w',newline='') as f: csv.DictWriter(f,delimiter='\t',fieldnames=fields).writeheader()

def append_tsv(p,fields,row):
    with open(p,'a',newline='') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=fields,extrasaction='ignore'); w.writerow(row); f.flush(); os.fsync(f.fileno())

def run(cmd,**kw): return subprocess.run(cmd,text=True,capture_output=True,**kw)

def count_fasta(zst,count_script):
    p1=subprocess.Popen(['zstd','-dc',str(zst)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    p2=subprocess.Popen([sys.executable,str(count_script)],stdin=p1.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    p1.stdout.close()
    out2,err2=p2.communicate(); err1=p1.stderr.read().decode(errors='replace'); rc1=p1.wait()
    if rc1!=0 or p2.returncode!=0: raise RuntimeError((err1+' '+err2).strip()[:500])
    j=json.loads(out2); return int(j['bp']),int(j['records'])

def link_or_keep(src,dst):
    src=Path(src).resolve(); dst=Path(dst)
    if dst.exists() or dst.is_symlink():
        try:
            if os.path.samefile(src,dst): return 'existing-samefile'
        except Exception: pass
        if dst.is_symlink() or dst.is_file(): dst.unlink()
        else: raise RuntimeError(f'cannot replace non-file {dst}')
    try:
        os.link(src,dst)
        return 'hardlink'
    except OSError:
        dst.symlink_to(src)
        return 'symlink'

def load_last_status(ledger):
    last={}
    for r in read_tsv(ledger): last[r['accession']]=r['status']
    return last

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--stage',required=True)
    ap.add_argument('--resume-from',default='')
    ap.add_argument('--max-accession-bp',type=int,default=int(os.environ.get('MAX_ACCESSION_BP','1500000000')))
    ap.add_argument('--min-zst-bytes',type=int,default=int(os.environ.get('MIN_ZST_BYTES','10000000')))
    ap.add_argument('--max-zst-bytes',type=int,default=int(os.environ.get('MAX_ZST_BYTES','2000000000')))
    ap.add_argument('--transport-fail-limit',type=int,default=int(os.environ.get('F3V2_TRANSPORT_FAIL_LIMIT','5')))
    ap.add_argument('--download-retries',type=int,default=int(os.environ.get('F3V2_DOWNLOAD_RETRIES','3')))
    a=ap.parse_args()
    S=Path(a.stage).resolve(); resume=Path(a.resume_from).resolve() if a.resume_from else None
    oldp=S/'inputs/OLD_F3_ACCESSIONS.txt'; old={x.strip() for x in open(oldp) if x.strip() and not x.startswith('#')}
    datadir=S/'data/contigs'; datadir.mkdir(parents=True,exist_ok=True)
    work=S/'work/download'; work.mkdir(parents=True,exist_ok=True)
    acceptedp=work/'accepted.partial.tsv'; ledgerp=work/'F3V2_ACQUISITION_LEDGER.tsv'
    init_tsv(acceptedp,FIELDS); init_tsv(ledgerp,LEDGER_FIELDS)
    event_id=sum(1 for _ in open(ledgerp))-1
    def log(row,status,detail='',compressed_bytes='',bp='',nr=''):
        nonlocal event_id
        event_id+=1
        append_tsv(ledgerp,LEDGER_FIELDS,{'event_id':event_id,'stratum':row.get('stratum',''),'deterministic_rank':row.get('deterministic_rank',''),'accession':row.get('accession',''),'status':status,'compressed_bytes':compressed_bytes,'nucleotide_bp':bp,'contig_n':nr,'detail':detail[:500]})

    # Import already accepted v1.0.0 rows exactly once, only when the actual contig exists and SHA validates.
    accepted=read_tsv(acceptedp); aset={r['accession'] for r in accepted}
    if resume and resume.is_dir():
        candidates=[resume/'frozen/F3V2_DELTA_SELECTED_ACCESSIONS.tsv',resume/'work/download/selected.tmp.tsv']
        srcmanifest=next((p for p in candidates if p.is_file() and p.stat().st_size>0),None)
        if srcmanifest:
            imported=0
            for r in read_tsv(srcmanifest):
                acc=r.get('accession','')
                if not acc or acc in aset: continue
                if acc in old: raise SystemExit(f'FATAL: v1.0.0 partial contains old-universe accession {acc}')
                src=resume/'data/contigs'/f'{acc}.contigs.fa.zst'
                if not src.is_file() or src.stat().st_size==0:
                    eprint(f'[R1 IMPORT SKIP] missing old contig {src}'); continue
                expected=r.get('zst_sha256','')
                actual=sha256(src)
                if expected and expected!=actual:
                    eprint(f'[R1 IMPORT SKIP] SHA mismatch {acc}'); continue
                dst=datadir/f'{acc}.contigs.fa.zst'; link_mode=link_or_keep(src,dst)
                row={k:r.get(k,'') for k in FIELDS}
                append_tsv(acceptedp,FIELDS,row); accepted.append(row); aset.add(acc); imported+=1
                log({'stratum':row['stratum'],'deterministic_rank':'','accession':acc},'ACCEPTED_IMPORTED_V100',f'SHA-verified existing v1.0.0 partial; {link_mode} reuse',row['compressed_bytes'],row['nucleotide_bp'],row['contig_n'])
            print(f'[R1 IMPORT] imported_v100_accepted={imported}',flush=True)

    # Deduplicate accepted state defensively and calculate current quotas.
    uniq=[]; seen=set()
    for r in accepted:
        if r['accession'] in seen: continue
        seen.add(r['accession']); uniq.append(r)
    if len(uniq)!=len(accepted):
        tmp=acceptedp.with_suffix('.tmp')
        with open(tmp,'w',newline='') as f:
            w=csv.DictWriter(f,delimiter='\t',fieldnames=FIELDS); w.writeheader(); w.writerows(uniq)
        os.replace(tmp,acceptedp)
    accepted=uniq; aset={r['accession'] for r in accepted}
    got={s:sum(int(r['nucleotide_bp']) for r in accepted if r['stratum']==s) for s in 'ABCD'}
    print('[R1 RESUME] '+ ' '.join(f'{s}={got[s]}/{TARGET[s]}' for s in 'ABCD'),flush=True)
    last=load_last_status(ledgerp)
    count_script=S/'scripts/count_fasta_stream.py'
    transport_fail=0

    for STR in 'ABCD':
        if got[STR]>=TARGET[STR]:
            print(f'[STRATUM {STR}] already complete {got[STR]}/{TARGET[STR]}',flush=True); continue
        cp=S/'inputs'/f'F3V2_DELTA_CANDIDATES_{STR}.tsv'
        if not cp.is_file(): raise SystemExit(f'FATAL missing candidate table {cp}; run 01 first')
        total=max(sum(1 for _ in open(cp))-1,0)
        print(f'[STRATUM {STR}] target_bp={TARGET[STR]} deterministic_candidates={total} resume_bp={got[STR]}',flush=True)
        local_seen=set()
        with open(cp,newline='') as f:
            for row in csv.DictReader(f,delimiter='\t'):
                if got[STR]>=TARGET[STR]: break
                acc=row['accession']; rank=row['deterministic_rank']
                if acc in local_seen:
                    log(row,'DUPLICATE_ACCESSION','duplicate accession inside deterministic candidate table'); continue
                local_seen.add(acc)
                if acc in old:
                    log(row,'ALREADY_IN_OLD_UNIVERSE','fatal invariant violation: old accession leaked into delta')
                    raise SystemExit(f'FATAL old accession leaked into delta: {acc}')
                if acc in aset: continue
                if last.get(acc) in PERMANENT: continue
                log(row,'SELECTED','candidate reached by deterministic acquisition stream')
                key=f'c/{acc}/{acc}.contigs.fa.zst'; uri=f's3://logan-pub/{key}'; z=datadir/f'{acc}.contigs.fa.zst'
                bytes_n=0
                if not z.is_file() or z.stat().st_size==0:
                    q=run(['aws','s3api','head-object','--bucket','logan-pub','--key',key,'--no-sign-request','--query','ContentLength','--output','text'])
                    if q.returncode!=0:
                        msg=(q.stderr or q.stdout or '').strip().replace('\n',' ')
                        if re.search(r'(^|\D)404(\D|$)|NoSuchKey|Not Found',msg,re.I):
                            log(row,'LOGAN_OBJECT_NOT_FOUND',msg); last[acc]='LOGAN_OBJECT_NOT_FOUND'; transport_fail=0; continue
                        log(row,'HEAD_FAILED_TRANSIENT',msg); last[acc]='HEAD_FAILED_TRANSIENT'; transport_fail+=1
                        eprint(f'[TRANSIENT HEAD] {acc} consecutive={transport_fail}/{a.transport_fail_limit}: {msg[:180]}')
                        if transport_fail>=a.transport_fail_limit: raise SystemExit('FATAL: repeated S3 transport/head failures; acquisition stopped fail-closed. Rerun the same R1 command after network recovery.')
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
                        msg=(d.stderr or d.stdout or '').strip().replace('\n',' '); part.unlink(missing_ok=True)
                        time.sleep(min(attempt*2,5))
                    if not ok:
                        log(row,'DOWNLOAD_FAILED',msg,bytes_n); last[acc]='DOWNLOAD_FAILED'; transport_fail+=1
                        if transport_fail>=a.transport_fail_limit: raise SystemExit('FATAL: repeated S3 download failures; acquisition stopped fail-closed. Rerun the same R1 command after network recovery.')
                        continue
                    transport_fail=0
                else:
                    bytes_n=z.stat().st_size
                try:
                    bp,nr=count_fasta(z,count_script)
                except Exception as ex:
                    if z.is_symlink(): z.unlink(missing_ok=True)
                    else: z.unlink(missing_ok=True)
                    log(row,'DECOMPRESS_FAILED',str(ex),bytes_n); last[acc]='DECOMPRESS_FAILED'; continue
                if bp<=0 or nr<=0:
                    z.unlink(missing_ok=True); log(row,'EMPTY_OR_INVALID_FASTA','zero bp or zero records',bytes_n,bp,nr); last[acc]='EMPTY_OR_INVALID_FASTA'; continue
                if bp>a.max_accession_bp:
                    z.unlink(missing_ok=True); log(row,'REALIZED_BP_OUT_OF_RANGE',f'bp>{a.max_accession_bp}',bytes_n,bp,nr); last[acc]='REALIZED_BP_OUT_OF_RANGE'; continue
                sh=sha256(z)
                rec={'stratum':STR,'accession':acc,'nucleotide_bp':bp,'contig_n':nr,'compressed_bytes':bytes_n,'zst_sha256':sh,'s3_uri':uri}
                append_tsv(acceptedp,FIELDS,rec); accepted.append(rec); aset.add(acc); got[STR]+=bp
                log(row,'ACCEPTED','accepted under frozen technical-size and realized-bp rules',bytes_n,bp,nr); last[acc]='ACCEPTED'
                print(f'[ACCEPT] {STR} rank={rank} {acc} bp={bp} cumulative={got[STR]}/{TARGET[STR]}',flush=True)
        if got[STR]<TARGET[STR]:
            raise SystemExit(f'FATAL: exhausted all eligible deterministic candidates for stratum {STR}: {got[STR]} < {TARGET[STR]}. This is a source-universe limitation, not a K ceiling.')

    # Finalize only after every frozen quota is met.
    accepted=read_tsv(acceptedp)
    overlap={r['accession'] for r in accepted}&old
    if overlap: raise SystemExit(f'FATAL old-universe overlap at freeze: {sorted(overlap)[:10]}')
    if len({r['accession'] for r in accepted})!=len(accepted): raise SystemExit('FATAL duplicate accepted accessions')
    by={s:sum(int(r['nucleotide_bp']) for r in accepted if r['stratum']==s) for s in 'ABCD'}
    if not all(by[s]>=TARGET[s] for s in 'ABCD'): raise SystemExit(f'FATAL quota invariant failed: {by}')
    final=S/'frozen/F3V2_DELTA_SELECTED_ACCESSIONS.tsv'
    with open(final,'w',newline='') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=FIELDS); w.writeheader(); w.writerows(accepted)
    bp=sum(by.values())
    gate=(50_000_000_000<=bp<=60_000_000_000)
    ledger_final=S/'frozen/F3V2_ACQUISITION_LEDGER.tsv'
    shutil.copy2(ledgerp,ledger_final)
    ledger_rows=read_tsv(ledgerp)
    imported_rows=[r for r in ledger_rows if r.get('status')=='ACCEPTED_IMPORTED_V100']
    imported_bp=sum(int(r.get('nucleotide_bp') or 0) for r in imported_rows)
    g={'stage':'Stage05F3V2-R1-source-freeze','gate':'PASS_F3V2_R1_54G_DELTA_SOURCE_BLIND_FREEZE' if gate else 'FAIL_F3V2_R1_TOTAL_BP_GATE','gate_pass':gate,
       'delta_accession_n':len(accepted),'delta_nucleotide_bp':bp,'stratum_bp':by,'stratum_targets_bp':TARGET,'old_accession_overlap_n':0,
       'source_selection_used_rt_or_candidate_information':False,'selection_inputs':'SRA accession/assay/organism/taxonomy plus Logan technical file size and realized nucleotide bp only',
       'contigs_primary':True,'unitigs_used':False,'selected_manifest_sha256':sha256(final),'acquisition_ledger_sha256':sha256(ledger_final),
       'resume_from_v100':str(resume) if resume else None,'v100_imported_acceptance_n':len(imported_rows),'v100_imported_nucleotide_bp':imported_bp,
       'r1_change':'acquisition-only: removed fixed K=5000 ceiling; deterministic v1.0.0 hash order retained; SHA-verified accepted v1.0.0 partial is resumed without redownload; transport errors fail closed',
       'claim_boundary':'This freezes a post-discovery external delta universe before any F3V2 HMM result is inspected; it does not modify the original 6.12-Gb universe or any scan/novelty threshold.'}
    (S/'frozen/F3V2_DELTA_SOURCE_FREEZE.json').write_text(json.dumps(g,indent=2)+'\n')
    print(json.dumps(g,indent=2),flush=True)
    if not gate: raise SystemExit(2)
    try:
        final.chmod(0o444); ledger_final.chmod(0o444); (S/'frozen/F3V2_DELTA_SOURCE_FREEZE.json').chmod(0o444)
    except Exception: pass
    print('PASS_F3V2_R1_54G_DELTA_DOWNLOAD_AND_SOURCE_FREEZE_COMPLETE',flush=True)
if __name__=='__main__': main()
