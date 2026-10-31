#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, os, re, shutil, subprocess, sys, time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

EXPECTED_FORMAL=88487
EXPECTED_UNIQUE=85503
EXPECTED_FAMILIES=26396
EXPECTED_DELTA_ONLY=46323
EXPECTED_MIXED=42164
KSET=[8,10,12,14,16]
STARTS={"ATG","GTG","TTG"}; STOPS={"TAA","TAG","TGA"}
CODON={
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
COMP=str.maketrans('ACGTacgt','TGCAtgca')

def rc(s): return s.translate(COMP)[::-1].upper()
def translate(nt):
    nt=nt.upper(); return ''.join(CODON.get(nt[i:i+3],'X') for i in range(0,len(nt)-2,3))
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def write_json(p,x): Path(p).write_text(json.dumps(x,indent=2,sort_keys=True)+"\n")
def tsv_read(p):
    with open(p,newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
def tsv_write(p,rows,fields):
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
def parse_cid(cid):
    z=cid.split('|')
    if len(z)<7: raise ValueError(f'candidate_id format: {cid}')
    acc,contig,strand,frame,coord,aalen,ordv=z[:7]
    m=re.fullmatch(r'(\d+)-(\d+)',coord); ma=re.fullmatch(r'aa(\d+)',aalen)
    if not m or not ma or strand not in ['+','-']: raise ValueError(cid)
    a,b=map(int,m.groups()); aa=int(ma.group(1))
    return dict(candidate_id=cid,accession=acc,contig_id=contig,strand=strand,frame=frame,start0=a,end0=b,aa_length=aa,ordinal=ordv)

def iter_fasta(path):
    p=Path(path); proc=None
    if str(p).endswith('.zst'):
        proc=subprocess.Popen(['zstdcat',str(p)],stdout=subprocess.PIPE,text=True,bufsize=1); fh=proc.stdout
    elif str(p).endswith('.gz'): fh=gzip.open(p,'rt')
    else: fh=open(p)
    name=None; buf=[]
    try:
        for line in fh:
            if line.startswith('>'):
                if name is not None: yield name,''.join(buf).upper()
                name=line[1:].strip().split()[0]; buf=[]
            else: buf.append(line.strip())
        if name is not None: yield name,''.join(buf).upper()
    finally:
        try: fh.close()
        except Exception: pass
        if proc is not None:
            if proc.poll() is None:
                try: proc.terminate()
                except Exception: pass
            try: proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill(); proc.wait()

def fasta_count(path):
    n=0
    opener=gzip.open if str(path).endswith('.gz') else open
    if str(path).endswith('.zst'):
        p=subprocess.Popen(['zstdcat',str(path)],stdout=subprocess.PIPE,text=True)
        try:
            for line in p.stdout:
                if line.startswith('>'): n+=1
        finally:
            p.stdout.close(); p.wait()
        return n
    with opener(path,'rt' if str(path).endswith('.gz') else 'r') as f:
        for line in f:
            if line.startswith('>'): n+=1
    return n

def find_selected_ledger(P):
    hits=list(P.glob('**/frozen/F3V2_DELTA_SELECTED_ACCESSIONS.tsv'))
    scored=[]
    for h in hits:
        try: n=sum(1 for _ in open(h))-1
        except: continue
        s=(1000 if 'DISTRIBUTED_ACQUISITION_R2' in str(h) else 0)+(100 if n==254 else 0)
        scored.append((s,str(h),h))
    if not scored: raise RuntimeError('F3V2_DELTA_SELECTED_ACCESSIONS.tsv not found')
    return sorted(scored,reverse=True)[0][2]

def ledger_map(ledger):
    out={}
    for r in tsv_read(ledger):
        acc=r.get('accession') or r.get('run_accession') or r.get('sra_accession')
        if acc: out[acc]=r
    return out

def find_source_file(source_root,acc,row,ledger_path):
    for k,v in row.items():
        if v and ('path' in k.lower() or 'file' in k.lower()):
            q=Path(v)
            for qq in ([q] if q.is_absolute() else [source_root/q,ledger_path.parent/q]):
                if qq.exists() and qq.is_file(): return qq.resolve()
    pats=[f'{acc}.contigs.fa.zst',f'{acc}.fa.zst',f'*{acc}*contigs*.zst']
    for pat in pats:
        hs=list(source_root.glob('**/'+pat))
        if hs: return hs[0].resolve()
    raise FileNotFoundError(f'no source FASTA for {acc}')

def discover_formal_faa(P,override=None):
    if override:
        p=Path(override).resolve(); n=fasta_count(p)
        if n!=EXPECTED_FORMAL: raise RuntimeError(f'override formal FASTA count {n} != {EXPECTED_FORMAL}')
        return p
    candidates=[]
    for root in P.iterdir():
        if not root.is_dir() or 'STAGE05F3V2' not in root.name: continue
        for fp in root.glob('**/*'):
            if not fp.is_file() or fp.suffix.lower() not in ('.faa','.fa','.fasta'): continue
            low=fp.name.lower()
            if not any(x in low for x in ('formal','rvt','candidate')): continue
            try:
                n=fasta_count(fp)
            except Exception:
                continue
            if n==EXPECTED_FORMAL:
                score=(100 if 'formal' in low else 0)+(50 if 'stage03' in str(fp).lower() else 0)-(20 if 'b4' in str(fp).lower() else 0)
                candidates.append((score,str(fp),fp))
    if not candidates: raise RuntimeError('HOLD: could not auto-discover an 88,487-record formal FASTA; pass --formal-faa PATH')
    return sorted(candidates,reverse=True)[0][2].resolve()

def read_fasta_dict(p): return {n:s for n,s in iter_fasta(p)}

def discover_delta_only_table(P,formal_ids):
    best=None
    for root in P.iterdir():
        if not root.is_dir() or 'STAGE04B2' not in root.name: continue
        for fp in root.glob('**/*.tsv'):
            try:
                with open(fp) as f: hdr=f.readline().rstrip('\n').split('\t')
            except: continue
            if 'candidate_id' not in hdr or 'family_rep_id' not in hdr: continue
            try: rows=tsv_read(fp)
            except: continue
            rr=[r for r in rows if r.get('candidate_id') in formal_ids]
            n=len({r['candidate_id'] for r in rr})
            if n and (best is None or n>best[0]): best=(n,fp,rr)
    if best is None or best[0]!=EXPECTED_DELTA_ONLY:
        raise RuntimeError(f'HOLD: delta-only characterization unresolved; best={None if best is None else (best[0],str(best[1]))}, expected={EXPECTED_DELTA_ONLY}')
    mp={r['candidate_id']:r['family_rep_id'] for r in best[2]}
    return mp,best[1].resolve()

def repeat_arrays(seq,rt0,rt1):
    found={}; n=len(seq)
    for k in KSET:
        pos=defaultdict(list)
        for i in range(0,n-k+1): pos[seq[i:i+k]].append(i)
        for seed,pp in pos.items():
            if len(pp)<3: continue
            pp=sorted(pp); run=[pp[0]]; runs=[]
            for x in pp[1:]:
                if 60 <= x-run[-1] <= 360: run.append(x)
                else:
                    if len(run)>=3: runs.append(run)
                    run=[x]
            if len(run)>=3: runs.append(run)
            for starts in runs:
                L=0
                while min(starts)-L-1>=0 and k+L<80:
                    ch=seq[starts[0]-L-1]
                    if all(seq[s-L-1]==ch for s in starts): L+=1
                    else: break
                R=0
                while max(starts)+k+R<n and k+L+R<80:
                    ch=seq[starts[0]+k+R]
                    if all(seq[s+k+R]==ch for s in starts): R+=1
                    else: break
                ss=tuple(s-L for s in starts); motif=seq[ss[0]:ss[0]+k+L+R]
                prev=found.get(ss)
                if prev is None or len(motif)>len(prev['motif']):
                    gaps=[ss[i+1]-ss[i] for i in range(len(ss)-1)]
                    spacers=[g-len(motif) for g in gaps]
                    a0=ss[0]; a1=ss[-1]+len(motif)
                    orient='UPSTREAM' if a1<=rt0 else ('DOWNSTREAM' if a0>=rt1 else 'OVERLAP_RT')
                    strict=(orient=='UPSTREAM' and 3<=len(ss)<=21 and all(120<=z<=220 for z in spacers))
                    found[ss]=dict(motif=motif,motif_len=len(motif),copy_n=len(ss),starts=list(ss),start_gaps=gaps,spacers=spacers,array_start0=a0,array_end0=a1,orientation=orient,strict_published_art_geometry=strict)
    return sorted(found.values(),key=lambda x:(x['array_start0'],-x['motif_len'],-x['copy_n']))

def find_orfs_downstream(seq,rt1,maxbp=10000):
    sub=seq[rt1:min(len(seq),rt1+maxbp)]; ans=[]
    for frame in range(3):
        for i in range(frame,len(sub)-2,3):
            if sub[i:i+3] not in STARTS: continue
            for j in range(i+3,len(sub)-2,3):
                if sub[j:j+3] in STOPS:
                    aa=translate(sub[i:j])
                    if len(aa)>=30: ans.append((i,j+3,'+',aa))
                    break
    rsub=rc(sub); L=len(sub)
    for frame in range(3):
        for i in range(frame,len(rsub)-2,3):
            if rsub[i:i+3] not in STARTS: continue
            for j in range(i+3,len(rsub)-2,3):
                if rsub[j:j+3] in STOPS:
                    aa=translate(rsub[i:j])
                    if len(aa)>=30: ans.append((L-(j+3),L-i,'-',aa))
                    break
    ans.sort(key=lambda z:(z[0],z[1]-z[0],z[2])); return ans

def preflight(args):
    H=Path(args.home).resolve(); P=Path(args.project_root).resolve(); W=H/'work'; W.mkdir(exist_ok=True)
    ledger=find_selected_ledger(P); lm=ledger_map(ledger)
    if len(lm)!=254: raise SystemExit(f'HOLD selected ledger accession count {len(lm)} != 254')
    SR=ledger.parents[1]
    formal=discover_formal_faa(P,args.formal_faa)
    aa=read_fasta_dict(formal)
    if len(aa)!=EXPECTED_FORMAL: raise SystemExit(f'HOLD unique formal headers {len(aa)} != {EXPECTED_FORMAL}')
    formal_ids=set(aa)
    delta_only_map,delta_only_src=discover_delta_only_table(P,formal_ids)
    cls=Counter('DELTA_ONLY' if x in delta_only_map else 'MIXED' for x in formal_ids)
    if cls['DELTA_ONLY']!=EXPECTED_DELTA_ONLY or cls['MIXED']!=EXPECTED_MIXED:
        raise SystemExit(f'HOLD class counts {dict(cls)} expected DELTA_ONLY={EXPECTED_DELTA_ONLY} MIXED={EXPECTED_MIXED}')
    # exact AA dedupe and per-accession manifests
    seqrep={}; cand_to_unique={}
    for cid,seq in aa.items():
        rep=seqrep.setdefault(seq,cid); cand_to_unique[cid]=rep
    if len(seqrep)!=EXPECTED_UNIQUE: raise SystemExit(f'HOLD exact unique AA {len(seqrep)} != {EXPECTED_UNIQUE}')
    inp=W/'accession_inputs'; inp.mkdir(exist_ok=True)
    rows_by_acc=defaultdict(list)
    for cid,seq in aa.items():
        pc=parse_cid(cid); acc=pc['accession']
        if acc not in lm: raise SystemExit(f'HOLD formal accession {acc} absent selected ledger')
        rows_by_acc[acc].append({**pc,'aa_sequence':seq,'candidate_class':('DELTA_ONLY' if cid in delta_only_map else 'MIXED'),'combined_delta_only_family_id':delta_only_map.get(cid,''),'unique_aa_rep_id':cand_to_unique[cid]})
    if sum(map(len,rows_by_acc.values()))!=EXPECTED_FORMAL: raise SystemExit('HOLD row accounting mismatch')
    # write unique AA FASTA for exact frozen family clustering
    uf=W/'FORMAL_UNIQUE_AA_85503.faa'
    with open(uf,'w') as f:
        for seq,cid in seqrep.items(): f.write(f'>{cid}\n{seq}\n')
    acc_manifest=[]
    fields=['candidate_id','accession','contig_id','strand','frame','start0','end0','aa_length','ordinal','aa_sequence','candidate_class','combined_delta_only_family_id','unique_aa_rep_id']
    for acc in sorted(rows_by_acc):
        sf=find_source_file(SR,acc,lm[acc],ledger)
        ap=inp/f'{acc}.tsv'; tsv_write(ap,rows_by_acc[acc],fields)
        acc_manifest.append(dict(accession=acc,stratum=lm[acc].get('stratum',''),candidate_n=len(rows_by_acc[acc]),source_file=str(sf),ledger_zst_sha256=lm[acc].get('zst_sha256',''),input_manifest=str(ap)))
    tsv_write(W/'ACCESSION_SCAN_MANIFEST.tsv',acc_manifest,['accession','stratum','candidate_n','source_file','ledger_zst_sha256','input_manifest'])
    # Selected accessions with zero frozen formal RT-like candidates are valid zero-results.
    # They have nothing to scan and therefore intentionally have no accession_results directory.
    nonzero_acc={r['accession'] for r in acc_manifest}
    zero_rows=[]
    for acc in sorted(set(lm)-nonzero_acc):
        zero_rows.append(dict(accession=acc,stratum=lm[acc].get('stratum',''),formal_candidate_n=0,reason='ZERO_FROZEN_FORMAL_RT_CANDIDATES'))
    if len(acc_manifest)+len(zero_rows)!=len(lm):
        raise SystemExit('HOLD selected/nonzero/zero accession accounting mismatch')
    tsv_write(W/'ZERO_FORMAL_ACCESSIONS.tsv',zero_rows,['accession','stratum','formal_candidate_n','reason'])
    # candidate class map
    cmap=[]
    for cid in sorted(formal_ids):
        cmap.append(dict(candidate_id=cid,candidate_class=('DELTA_ONLY' if cid in delta_only_map else 'MIXED'),combined_delta_only_family_id=delta_only_map.get(cid,''),unique_aa_rep_id=cand_to_unique[cid]))
    tsv_write(W/'FORMAL_CANDIDATE_CLASS_MAP.tsv',cmap,['candidate_id','candidate_class','combined_delta_only_family_id','unique_aa_rep_id'])
    gate={'gate':'PASS_STAGE04B5_PREFLIGHT_88487_FORMAL_FROZEN','gate_pass':True,'formal_faa':str(formal),'formal_faa_sha256':sha(formal),'formal_instance_n':len(aa),'exact_aa_unique_n':len(seqrep),'delta_only_instance_n':cls['DELTA_ONLY'],'mixed_instance_n':cls['MIXED'],'selected_accession_n':len(lm),'formal_bearing_accession_n':len(acc_manifest),'zero_formal_accession_n':len(zero_rows),'selected_ledger':str(ledger),'delta_only_source':str(delta_only_src),'source_root':str(SR),'source_reselection_performed':False,'retuning_performed':False}
    write_json(W/'PREFLIGHT_GATE.json',gate); print(json.dumps(gate,indent=2))

def scan_accession_job(job):
    acc=job['accession']; inp=Path(job['input_manifest']); sf=Path(job['source_file']); out=Path(job['outdir'])
    out.mkdir(parents=True,exist_ok=True); done=out/'DONE.json'
    if done.exists():
        try:
            d=json.load(open(done))
            if d.get('gate_pass') and d.get('candidate_n')==int(job['candidate_n']): return d
        except Exception: pass
    rows=tsv_read(inp); bycontig=defaultdict(list)
    for r in rows: bycontig[r['contig_id']].append(r)
    want=set(bycontig); found=set(); contexts=[]; arrays=[]; partners=[]; pfa=[]
    processed=0
    for name,genome in iter_fasta(sf):
        if name not in want: continue
        found.add(name)
        for r in bycontig[name]:
            cid=r['candidate_id']; a=int(r['start0']); b=int(r['end0']); aa_len=int(r['aa_length']); strand=r['strand']; aa=r['aa_sequence']
            if not (0<=a<b<=len(genome)): raise RuntimeError(f'{acc} coords out of range {cid} contig_len={len(genome)}')
            nt=genome[a:b]; ori=nt if strand=='+' else rc(nt)
            if len(ori)!=aa_len*3: raise RuntimeError(f'{acc} nt/aa length mismatch {cid}')
            tr=translate(ori).rstrip('*')
            if tr!=aa.rstrip('*'): raise RuntimeError(f'{acc} exact peptide replay mismatch {cid}')
            left=max(0,a-10000); right=min(len(genome),b+10000); rawctx=genome[left:right]
            if strand=='+': ctx=rawctx; rt0=a-left; rt1=b-left
            else: ctx=rc(rawctx); rt0=right-b; rt1=right-a
            up=rt0; down=len(ctx)-rt1
            contexts.append(dict(candidate_id=cid,accession=acc,stratum=job.get('stratum',''),candidate_class=r['candidate_class'],combined_delta_only_family_id=r['combined_delta_only_family_id'],unique_aa_rep_id=r['unique_aa_rep_id'],contig_id=name,contig_length_bp=len(genome),strand=strand,aa_length=aa_len,current_context_bp=len(ctx),upstream_bp=up,downstream_bp=down,strict_array_evaluable=str(up>=264),partner_evaluable=str(down>=93),both_flanks_500bp=str(up>=500 and down>=500),both_flanks_2000bp=str(up>=2000 and down>=2000),both_flanks_5000bp=str(up>=5000 and down>=5000),both_flanks_10000bp=str(up>=10000 and down>=10000)))
            for ix,z in enumerate(repeat_arrays(ctx,rt0,rt1),1):
                arrays.append(dict(candidate_id=cid,accession=acc,stratum=job.get('stratum',''),candidate_class=r['candidate_class'],combined_delta_only_family_id=r['combined_delta_only_family_id'],array_index=ix,motif=z['motif'],motif_len=z['motif_len'],copy_n=z['copy_n'],repeat_starts=','.join(map(str,z['starts'])),start_gaps=','.join(map(str,z['start_gaps'])),inferred_spacers=','.join(map(str,z['spacers'])),orientation=z['orientation'],strict_published_art_geometry=str(z['strict_published_art_geometry'])))
            orfs=find_orfs_downstream(ctx,rt1)
            if orfs:
                oa,ob,ost,opa=orfs[0]; pid=f'{acc}__P{len(partners)+1:07d}'
                partners.append(dict(partner_id=pid,candidate_id=cid,accession=acc,stratum=job.get('stratum',''),candidate_class=r['candidate_class'],combined_delta_only_family_id=r['combined_delta_only_family_id'],orf_strand=ost,orf_start0_downstream=oa,orf_end0_downstream=ob,aa_length=len(opa),rt_to_partner_gap_bp=oa))
                pfa.append((pid,opa))
            processed+=1
        if found==want: break
    missing=sorted(want-found)
    if missing: raise RuntimeError(f'{acc} missing {len(missing)} candidate contigs first={missing[:3]}')
    if processed!=len(rows): raise RuntimeError(f'{acc} processed {processed} != expected {len(rows)}')
    cf=['candidate_id','accession','stratum','candidate_class','combined_delta_only_family_id','unique_aa_rep_id','contig_id','contig_length_bp','strand','aa_length','current_context_bp','upstream_bp','downstream_bp','strict_array_evaluable','partner_evaluable','both_flanks_500bp','both_flanks_2000bp','both_flanks_5000bp','both_flanks_10000bp']
    af=['candidate_id','accession','stratum','candidate_class','combined_delta_only_family_id','array_index','motif','motif_len','copy_n','repeat_starts','start_gaps','inferred_spacers','orientation','strict_published_art_geometry']
    pf=['partner_id','candidate_id','accession','stratum','candidate_class','combined_delta_only_family_id','orf_strand','orf_start0_downstream','orf_end0_downstream','aa_length','rt_to_partner_gap_bp']
    tsv_write(out/'contexts.tsv',contexts,cf); tsv_write(out/'arrays.tsv',arrays,af); tsv_write(out/'partners.tsv',partners,pf)
    with open(out/'partners.faa','w') as f:
        for n,s in pfa: f.write(f'>{n}\n{s}\n')
    d={'gate':'PASS_STAGE04B5_ACCESSION_SCAN','gate_pass':True,'accession':acc,'candidate_n':processed,'target_contig_n':len(want),'repeat_array_n':len(arrays),'partner_orf_n':len(partners),'source_file':str(sf),'source_file_size':sf.stat().st_size}
    write_json(done,d); return d

def run_scan(args):
    H=Path(args.home).resolve(); W=H/'work'; M=tsv_read(W/'ACCESSION_SCAN_MANIFEST.tsv')
    ad=W/'accession_results'; ad.mkdir(exist_ok=True)
    jobs=[]
    for r in M: jobs.append({**r,'outdir':str(ad/r['accession'])})
    print(f'[STAGE04B5] accessions={len(jobs)} formal={sum(int(x["candidate_n"]) for x in jobs)} workers={args.workers}',flush=True)
    done=[]; t0=time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        fut={ex.submit(scan_accession_job,j):j['accession'] for j in jobs}
        for i,f in enumerate(as_completed(fut),1):
            acc=fut[f]
            try: d=f.result()
            except Exception as e: raise SystemExit(f'HOLD accession {acc} scan failure: {e}')
            done.append(d); print(f'[PASS ACCESSION {i}/{len(jobs)}] {acc} candidates={d["candidate_n"]} arrays={d["repeat_array_n"]} partners={d["partner_orf_n"]} elapsed={time.time()-t0:.1f}s',flush=True)
    if sum(x['candidate_n'] for x in done)!=EXPECTED_FORMAL: raise SystemExit('HOLD all-accession candidate count mismatch')
    zero=tsv_read(W/'ZERO_FORMAL_ACCESSIONS.tsv')
    selected_n=len(done)+len(zero)
    if selected_n!=254: raise SystemExit(f'HOLD selected accession accounting {len(done)}+{len(zero)} != 254')
    write_json(W/'SCAN_GATE.json',{'gate':'PASS_STAGE04B5_ALL_FORMAL_ACCESSION_SCAN_COMPLETE','gate_pass':True,'selected_accession_n':selected_n,'formal_bearing_accession_n':len(done),'zero_formal_accession_n':len(zero),'formal_instance_n':sum(x['candidate_n'] for x in done),'elapsed_seconds':time.time()-t0})
    print(f'PASS_STAGE04B5_ALL_FORMAL_ACCESSION_SCAN_COMPLETE formal_bearing={len(done)} zero_formal={len(zero)} selected={selected_n}')

def concat_tsv(files,out):
    header=None
    with open(out,'w',newline='') as fo:
        for p in files:
            with open(p) as fi:
                h=fi.readline()
                if header is None: header=h; fo.write(h)
                elif h!=header: raise RuntimeError(f'header mismatch {p}')
                shutil.copyfileobj(fi,fo)

def run_mmseqs_cluster(mm,faa,prefix,tmp,threads):
    if tmp.exists(): shutil.rmtree(tmp)
    subprocess.run([mm,'easy-cluster',str(faa),str(prefix),str(tmp),'--min-seq-id','0.3','-c','0.7','--cov-mode','0','--cluster-mode','0','-s','7.5','--threads',str(threads)],check=True)
    c=Path(str(prefix)+'_cluster.tsv')
    if not c.exists(): raise RuntimeError(f'mmseqs cluster TSV missing {c}')
    mp={}
    with open(c) as f:
        for line in f:
            rep,mem=line.rstrip('\n').split('\t')[:2]; mp[mem]=rep
    return mp

def cluster_and_close(args):
    H=Path(args.home).resolve(); W=H/'work'; O=Path(args.outdir).resolve(); O.mkdir(parents=True,exist_ok=True)
    ad=W/'accession_results'; dirs=sorted(x for x in ad.iterdir() if x.is_dir())
    manifest=tsv_read(W/'ACCESSION_SCAN_MANIFEST.tsv')
    expected_acc={r['accession'] for r in manifest}
    observed_acc={d.name for d in dirs}
    if observed_acc!=expected_acc:
        miss=sorted(expected_acc-observed_acc); extra=sorted(observed_acc-expected_acc)
        raise SystemExit(f'HOLD formal-bearing accession results mismatch observed={len(observed_acc)} expected={len(expected_acc)} missing={miss[:5]} extra={extra[:5]}')
    zero=tsv_read(W/'ZERO_FORMAL_ACCESSIONS.tsv')
    if len(expected_acc)+len(zero)!=254:
        raise SystemExit(f'HOLD accession accounting formal-bearing={len(expected_acc)} zero={len(zero)} != 254')
    shutil.copy2(W/'ZERO_FORMAL_ACCESSIONS.tsv',O/'ZERO_FORMAL_ACCESSIONS.tsv')
    concat_tsv([d/'contexts.tsv' for d in dirs],O/'DELTA88487_CONTEXT_EVALUABILITY.tsv')
    concat_tsv([d/'arrays.tsv' for d in dirs],O/'DELTA88487_REPEAT_ARRAYS.tsv')
    concat_tsv([d/'partners.tsv' for d in dirs],O/'DELTA88487_NEAREST_DOWNSTREAM_ORFS.tsv')
    with open(O/'DELTA88487_NEAREST_DOWNSTREAM_ORFS.faa','w') as fo:
        for d in dirs:
            p=d/'partners.faa'
            if p.exists(): fo.write(p.read_text())
    ctx=tsv_read(O/'DELTA88487_CONTEXT_EVALUABILITY.tsv')
    if len(ctx)!=EXPECTED_FORMAL: raise SystemExit(f'HOLD concatenated context n={len(ctx)}')
    mm=shutil.which('mmseqs')
    if not mm: raise SystemExit('HOLD mmseqs not found')
    # candidate family clustering on exact-unique AA representatives
    fam_map_unique=run_mmseqs_cluster(mm,W/'FORMAL_UNIQUE_AA_85503.faa',O/'mmseqs_delta_family',O/'mmseqs_delta_family_tmp',args.threads)
    nfamily=len(set(fam_map_unique.values()))
    if nfamily!=EXPECTED_FAMILIES: raise SystemExit(f'HOLD reproduced delta family n={nfamily} != frozen {EXPECTED_FAMILIES}')
    cmap={r['candidate_id']:r for r in tsv_read(W/'FORMAL_CANDIDATE_CLASS_MAP.tsv')}
    cand_family={cid:fam_map_unique[r['unique_aa_rep_id']] for cid,r in cmap.items()}
    # partner clustering
    partners=tsv_read(O/'DELTA88487_NEAREST_DOWNSTREAM_ORFS.tsv')
    pmap={}
    if partners:
        pmap=run_mmseqs_cluster(mm,O/'DELTA88487_NEAREST_DOWNSTREAM_ORFS.faa',O/'mmseqs_partner',O/'mmseqs_partner_tmp',args.threads)
    pfields=list(partners[0].keys())+['partner_cluster'] if partners else ['partner_id','candidate_id','accession','partner_cluster']
    for r in partners: r['partner_cluster']=pmap.get(r['partner_id'],'UNCLUSTERED')
    tsv_write(O/'DELTA88487_NEAREST_DOWNSTREAM_ORFS_WITH_CLUSTER.tsv',partners,pfields)
    # candidate -> family/class map
    maprows=[]
    for cid in sorted(cmap):
        maprows.append(dict(candidate_id=cid,delta_operational_family_id=cand_family[cid],candidate_class=cmap[cid]['candidate_class'],combined_delta_only_family_id=cmap[cid]['combined_delta_only_family_id']))
    tsv_write(O/'DELTA88487_CANDIDATE_FAMILY_MAP.tsv',maprows,['candidate_id','delta_operational_family_id','candidate_class','combined_delta_only_family_id'])
    arrays=tsv_read(O/'DELTA88487_REPEAT_ARRAYS.tsv')
    arr_by=defaultdict(list)
    for r in arrays: arr_by[r['candidate_id']].append(r)
    part_by={r['candidate_id']:r for r in partners}
    ctx_by={r['candidate_id']:r for r in ctx}
    fam_members=defaultdict(list)
    for cid,fam in cand_family.items(): fam_members[fam].append(cid)
    summ=[]
    for fam in sorted(fam_members):
        cids=fam_members[fam]; classes=Counter(cmap[c]['candidate_class'] for c in cids)
        arm='DELTA_ONLY_NOVELTY_ARM' if classes['MIXED']==0 else 'MIXED_REPLICATION_ARM'
        accs={ctx_by[c]['accession'] for c in cids}
        strict_loci=[]; upstream_loci=[]; any_loci=[]
        for c in cids:
            aa=arr_by.get(c,[])
            if aa: any_loci.append(c)
            if any(x['orientation']=='UPSTREAM' for x in aa): upstream_loci.append(c)
            if any(x['strict_published_art_geometry']=='True' for x in aa): strict_loci.append(c)
        clusters=defaultdict(set)
        for c in strict_loci:
            pr=part_by.get(c)
            if pr and pr.get('partner_cluster') and pr['partner_cluster']!='UNCLUSTERED': clusters[pr['partner_cluster']].add(pr['accession'])
        recurrent={k:len(v) for k,v in clusters.items() if len(v)>=2}
        strict_eval=sum(ctx_by[c]['strict_array_evaluable']=='True' for c in cids)
        partner_eval=sum(ctx_by[c]['partner_evaluable']=='True' for c in cids)
        if strict_loci and recurrent: status='STRICT_PLUS_RECURRENT_PARTNER'
        elif strict_loci: status='STRICT_ART_GEOMETRY_NO_RECURRENT_PARTNER'
        elif strict_eval>0: status='EVALUABLE_NO_STRICT_ART_GEOMETRY'
        else: status='CONTEXT_CENSORED_NO_STRICT_CALL'
        dfams=sorted({cmap[c]['combined_delta_only_family_id'] for c in cids if cmap[c]['combined_delta_only_family_id']})
        summ.append(dict(delta_operational_family_id=fam,arm=arm,member_locus_n=len(cids),accession_n=len(accs),delta_only_member_n=classes['DELTA_ONLY'],mixed_member_n=classes['MIXED'],strict_array_evaluable_locus_n=strict_eval,partner_evaluable_locus_n=partner_eval,any_repeat_locus_n=len(set(any_loci)),upstream_repeat_locus_n=len(set(upstream_loci)),strict_art_geometry_locus_n=len(set(strict_loci)),recurrent_partner_cluster_n=len(recurrent),strict_plus_recurrent_partner=str(bool(strict_loci and recurrent)),architecture_status=status,combined_delta_only_family_ids=';'.join(dfams),recurrent_partner_clusters=';'.join(f'{k}:{v}' for k,v in sorted(recurrent.items()))))
    sf=list(summ[0].keys()); tsv_write(O/'DELTA26396_FAMILY_ARCHITECTURE_SUMMARY.tsv',summ,sf)
    # arm and evaluability summaries
    armrows=[]
    for arm in ['MIXED_REPLICATION_ARM','DELTA_ONLY_NOVELTY_ARM']:
        xs=[x for x in summ if x['arm']==arm]
        armrows.append(dict(arm=arm,family_n=len(xs),locus_n=sum(int(x['member_locus_n']) for x in xs),strict_array_evaluable_locus_n=sum(int(x['strict_array_evaluable_locus_n']) for x in xs),partner_evaluable_locus_n=sum(int(x['partner_evaluable_locus_n']) for x in xs),family_with_any_repeat_n=sum(int(x['any_repeat_locus_n'])>0 for x in xs),family_with_upstream_repeat_n=sum(int(x['upstream_repeat_locus_n'])>0 for x in xs),family_with_strict_art_geometry_n=sum(int(x['strict_art_geometry_locus_n'])>0 for x in xs),family_with_strict_plus_recurrent_partner_n=sum(x['strict_plus_recurrent_partner']=='True' for x in xs)))
    tsv_write(O/'ARCHITECTURE_BY_ARM.tsv',armrows,list(armrows[0].keys()))
    def ntrue(k): return sum(r[k]=='True' for r in ctx)
    evalsum={'formal_locus_n':len(ctx),'strict_array_evaluable_n':ntrue('strict_array_evaluable'),'partner_evaluable_n':ntrue('partner_evaluable'),'both_flanks_500bp_n':ntrue('both_flanks_500bp'),'both_flanks_2000bp_n':ntrue('both_flanks_2000bp'),'both_flanks_5000bp_n':ntrue('both_flanks_5000bp'),'both_flanks_10000bp_n':ntrue('both_flanks_10000bp')}
    write_json(O/'CONTEXT_EVALUABILITY_SUMMARY.json',evalsum)
    gate={'stage':'Stage05F3V2-Stage04B5','gate':'PASS_STAGE04B5_54G_MATCHED_ARCHITECTURE_REPLAY_COMPLETE','gate_pass':True,'selected_accession_n':254,'formal_bearing_accession_n':len(expected_acc),'zero_formal_accession_n':len(zero),'formal_locus_n':len(ctx),'exact_aa_unique_n':EXPECTED_UNIQUE,'delta_operational_family_n':nfamily,'delta_only_instance_n':sum(r['candidate_class']=='DELTA_ONLY' for r in ctx),'mixed_instance_n':sum(r['candidate_class']=='MIXED' for r in ctx),'context_evaluability':evalsum,'architecture_by_arm':armrows,'rules_changed':False,'retuning_performed':False,'source_reselection_performed':False,'structure_used_for_selection':False,'claim_boundary':'Architecture absence is interpretable only within evaluable loci. This stage scans all 88,487 frozen formal delta loci before any structural shortlist. Selected accessions with zero frozen formal RT-like candidates are retained as valid zero-results and require no sequence scan. Broad ART-like sensitivity is not redefined; strict published geometry plus recurrent partner is the matched primary endpoint.'}
    write_json(O/'STAGE04B5_CLOSEOUT_GATE.json',gate)
    # concise architecture-selected table
    sel=[x for x in summ if x['strict_plus_recurrent_partner']=='True']
    tsv_write(O/'ARCHITECTURE_SELECTED_FAMILIES.tsv',sel,sf)
    # hashes use final relative output names, not staging names
    with open(O/'SHA256SUMS.txt','w') as f:
        for p in sorted(x for x in O.iterdir() if x.is_file() and x.name!='SHA256SUMS.txt'):
            f.write(f'{sha(p)}  {p.name}\n')
    print(json.dumps(gate,indent=2)); print('PASS_STAGE04B5_FAMILY_AND_PARTNER_CLUSTERING_COMPLETE'); print('PASS_STAGE04B5_54G_MATCHED_ARCHITECTURE_REPLAY_COMPLETE')

def main():
    ap=argparse.ArgumentParser(); sub=ap.add_subparsers(dest='cmd',required=True)
    p=sub.add_parser('preflight'); p.add_argument('--project-root',required=True); p.add_argument('--home',required=True); p.add_argument('--formal-faa',default=None)
    s=sub.add_parser('scan'); s.add_argument('--home',required=True); s.add_argument('--workers',type=int,default=4)
    c=sub.add_parser('close'); c.add_argument('--home',required=True); c.add_argument('--outdir',required=True); c.add_argument('--threads',type=int,default=32)
    a=ap.parse_args()
    if a.cmd=='preflight': preflight(a)
    elif a.cmd=='scan': run_scan(a)
    else: cluster_and_close(a)
if __name__=='__main__': main()
