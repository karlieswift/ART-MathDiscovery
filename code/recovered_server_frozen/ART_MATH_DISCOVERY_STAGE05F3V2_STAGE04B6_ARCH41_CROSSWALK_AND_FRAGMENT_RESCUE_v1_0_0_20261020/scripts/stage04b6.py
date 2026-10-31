#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,os,re,shutil,subprocess,sys,time,urllib.request
from collections import defaultdict,Counter
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed

EXPECTED_DELTA_FORMAL=88487
EXPECTED_OLD_FORMAL=10598
EXPECTED_SELECTED=41
EXPECTED_MIXED=23
EXPECTED_DELTA=18
HIST7={'F3BI37356863','F3BI37262050','F3BI55915933','F3BI27373723','F3BI55368798','F3BI52178020','F3BI05882631'}
HIST5={'F3BI37356863','F3BI37262050','F3BI55915933','F3BI27373723','F3BI55368798'}
CODON={
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
COMP=str.maketrans('ACGTacgt','TGCAtgca')

def rc(s): return s.translate(COMP)[::-1].upper()
def tr(nt):
    nt=nt.upper(); return ''.join(CODON.get(nt[i:i+3],'X') for i in range(0,len(nt)-2,3))
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def tsv_read(p):
    with open(p,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def tsv_write(p,rows,fields):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
def write_json(p,x): Path(p).write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')
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
        except: pass
        if proc is not None:
            if proc.poll() is None:
                try: proc.terminate()
                except: pass
            try: proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill(); proc.wait()
def fasta_count(p): return sum(1 for n,s in iter_fasta(p))
def fasta_dict(p): return {n:s for n,s in iter_fasta(p)}

def parse_cid(cid):
    z=cid.split('|')
    if len(z)<7: raise ValueError(cid)
    acc,contig,strand,frame,coord,aalen,ordv=z[:7]
    m=re.fullmatch(r'(\d+)-(\d+)',coord); ma=re.fullmatch(r'aa(\d+)',aalen)
    if not m or not ma or strand not in '+-': raise ValueError(cid)
    a,b=map(int,m.groups()); aa=int(ma.group(1))
    return dict(accession=acc,contig_id=contig,strand=strand,start0=a,end0=b,aa_length=aa)

def discover_fasta(P,expected,prefer=()):
    cand=[]
    for root in P.iterdir():
        if not root.is_dir(): continue
        nm=root.name
        if 'ART_MATH_DISCOVERY' not in nm: continue
        for fp in root.glob('**/*'):
            if not fp.is_file() or fp.suffix.lower() not in ('.faa','.fa','.fasta'): continue
            low=(str(fp)).lower()
            if expected==EXPECTED_DELTA_FORMAL and 'f3v2' not in low: continue
            if expected==EXPECTED_OLD_FORMAL and ('f3v2' in low or 'delta' in low): continue
            if not any(k in low for k in ('formal','rvt','candidate','query')): continue
            try:n=fasta_count(fp)
            except: continue
            if n!=expected: continue
            score=sum(20 for k in prefer if k.lower() in low)+ (50 if 'formal' in low else 0)
            cand.append((score,str(fp),fp))
    if not cand: raise RuntimeError(f'HOLD no FASTA with {expected} records')
    return sorted(cand,reverse=True)[0][2].resolve()

def discover_selected_ledger(P):
    hs=list(P.glob('**/frozen/F3V2_DELTA_SELECTED_ACCESSIONS.tsv'))
    best=[]
    for h in hs:
        try:n=sum(1 for _ in open(h))-1
        except:continue
        best.append(((100 if n==254 else 0)+(1000 if 'DISTRIBUTED_ACQUISITION_R2' in str(h) else 0),str(h),h))
    if not best: raise RuntimeError('HOLD selected ledger not found')
    return sorted(best,reverse=True)[0][2].resolve()
def ledger_map(p):
    d={}
    for r in tsv_read(p):
        a=r.get('accession') or r.get('run_accession');
        if a:d[a]=r
    return d

def find_source_file(root,acc,row,ledger):
    for k,v in row.items():
        if v and ('path' in k.lower() or 'file' in k.lower()):
            q=Path(v); opts=[q] if q.is_absolute() else [root/q,ledger.parent/q]
            for x in opts:
                if x.exists() and x.is_file(): return x.resolve()
    for pat in [f'{acc}.contigs.fa.zst',f'{acc}.fa.zst',f'*{acc}*contigs*.zst']:
        hs=list(root.glob('**/'+pat))
        if hs:return hs[0].resolve()
    raise FileNotFoundError(acc)

def discover_old_family_map(P,old_ids):
    best=None
    dirs=[d for d in P.iterdir() if d.is_dir() and 'STAGE05F3B' in d.name and 'F3V2' not in d.name]
    for d in dirs:
      for fp in d.glob('**/*.tsv'):
        try:
            with open(fp) as f: hdr=f.readline().rstrip('\n').split('\t')
        except: continue
        ccols=[x for x in hdr if any(k in x.lower() for k in ('candidate_id','instance_id','query_id','member_id'))]
        fcols=[x for x in hdr if 'family' in x.lower()]
        if not ccols or not fcols: continue
        try: rows=tsv_read(fp)
        except: continue
        for cc in ccols:
          for fc in fcols:
            mp={}
            for r in rows:
                c=r.get(cc,''); fam=r.get(fc,'')
                if c in old_ids and fam.startswith('F3BI'): mp[c]=fam
            cov=len(mp)
            if cov and (best is None or cov>best[0]): best=(cov,fp,cc,fc,mp)
    if best is None: return {},None,None,None
    return best[4],best[1].resolve(),best[2],best[3]

def motif_hits(seq):
    return [(m.start(),m.group()) for m in re.finditer(r'[YF].DD',seq)]

def run_mmseqs_search(qfaa,tfaa,out,tmp,threads):
    mm=shutil.which('mmseqs')
    if not mm: raise RuntimeError('HOLD mmseqs not found')
    if Path(tmp).exists(): shutil.rmtree(tmp)
    cmd=[mm,'easy-search',str(qfaa),str(tfaa),str(out),str(tmp),'--min-seq-id','0.3','-c','0.7','--cov-mode','0','-s','7.5','--max-seqs','1000','--threads',str(threads),'--format-output','query,target,fident,qcov,tcov,evalue,bits']
    subprocess.run(cmd,check=True)

def prepare(args):
    H=Path(args.home).resolve(); P=Path(args.project_root).resolve(); O=Path(args.outdir).resolve(); W=H/'work'; W.mkdir(exist_ok=True); O.mkdir(parents=True,exist_ok=True)
    F=H/'frozen'; sel=tsv_read(F/'ARCHITECTURE_SELECTED_FAMILIES.tsv'); mem=tsv_read(F/'B5_SELECTED41_MEMBERS.tsv')
    if len(sel)!=EXPECTED_SELECTED: raise SystemExit(f'HOLD selected {len(sel)} !=41')
    if Counter(r['arm'] for r in sel)!=Counter({'MIXED_REPLICATION_ARM':23,'DELTA_ONLY_NOVELTY_ARM':18}): raise SystemExit('HOLD arm counts')
    df=discover_fasta(P,EXPECTED_DELTA_FORMAL,('stage03','delta'))
    daa=fasta_dict(df)
    by=defaultdict(list)
    for r in mem:
        cid=r['candidate_id']
        if cid not in daa: raise SystemExit(f'HOLD delta formal missing {cid}')
        z=dict(r); z['aa_sequence']=daa[cid]; z['motif_c_like_n']=len(motif_hits(daa[cid])); z['motif_c_like']=';'.join(x[1] for x in motif_hits(daa[cid])); by[r['delta_operational_family_id']].append(z)
    reps=[]
    for s in sel:
        fam=s['delta_operational_family_id']; xs=by[fam]
        # longest existing member, deterministic candidate-id tie-break
        best=sorted(xs,key=lambda r:(-int(r['aa_length']),r['candidate_id']))[0]
        reps.append({**s,'selected_candidate_id':best['candidate_id'],'selected_accession':best['accession'],'selected_aa_length':best['aa_length'],'selected_candidate_class':best['candidate_class'],'selected_motif_c_like_n':best['motif_c_like_n'],'selected_motif_c_like':best['motif_c_like']})
    rf=list(reps[0].keys()); tsv_write(O/'ARCH41_LONGEST_EXISTING_REPRESENTATIVES.tsv',reps,rf)
    mix=[r for r in reps if r['arm']=='MIXED_REPLICATION_ARM']; nov=[r for r in reps if r['arm']=='DELTA_ONLY_NOVELTY_ARM']
    # fastas
    with open(O/'MIXED23_LONGEST_EXISTING_REPRESENTATIVES.faa','w') as f:
        for r in mix:f.write(f">{r['delta_operational_family_id']}|selected={r['selected_candidate_id']}\n{daa[r['selected_candidate_id']]}\n")
    with open(O/'DELTA18_LONGEST_EXISTING_REPRESENTATIVES.faa','w') as f:
        for r in nov:f.write(f">{r['delta_operational_family_id']}|selected={r['selected_candidate_id']}\n{daa[r['selected_candidate_id']]}\n")
    # family-wide motif audit novel18
    ma=[]
    for s in nov:
        fam=s['delta_operational_family_id']
        xs=by[fam]
        ma.append(dict(delta_operational_family_id=fam,member_n=len(xs),max_aa_length=max(int(x['aa_length']) for x in xs),member_with_motif_c_like_n=sum(int(x['motif_c_like_n'])>0 for x in xs),all_motif_c_like=';'.join(sorted({m for x in xs for _,m in motif_hits(x['aa_sequence'])}))))
    tsv_write(O/'DELTA18_FAMILY_WIDE_SEQUENCE_AUDIT.tsv',ma,list(ma[0].keys()))
    # mixed historical crosswalk
    oldf=discover_fasta(P,EXPECTED_OLD_FORMAL,('f3b1','formal'))
    oldaa=fasta_dict(oldf); oldmap,mapsrc,ccol,fcol=discover_old_family_map(P,set(oldaa))
    qfaa=W/'mixed23_query.faa'
    with open(qfaa,'w') as f:
        for r in mix:f.write(f">{r['delta_operational_family_id']}\n{daa[r['selected_candidate_id']]}\n")
    m8=W/'mixed23_vs_old10598.m8'; run_mmseqs_search(qfaa,oldf,m8,W/'mmseqs_tmp',args.threads)
    hits=defaultdict(list)
    with open(m8) as f:
        for line in f:
            q,t,ident,qcov,tcov,e,bits=line.rstrip('\n').split('\t')[:7]
            fam=oldmap.get(t,'')
            hits[q].append(dict(target=t,old_family=fam,fident=float(ident),qcov=float(qcov),tcov=float(tcov),evalue=float(e),bits=float(bits)))
    cross=[]
    for r in mix:
        q=r['delta_operational_family_id']; hs=sorted(hits.get(q,[]),key=lambda x:(-x['bits'],-x['qcov'],-x['tcov'],x['target']))
        mapped=[x for x in hs if x['old_family']]
        fams=[]
        for x in mapped:
            if x['old_family'] not in fams:fams.append(x['old_family'])
        best=mapped[0] if mapped else (hs[0] if hs else None)
        any7=sorted(set(fams)&HIST7); any5=sorted(set(fams)&HIST5)
        cross.append(dict(delta_operational_family_id=q,selected_candidate_id=r['selected_candidate_id'],selected_aa_length=r['selected_aa_length'],direct_old_hit_n=len(hs),mapped_old_family_n=len(fams),best_old_family='' if best is None else best.get('old_family',''),best_old_target='' if best is None else best['target'],best_fident='' if best is None else f"{best['fident']:.6g}",best_qcov='' if best is None else f"{best['qcov']:.6g}",best_tcov='' if best is None else f"{best['tcov']:.6g}",historical7_overlap=';'.join(any7),historical5_overlap=';'.join(any5),all_mapped_old_families=';'.join(fams),crosswalk_status=('HISTORICAL5_REPLICATION' if any5 else 'HISTORICAL7_COMPARATOR_REPLICATION' if any7 else 'OLD_FAMILY_REPLICATION' if fams else 'DIRECT_OLD_HIT_UNMAPPED' if hs else 'NO_DIRECT_OLD_HIT')))
    tsv_write(O/'MIXED23_HISTORICAL_CROSSWALK.tsv',cross,list(cross[0].keys()))
    # unitig preflight for novel18
    seen={}; pref=[]
    for r in nov:
        acc=r['selected_accession']; url=f'https://s3.amazonaws.com/logan-pub/u/{acc}/{acc}.unitigs.fa.zst'
        if acc not in seen:
            status=0; size=0; err=''
            try:
                req=urllib.request.Request(url,method='HEAD')
                with urllib.request.urlopen(req,timeout=20) as resp:
                    status=getattr(resp,'status',200); size=int(resp.headers.get('Content-Length','0') or 0)
            except Exception as e: err=str(e)
            seen[acc]=(status,size,err)
        status,size,err=seen[acc]
        pref.append(dict(delta_operational_family_id=r['delta_operational_family_id'],selected_candidate_id=r['selected_candidate_id'],accession=acc,selected_aa_length=r['selected_aa_length'],unitig_url=url,remote_content_length_bytes=size,head_status=status,head_error=err))
    tsv_write(O/'DELTA18_UNITIG_PREFLIGHT.tsv',pref,list(pref[0].keys()))
    # output frozen rescue target manifest
    tsv_write(W/'DELTA18_RESCUE_TARGETS.tsv',[dict(delta_operational_family_id=r['delta_operational_family_id'],candidate_id=r['selected_candidate_id'],accession=r['selected_accession'],aa_length=r['selected_aa_length']) for r in nov],['delta_operational_family_id','candidate_id','accession','aa_length'])
    gate={'gate':'PASS_STAGE04B6_PREPARE_AND_MIXED23_CROSSWALK_COMPLETE','gate_pass':True,'target_family_n':41,'mixed_family_n':23,'delta_only_family_n':18,'delta_formal_faa':str(df),'old_formal_faa':str(oldf),'old_family_map_source':None if mapsrc is None else str(mapsrc),'old_family_map_coverage_n':len(oldmap),'mixed23_with_any_old_family_n':sum(bool(r['all_mapped_old_families']) for r in cross),'mixed23_historical7_overlap_n':sum(bool(r['historical7_overlap']) for r in cross),'mixed23_historical5_overlap_n':sum(bool(r['historical5_overlap']) for r in cross),'delta18_unique_selected_accession_n':len(set(r['selected_accession'] for r in nov)),'science_change':False,'retuning_performed':False}
    write_json(O/'STAGE04B6_PREPARE_GATE.json',gate); print(json.dumps(gate,indent=2)); print('PASS_STAGE04B6_PREPARE_AND_MIXED23_CROSSWALK_COMPLETE')

def download(args):
    H=Path(args.home).resolve(); O=Path(args.outdir).resolve(); U=Path(args.unitig_dir).resolve(); U.mkdir(parents=True,exist_ok=True)
    rows=tsv_read(O/'DELTA18_UNITIG_PREFLIGHT.tsv'); urls={r['accession']:r['unitig_url'] for r in rows}
    def one(item):
        acc,url=item; dst=U/f'{acc}.unitigs.fa.zst'
        if dst.exists() and dst.stat().st_size>0:return acc,str(dst),dst.stat().st_size,'EXISTS'
        if shutil.which('wget'):
            subprocess.run(['wget','-c',url,'-O',str(dst)],check=True)
        elif shutil.which('curl'):
            subprocess.run(['curl','-L','-C','-','-o',str(dst),url],check=True)
        else: raise RuntimeError('wget/curl unavailable')
        return acc,str(dst),dst.stat().st_size,'DOWNLOADED'
    out=[]
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        fut={ex.submit(one,it):it[0] for it in sorted(urls.items())}
        for i,f in enumerate(as_completed(fut),1):
            acc=fut[f]
            try:a,p,s,st=f.result()
            except Exception as e: raise SystemExit(f'HOLD unitig download {acc}: {e}')
            print(f'[PASS UNITIG {i}/{len(fut)}] {a} bytes={s} {st}',flush=True); out.append(dict(accession=a,local_unitig_file=p,bytes=s,status=st))
    tsv_write(H/'work/DELTA18_UNITIG_LOCAL_FILES.tsv',sorted(out,key=lambda r:r['accession']),['accession','local_unitig_file','bytes','status'])
    print('PASS_STAGE04B6_DELTA18_UNITIG_DOWNLOAD_COMPLETE')

def get_contig(sf,name):
    for n,s in iter_fasta(sf):
        if n==name:return s
    raise RuntimeError(f'contig not found {name}')

def recover_maximal(unitig,cand_nt):
    matches=[]
    for name,s in iter_fasta(unitig):
        start=0
        while True:
            i=s.find(cand_nt,start)
            if i<0:break
            matches.append((name,s,i,'+')); start=i+1
        rnt=rc(cand_nt); start=0
        while True:
            i=s.find(rnt,start)
            if i<0:break
            matches.append((name,rc(s),len(s)-(i+len(rnt)),'-')); start=i+1
    if len(matches)!=1:return None,len(matches)
    name,s,j,ori=matches[0]; L=len(cand_nt)
    # candidate is forward in oriented s at j
    left_stop=False; right_stop=False; a=j
    k=j-3
    while k>=0:
        if s[k:k+3] in ('TAA','TAG','TGA'):
            a=k+3; left_stop=True; break
        k-=3
    if not left_stop: a=j%3
    b=j+L; k=b
    while k+3<=len(s):
        if s[k:k+3] in ('TAA','TAG','TGA'):
            b=k; right_stop=True; break
        k+=3
    if not right_stop: b=len(s)-((len(s)-j)%3)
    nt=s[a:b]; aa=tr(nt)
    return dict(unitig_id=name,unitig_orientation=ori,match_start0=j,recovered_nt=nt,recovered_aa=aa,left_stop_bounded=left_stop,right_stop_bounded=right_stop,left_extension_aa=(j-a)//3,right_extension_aa=(b-(j+L))//3,recovered_aa_length=len(aa),unitig_length_bp=len(s)),1

def rescue(args):
    H=Path(args.home).resolve(); P=Path(args.project_root).resolve(); O=Path(args.outdir).resolve(); U=Path(args.unitig_dir).resolve()
    targets=tsv_read(H/'work/DELTA18_RESCUE_TARGETS.tsv'); df=discover_fasta(P,EXPECTED_DELTA_FORMAL,('stage03','delta')); daa=fasta_dict(df)
    ledger=discover_selected_ledger(P); lm=ledger_map(ledger); root=ledger.parents[1]
    sf_cache={}; contig_cache={}; out=[]; fa=[]
    for i,r in enumerate(targets,1):
        cid=r['candidate_id']; pc=parse_cid(cid); acc=pc['accession']; sf=sf_cache.get(acc)
        if sf is None: sf=find_source_file(root,acc,lm[acc],ledger); sf_cache[acc]=sf
        key=(str(sf),pc['contig_id'])
        if key not in contig_cache: contig_cache[key]=get_contig(sf,pc['contig_id'])
        genome=contig_cache[key]; nt=genome[pc['start0']:pc['end0']]; cand_nt=nt if pc['strand']=='+' else rc(nt); aa=tr(cand_nt).rstrip('*')
        if aa!=daa[cid].rstrip('*'): raise SystemExit(f'HOLD exact peptide replay {cid}')
        upath=U/f'{acc}.unitigs.fa.zst'
        if not upath.exists():
            out.append(dict(delta_operational_family_id=r['delta_operational_family_id'],candidate_id=cid,accession=acc,existing_aa_length=len(aa),unitig_match_n=0,rescue_status='UNITIG_FILE_MISSING',recovered_aa_length='',left_extension_aa='',right_extension_aa='',left_stop_bounded='',right_stop_bounded='',recovered_motif_c_like='',structure_interpretation='TECHNICAL_RESCUE_PENDING'))
            continue
        rec,nmatch=recover_maximal(upath,cand_nt)
        if rec is None:
            out.append(dict(delta_operational_family_id=r['delta_operational_family_id'],candidate_id=cid,accession=acc,existing_aa_length=len(aa),unitig_match_n=nmatch,rescue_status=('ABSENT' if nmatch==0 else 'AMBIGUOUS'),recovered_aa_length='',left_extension_aa='',right_extension_aa='',left_stop_bounded='',right_stop_bounded='',recovered_motif_c_like='',structure_interpretation='FRAGMENT_CONTEXT_UNRESOLVED'))
        else:
            # original peptide must occur in recovered aa at expected offset
            off=rec['left_extension_aa']
            if rec['recovered_aa'][off:off+len(aa)]!=aa: raise SystemExit(f'HOLD rescued peptide containment {cid}')
            mh=';'.join(m for _,m in motif_hits(rec['recovered_aa']))
            complete=rec['left_stop_bounded'] and rec['right_stop_bounded']
            status='UNIQUE_STOP_BOUNDED_BOTH' if complete else 'UNIQUE_EDGE_CENSORED'
            interp='STRUCTURE_READY_COMPLETE_STOP_TO_STOP' if complete else 'STRUCTURE_ALLOWED_BUT_NEGATIVE_NOT_AUTHORIZED'
            out.append(dict(delta_operational_family_id=r['delta_operational_family_id'],candidate_id=cid,accession=acc,existing_aa_length=len(aa),unitig_match_n=1,rescue_status=status,recovered_aa_length=rec['recovered_aa_length'],left_extension_aa=rec['left_extension_aa'],right_extension_aa=rec['right_extension_aa'],left_stop_bounded=rec['left_stop_bounded'],right_stop_bounded=rec['right_stop_bounded'],recovered_motif_c_like=mh,structure_interpretation=interp))
            fa.append((r['delta_operational_family_id'],cid,rec['recovered_aa']))
        print(f'[RESCUE {i}/{len(targets)}] {r["delta_operational_family_id"]}',flush=True)
    tsv_write(O/'DELTA18_UNITIG_RESCUE.tsv',out,list(out[0].keys()))
    with open(O/'DELTA18_RESCUED_MAXIMAL_STOP_TO_STOP.faa','w') as f:
        for fam,cid,s in fa:f.write(f'>{fam}|source_candidate={cid}\n{s}\n')
    write_json(O/'STAGE04B6_RESCUE_GATE.json',{'gate':'PASS_STAGE04B6_DELTA18_UNITIG_FRAGMENT_RESCUE_COMPLETE','gate_pass':True,'target_n':18,'unique_match_n':sum(int(r['unitig_match_n'])==1 for r in out),'stop_bounded_both_n':sum(r['rescue_status']=='UNIQUE_STOP_BOUNDED_BOTH' for r in out),'edge_censored_unique_n':sum(r['rescue_status']=='UNIQUE_EDGE_CENSORED' for r in out),'ambiguous_n':sum(r['rescue_status']=='AMBIGUOUS' for r in out),'absent_n':sum(r['rescue_status']=='ABSENT' for r in out),'science_change':False,'retuning_performed':False})
    print('PASS_STAGE04B6_DELTA18_UNITIG_FRAGMENT_RESCUE_COMPLETE')

def finalize(args):
    O=Path(args.outdir).resolve(); H=Path(args.home).resolve()
    pg=json.load(open(O/'STAGE04B6_PREPARE_GATE.json')); rg=json.load(open(O/'STAGE04B6_RESCUE_GATE.json')) if (O/'STAGE04B6_RESCUE_GATE.json').exists() else None
    if not pg.get('gate_pass'): raise SystemExit('HOLD prepare gate')
    if rg is None or not rg.get('gate_pass'): raise SystemExit('HOLD rescue gate missing')
    gate={'stage':'Stage05F3V2-Stage04B6','gate':'PASS_STAGE04B6_ARCH41_CROSSWALK_AND_FRAGMENT_RESCUE_COMPLETE','gate_pass':True,'mixed23_crosswalk_complete':True,'delta18_fragment_rescue_complete':True,'prepare':pg,'rescue':rg,'structure_run_performed':False,'rules_changed':False,'retuning_performed':False,'next_step':'Run frozen structural compatibility only on rescued delta-only representatives, while preserving fragment-censored negative interpretation rules.'}
    write_json(O/'STAGE04B6_CLOSEOUT_GATE.json',gate)
    with open(O/'SHA256SUMS.txt','w') as f:
        for p in sorted(x for x in O.iterdir() if x.is_file() and x.name!='SHA256SUMS.txt'): f.write(f'{sha(p)}  {p.name}\n')
    print(json.dumps(gate,indent=2)); print('PASS_STAGE04B6_ARCH41_CROSSWALK_AND_FRAGMENT_RESCUE_COMPLETE')

def main():
    ap=argparse.ArgumentParser(); sp=ap.add_subparsers(dest='cmd',required=True)
    p=sp.add_parser('prepare'); p.add_argument('--project-root',required=True); p.add_argument('--home',required=True); p.add_argument('--outdir',required=True); p.add_argument('--threads',type=int,default=32)
    d=sp.add_parser('download'); d.add_argument('--home',required=True); d.add_argument('--outdir',required=True); d.add_argument('--unitig-dir',required=True); d.add_argument('--workers',type=int,default=4)
    r=sp.add_parser('rescue'); r.add_argument('--project-root',required=True); r.add_argument('--home',required=True); r.add_argument('--outdir',required=True); r.add_argument('--unitig-dir',required=True)
    f=sp.add_parser('finalize'); f.add_argument('--home',required=True); f.add_argument('--outdir',required=True)
    a=ap.parse_args()
    if a.cmd=='prepare':prepare(a)
    elif a.cmd=='download':download(a)
    elif a.cmd=='rescue':rescue(a)
    else:finalize(a)
if __name__=='__main__':main()
