#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, os, re, shutil, subprocess, tarfile, tempfile
from collections import defaultdict
from pathlib import Path

EXPECTED_SOURCE_SHA = "5eec6d7ee00584c6266b8a953b6112540148b949341e400099b285ce30d235f0"
LOCAL_FLANK = 3000
MIN_AA = 25
LOCAL_E = 0.1
MAX_JOIN_AA = 900
RC = str.maketrans('ACGTN','TGCAN')
CODON_TABLE_11 = {
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

def sha256_file(p: Path):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20), b''): h.update(b)
    return h.hexdigest()

def open_text(p: Path):
    return gzip.open(p,'rt') if p.suffix=='.gz' else p.open('r', newline='')

def read_tsv(p: Path):
    with open_text(p) as f: return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(p: Path, rows, cols=None):
    rows=list(rows)
    if cols is None: cols=list(rows[0]) if rows else []
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(rows)

def read_fasta(path: Path):
    out={}; name=None; seq=[]
    with open_text(path) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if name is not None: out[name]=''.join(seq)
                name=line[1:].split()[0]; seq=[]
            else: seq.append(line.upper())
    if name is not None: out[name]=''.join(seq)
    return out

def read_fasta_one(path: Path):
    d=read_fasta(path)
    if not d: raise RuntimeError(f'no FASTA sequence in {path}')
    k=next(iter(d)); return k,d[k]

def materialize_input(src: Path, work: Path, label: str):
    if src.is_dir(): return src.resolve(), 'live_directory'
    if src.is_file() and src.name.endswith(('.tar.gz','.tgz')):
        out=work/f'unpacked_{label}'; out.mkdir(parents=True,exist_ok=True)
        with tarfile.open(src,'r:gz') as tf:
            for m in tf.getmembers():
                if not m.isfile(): continue
                dest=(out/m.name).resolve()
                if out.resolve() not in dest.parents and dest!=out.resolve(): raise RuntimeError('unsafe tar member')
            tf.extractall(out)
        # archive may unpack flat or nested
        candidates=[out]+[p for p in out.iterdir() if p.is_dir()]
        for p in candidates:
            if (p/'outputs').exists(): return p.resolve(), f'archive:{src.name}'
        return out.resolve(), f'archive:{src.name}'
    raise RuntimeError(f'unsupported input {src}')

def find_required(base: Path, name: str):
    p=base/'outputs'/name
    if p.exists(): return p
    hits=list(base.rglob(name))
    if len(hits)==1: return hits[0]
    raise RuntimeError(f'cannot resolve {name} under {base}; hits={len(hits)}')

def resolve_source_manifest(root: Path):
    override=os.environ.get('STAGE05F2B_SOURCE_MANIFEST','').strip()
    if override:
        p=Path(override)
        if sha256_file(p)!=EXPECTED_SOURCE_SHA: raise RuntimeError('source manifest override hash mismatch')
        return p.resolve(),'override'
    pats=['*SOURCE*MANIFEST*.tsv','*source*manifest*.tsv','*MANIFEST*.tsv','*manifest*.tsv']
    seen=set()
    for pat in pats:
        for p in root.rglob(pat):
            if not p.is_file(): continue
            rp=str(p.resolve())
            if rp in seen: continue
            seen.add(rp)
            try:
                if p.stat().st_size < 100_000_000 and sha256_file(p)==EXPECTED_SOURCE_SHA:
                    return p.resolve(),'live_exact_sha256'
            except: pass
    raise RuntimeError('could not resolve frozen Stage05B source manifest by exact SHA256')

def resolve_fasta_manifest(root: Path, manifest: Path):
    rows=read_tsv(manifest)
    if not rows: raise RuntimeError('empty source manifest')
    cols=list(rows[0])
    acc_col=next((x for x in ['sequence_id','accession','seq_id','genome_id'] if x in cols),None)
    path_col=next((x for x in cols if 'fasta' in x.lower() and ('path' in x.lower() or 'file' in x.lower() or 'frozen' in x.lower())),None)
    if not acc_col or not path_col: raise RuntimeError(f'cannot infer FASTA columns from {cols}')
    out={}
    for r in rows:
        acc=r.get(acc_col,'').strip(); q=r.get(path_col,'').strip()
        if not acc or not q: continue
        cands=[Path(q),manifest.parent/q,root/q]
        hit=next((p for p in cands if p.exists() and p.is_file()),None)
        if hit is None:
            hits=list(root.glob(f'**/{Path(q).name}'))[:5]
            hit=next((p for p in hits if p.exists() and p.is_file()),None)
        if hit: out[acc]=hit.resolve()
    return out

def translate(nt):
    return ''.join(CODON_TABLE_11.get(nt[i:i+3],'X') for i in range(0,len(nt)-2,3))

def segments(acc, seq, ws, we, min_aa=MIN_AA):
    sub=seq[ws:we]
    out=[]
    for strand in ['+','-']:
        nt=sub if strand=='+' else sub.translate(RC)[::-1]
        for frame in (0,1,2):
            aa=translate(nt[frame:]); pos=0
            parts=aa.split('*')
            for j,part in enumerate(parts):
                end=pos+len(part)
                if len(part)>=min_aa:
                    n0=frame+pos*3; n1=frame+end*3
                    if strand=='+': gs,ge=ws+n0,ws+n1
                    else: gs,ge=we-n1,we-n0
                    out.append({'sequence_id':acc,'start0':gs,'end0':ge,'strand':strand,'frame':frame,'aa_sequence':part,'aa_length':len(part),'left_edge_partial':pos==0,'right_edge_partial':j==len(parts)-1})
                pos=end+1
    return out

def parse_domtbl(path: Path):
    rows=[]
    if not path.exists(): return rows
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            x=line.rstrip('\n').split(maxsplit=22)
            if len(x)<23: continue
            try:
                rows.append({
                    'target_name':x[0],'target_len':int(x[2]),'query_name':x[3],'query_len':int(x[5]),
                    'full_evalue':float(x[6]),'full_score':float(x[7]),'domain_i_evalue':float(x[12]),'domain_score':float(x[13]),
                    'hmm_from':int(x[15]),'hmm_to':int(x[16]),'ali_from':int(x[17]),'ali_to':int(x[18]),'env_from':int(x[19]),'env_to':int(x[20]),'domain_acc':float(x[21]),'description':x[22]
                })
            except Exception: continue
    return rows

def union_len(intervals):
    arr=sorted((min(a,b),max(a,b)) for a,b in intervals)
    if not arr: return 0
    total=0; s,e=arr[0]
    for a,b in arr[1:]:
        if a<=e+1: e=max(e,b)
        else: total+=e-s+1; s,e=a,b
    return total+(e-s+1)

def boolv(x): return str(x).lower() in ('true','1','yes')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('f1_input'); ap.add_argument('f2a_input'); ap.add_argument('project_root'); ap.add_argument('--threads',type=int,default=16); ap.add_argument('--outdir',default='outputs')
    a=ap.parse_args(); root=Path(a.project_root).resolve(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True); work=Path('work'); work.mkdir(exist_ok=True)
    f1,f1route=materialize_input(Path(a.f1_input),work,'f1'); f2a,f2aroute=materialize_input(Path(a.f2a_input),work,'f2a')
    f1gate=json.loads(find_required(f1,'STAGE05F1_GATE.json').read_text())
    f2agate=json.loads(find_required(f2a,'STAGE05F2A_GATE.json').read_text())
    if not f1gate.get('gate_pass'): raise RuntimeError('F1 gate not PASS')
    if not f2agate.get('gate_pass'): raise RuntimeError('F2A gate not PASS')
    f1inst=read_tsv(find_required(f1,'STAGE05F1_SIXFRAME_PEPTIDE_INSTANCES.tsv'))
    f1pri=read_tsv(find_required(f1,'STAGE05F1_RESCUE_PRIORITY.tsv'))
    f1fa=read_fasta(find_required(f1,'STAGE05F1_UNIQUE_SIXFRAME_PEPTIDES.faa'))
    strict=read_tsv(find_required(f2a,'STAGE05F2A_STRICT_FORMAL_CALIBRATED.tsv'))
    # exact AA from query id
    aa_by_q=f1fa
    hash_to_q={r['aa_sha256']:r['query_id'] for r in f1pri}
    inst_by_hash=defaultdict(list)
    for r in f1inst: inst_by_hash[r['aa_sha256']].append(r)
    target_rows=[]
    for r in strict:
        prov=r.get('f2a_provenance_class','')
        queue=r.get('compute_queue','')
        selected=(prov=='KNOWN_CAPTURED_RT_CALIBRATION' or prov=='PRIOR_D1_FORMAL_CONTEXT' or (prov=='NOVELTY_FIRST_STRICT_RESCUE' and queue in {'P1_YXDD_PURE_MISSED_ORF','P2_CROSS_GENOME_REPLICATED','P3_HIGHCOV_SINGLETON'}))
        if not selected: continue
        role='NOVELTY_COMPUTE_PRIORITY' if prov=='NOVELTY_FIRST_STRICT_RESCUE' else ('KNOWN_RT_CALIBRATION' if prov=='KNOWN_CAPTURED_RT_CALIBRATION' else 'PRIOR_D1_BRIDGE')
        for inst in inst_by_hash.get(r['aa_sha256'],[]):
            z=dict(r); z.update({f'target_{k}':inst.get(k,'') for k in ['peptide_instance_id','sequence_id','start0','end0','strand','frame','edge_partial','prodigal_relation','candidate_keys']}); z['f2b_role']=role
            target_rows.append(z)
    if not target_rows: raise RuntimeError('no target instances resolved')
    # source FASTAs
    mani,maniroute=resolve_source_manifest(root); fmap=resolve_fasta_manifest(root,mani)
    missing=sorted({r['target_sequence_id'] for r in target_rows if r['target_sequence_id'] not in fmap})
    if missing: raise RuntimeError(f'missing genome FASTAs: {missing[:10]}')
    hmm=Path(f1gate['rvt_hmm'])
    if not hmm.exists():
        ov=os.environ.get('STAGE05F2B_RVT_HMM','').strip(); hmm=Path(ov) if ov else hmm
    if not hmm.exists(): raise RuntimeError('RVT HMM not found; set STAGE05F2B_RVT_HMM')
    # unique local windows per target instance (not merged, to preserve target-centric interpretation)
    local_frags=[]; tid_rows=[]; frag_counter=0
    for i,t in enumerate(target_rows,1):
        acc=t['target_sequence_id']; ts=int(t['target_start0']); te=int(t['target_end0']); strand=t['target_strand']
        _,seq=read_fasta_one(fmap[acc]); L=len(seq); ws=max(0,ts-LOCAL_FLANK); we=min(L,te+LOCAL_FLANK); tid=f'F2B_T{i:04d}'
        t['f2b_target_id']=tid; t['local_window_start0']=ws; t['local_window_end0']=we; t['genome_length']=L
        tid_rows.append(t)
        for p in segments(acc,seq,ws,we,MIN_AA):
            frag_counter+=1; fid=f'LF_{frag_counter:07d}'; aa=p.pop('aa_sequence'); h=hashlib.sha256(aa.encode()).hexdigest()
            row=dict(p); row.update({'f2b_target_id':tid,'fragment_id':fid,'aa_sha256':h,'aa_sequence':aa,'same_strand_as_target':p['strand']==strand,'overlaps_target':p['end0']>ts and p['start0']<te})
            local_frags.append(row)
    # dedupe AA for local HMM
    local_unique={}
    for r in local_frags: local_unique.setdefault(r['aa_sha256'],r['aa_sequence'])
    local_q={h:f'LQ_{i:07d}' for i,h in enumerate(sorted(local_unique),1)}
    lfa=work/'local_unique.faa'
    with lfa.open('w') as f:
        for h in sorted(local_unique): f.write(f'>{local_q[h]} sha256={h}\n{local_unique[h]}\n')
    ldom=out/'STAGE05F2B_LOCAL_RVT.domtblout'
    subprocess.run(['hmmscan','--cpu',str(a.threads),'-E',str(LOCAL_E),'--domE',str(LOCAL_E),'--domtblout',str(ldom),str(hmm),str(lfa)],check=True,stdout=(work/'local_hmmscan.stdout').open('w'),stderr=(work/'local_hmmscan.stderr').open('w'))
    dhits=parse_domtbl(ldom); hits_by_q=defaultdict(list)
    for h in dhits: hits_by_q[h['query_name']].append(h)
    frag_rows=[]
    for r in local_frags:
        q=local_q[r['aa_sha256']]; hs=sorted(hits_by_q.get(q,[]),key=lambda z:(z['full_evalue'],z['domain_i_evalue']))
        if not hs: continue
        for h in hs:
            z={k:v for k,v in r.items() if k!='aa_sequence'}; z.update(h); z['query_coverage']=(h['ali_to']-h['ali_from']+1)/max(1,h['query_len']); z['model_coverage']=(h['hmm_to']-h['hmm_from']+1)/max(1,h['target_len']); z['yxdd_like']=bool(re.search(r'[YF][A-Z]DD',local_unique[r['aa_sha256']]))
            frag_rows.append(z)
    write_tsv(out/'STAGE05F2B_LOCAL_FRAGMENT_HITS.tsv',frag_rows)
    # profile coverage / complementary fragments per target
    frag_by_tid=defaultdict(list)
    for r in frag_rows:
        if r['same_strand_as_target'] in (True,'True','true','1'): frag_by_tid[r['f2b_target_id']].append(r)
    cov_rows=[]; class_rows=[]
    target_by_tid={r['f2b_target_id']:r for r in tid_rows}
    for tid,t in target_by_tid.items():
        rows=frag_by_tid.get(tid,[]); by_model=defaultdict(list)
        for r in rows: by_model[r['target_name']].append(r)
        model_summ=[]
        for model,rr in by_model.items():
            tlen=max(int(x['target_len']) for x in rr); intervals=[(int(x['hmm_from']),int(x['hmm_to'])) for x in rr]
            union=union_len(intervals); best=max((int(x['hmm_to'])-int(x['hmm_from'])+1 for x in rr),default=0)
            # transcription-order monotonicity for fragments with non-overlapping genomic centers
            srt=sorted(rr,key=lambda x:(int(x['start0']),int(x['end0'])),reverse=(t['target_strand']=='-'))
            hmm_mids=[(int(x['hmm_from'])+int(x['hmm_to']))/2 for x in srt]
            mono=all(b>=a for a,b in zip(hmm_mids,hmm_mids[1:])) if len(hmm_mids)>1 else True
            gains=max(0,union-best)
            model_summ.append((union/max(1,tlen),gains/max(1,tlen),model,tlen,len(rr),mono,rr))
        model_summ.sort(reverse=True,key=lambda x:(x[0],x[1],x[4]))
        if model_summ:
            uc,gain,model,tlen,nfrag,mono,rr=model_summ[0]
            best_e=min(float(x['full_evalue']) for x in rr); yxdd=any(boolv(x['yxdd_like']) for x in rr)
        else:
            uc=gain=0.0; model=''; tlen=0; nfrag=0; mono=False; best_e=1.0; yxdd=False
        cov={'f2b_target_id':tid,'sequence_id':t['target_sequence_id'],'candidate_keys':t['target_candidate_keys'],'f2b_role':t['f2b_role'],'compute_queue':t.get('compute_queue',''),'target_query_id':t['query_id'],'target_aa_length':t['aa_length'],'target_rvt_profile':t['best_rvt_target'],'target_full_evalue':t['best_full_evalue'],'target_query_coverage':t['query_coverage'],'target_model_coverage':t['model_coverage'],'target_yxdd_like':t['yxdd_like'],'best_local_model':model,'best_local_model_len':tlen,'same_strand_local_hit_rows_n':len(rows),'best_model_fragment_n':nfrag,'best_model_union_coverage':uc,'best_model_incremental_coverage_over_single_fragment':gain,'profile_order_monotonic':mono,'any_local_yxdd_like':yxdd,'best_local_full_evalue':best_e}
        cov_rows.append(cov)
    write_tsv(out/'STAGE05F2B_LOCUS_PROFILE_COVERAGE.tsv',cov_rows)
    # one-stop same-frame joins within each target window: adjacent translated segments with exactly one codon stop between, same strand/frame
    # use all local fragments including <target> windows; join only neighboring genomic segments in translation order and total <=900 aa
    by_tid_sf=defaultdict(list)
    for r in local_frags:
        if not r['same_strand_as_target']: continue
        by_tid_sf[(r['f2b_target_id'],r['strand'],int(r['frame']))].append(r)
    joins=[]; joinseq={}; jc=0
    for key,rr in by_tid_sf.items():
        strand=key[1]; rr=sorted(rr,key=lambda x:int(x['start0']),reverse=(strand=='-'))
        for a1,a2 in zip(rr,rr[1:]):
            # adjacent stop-to-stop segments in same frame should be separated by ~3 nt in translation order; tolerate 3-6 nt due clipping arithmetic
            if strand=='+': gap=int(a2['start0'])-int(a1['end0'])
            else: gap=int(a1['start0'])-int(a2['end0'])
            if gap not in (3,6): continue
            aa=a1['aa_sequence']+'X'+a2['aa_sequence']
            if len(aa)>MAX_JOIN_AA: continue
            jc+=1; q=f'SJ_{jc:06d}'; joinseq[q]=aa
            joins.append({'join_id':q,'f2b_target_id':key[0],'sequence_id':a1['sequence_id'],'strand':strand,'frame':key[2],'left_fragment_id':a1['fragment_id'],'right_fragment_id':a2['fragment_id'],'genomic_gap_nt':gap,'joined_aa_length':len(aa),'contains_target_fragment':bool(a1['overlaps_target'] or a2['overlaps_target'])})
    sfa=work/'single_stop_joined.faa'
    with sfa.open('w') as f:
        for q in sorted(joinseq): f.write(f'>{q}\n{joinseq[q]}\n')
    jhits=[]
    if joinseq:
        jdom=out/'STAGE05F2B_SINGLE_STOP_RVT.domtblout'
        subprocess.run(['hmmscan','--cpu',str(a.threads),'-E',str(LOCAL_E),'--domE',str(LOCAL_E),'--domtblout',str(jdom),str(hmm),str(sfa)],check=True,stdout=(work/'join_hmmscan.stdout').open('w'),stderr=(work/'join_hmmscan.stderr').open('w'))
        hh=defaultdict(list)
        for h in parse_domtbl(jdom): hh[h['query_name']].append(h)
        for r in joins:
            for h in sorted(hh.get(r['join_id'],[]),key=lambda z:(z['full_evalue'],z['domain_i_evalue'])):
                z=dict(r); z.update(h); z['query_coverage']=(h['ali_to']-h['ali_from']+1)/max(1,h['query_len']); z['model_coverage']=(h['hmm_to']-h['hmm_from']+1)/max(1,h['target_len']); z['yxdd_like']=bool(re.search(r'[YF][A-Z]DD',joinseq[r['join_id']]))
                jhits.append(z)
    write_tsv(out/'STAGE05F2B_SINGLE_STOP_READTHROUGH_SENSITIVITY.tsv',jhits)
    # integrate one-stop evidence and classify
    join_by_tid=defaultdict(list)
    for r in jhits:
        if boolv(r['contains_target_fragment']): join_by_tid[r['f2b_target_id']].append(r)
    class_rows=[]
    for c in cov_rows:
        tid=c['f2b_target_id']; t=target_by_tid[tid]; jj=join_by_tid.get(tid,[])
        best_join=min(jj,key=lambda x:float(x['full_evalue'])) if jj else None
        join_improved=False
        if best_join:
            try:
                join_improved=(float(best_join['full_evalue']) < float(t['best_full_evalue'])/10.0) or (float(best_join['model_coverage']) >= float(t['model_coverage'])+0.15)
            except: pass
        gain=float(c['best_model_incremental_coverage_over_single_fragment'])
        union=float(c['best_model_union_coverage'])
        role=t['f2b_role']; q=t.get('compute_queue','')
        if role=='KNOWN_RT_CALIBRATION': cls='CALIBRATION_KNOWN_RT_RECOVERED'
        elif role=='PRIOR_D1_BRIDGE': cls='PRIOR_D1_CONTEXT_RECONSTRUCTION'
        elif gain>=0.15 and boolv(c['profile_order_monotonic']): cls='COMPLEMENTARY_SAME_PROFILE_FRAGMENT_CHAIN'
        elif join_improved: cls='SINGLE_STOP_BRIDGE_SUPPORT'
        elif int(t['aa_length'])>=180 and float(t['model_coverage'])>=0.40: cls='LONG_SINGLE_FRAGMENT_RT_LIKE'
        elif boolv(t['yxdd_like']) and t['strict_aa_status']=='PURE_NO_SAME_STRAND_OVERLAP_NONEDGE': cls='SHORT_YXDD_PROFILE_FRAGMENT_UNRESOLVED'
        elif q=='P2_CROSS_GENOME_REPLICATED': cls='REPLICATED_PROFILE_FRAGMENT_UNRESOLVED'
        else: cls='LOCAL_PROFILE_FRAGMENT_UNRESOLVED'
        z=dict(c); z.update({'f2b_locus_class':cls,'best_single_stop_profile':best_join['target_name'] if best_join else '','best_single_stop_full_evalue':best_join['full_evalue'] if best_join else '','best_single_stop_model_coverage':best_join['model_coverage'] if best_join else '','single_stop_improved_descriptively':join_improved})
        class_rows.append(z)
    write_tsv(out/'STAGE05F2B_LOCUS_CLASSIFICATION.tsv',class_rows)
    # follow-up queue novelty only, prioritized by stronger reconstruction evidence; no biological deletion
    priority_order={'COMPLEMENTARY_SAME_PROFILE_FRAGMENT_CHAIN':1,'SINGLE_STOP_BRIDGE_SUPPORT':2,'LONG_SINGLE_FRAGMENT_RT_LIKE':3,'SHORT_YXDD_PROFILE_FRAGMENT_UNRESOLVED':4,'REPLICATED_PROFILE_FRAGMENT_UNRESOLVED':5,'LOCAL_PROFILE_FRAGMENT_UNRESOLVED':6}
    follow=[]
    for r in class_rows:
        if r['f2b_role']!='NOVELTY_COMPUTE_PRIORITY': continue
        z=dict(r); z['f2c_priority']=priority_order.get(r['f2b_locus_class'],9); z['recommended_next']='FULL_LOCUS_ORF_RECONSTRUCTION_AND_CALIBRATED_STRUCTURE' if z['f2c_priority']<=4 else 'PROFILE_PROFILE_AND_CROSS_GENOME_CONTEXT_BEFORE_STRUCTURE'
        follow.append(z)
    follow.sort(key=lambda z:(int(z['f2c_priority']),-float(z['best_model_union_coverage']),float(z['best_local_full_evalue'])))
    write_tsv(out/'STAGE05F2B_F2C_FOLLOWUP_QUEUE.tsv',follow)
    write_tsv(out/'STAGE05F2B_TARGET_INSTANCES.tsv',tid_rows)
    # FASTA of target novelty AA for next stage
    qseen=set()
    with (out/'STAGE05F2B_F2C_TARGETS.faa').open('w') as f:
        for r in follow:
            qid=r['target_query_id']
            if qid in qseen or qid not in aa_by_q: continue
            qseen.add(qid); f.write(f'>{qid} class={r["f2b_locus_class"]} queue={r["compute_queue"]}\n{aa_by_q[qid]}\n')
    counts=defaultdict(int)
    for r in class_rows: counts[r['f2b_locus_class']]+=1
    novelty_counts=defaultdict(int)
    for r in class_rows:
        if r['f2b_role']=='NOVELTY_COMPUTE_PRIORITY': novelty_counts[r['f2b_locus_class']]+=1
    gate={'stage':'Stage05F2B','version':'v1.0.0','gate':'PASS_STAGE05F2B_TARGETED_LOCUS_FRAGMENT_RECONSTRUCTION','gate_pass':True,'f1_input_route':f1route,'f2a_input_route':f2aroute,'f1_gate':f1gate.get('gate'),'f2a_gate':f2agate.get('gate'),'source_manifest_sha256':sha256_file(mani),'source_manifest_route':maniroute,'target_instance_n':len(tid_rows),'target_unique_aa_n':len({r['aa_sha256'] for r in tid_rows}),'novelty_compute_priority_unique_aa_n':len({r['aa_sha256'] for r in tid_rows if r['f2b_role']=='NOVELTY_COMPUTE_PRIORITY'}),'novelty_compute_priority_instance_n':sum(1 for r in tid_rows if r['f2b_role']=='NOVELTY_COMPUTE_PRIORITY'),'known_rt_calibration_instance_n':sum(1 for r in tid_rows if r['f2b_role']=='KNOWN_RT_CALIBRATION'),'prior_d1_bridge_instance_n':sum(1 for r in tid_rows if r['f2b_role']=='PRIOR_D1_BRIDGE'),'local_unique_aa_n':len(local_unique),'local_hmm_hit_rows_n':len(frag_rows),'single_stop_join_n':len(joins),'single_stop_hmm_hit_rows_n':len(jhits),'locus_class_counts':dict(counts),'novelty_locus_class_counts':dict(novelty_counts),'local_flank_bp':LOCAL_FLANK,'local_min_aa':MIN_AA,'local_hmmscan_threshold':LOCAL_E,'dna_score_recomputed':False,'dna_rank_recomputed':False,'candidate_membership_changed':False,'claim_boundary':'Targeted locus-level reconstruction only. Profile-fragment chains and single-stop sensitivity are hypotheses for F2C, not biochemical RT or ART proof.'}
    (out/'STAGE05F2B_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    notes=['# Stage05F2B manuscript-facing notes\n\n','F2B performs targeted local reconstruction around the frozen F2A P1/P2/P3 novelty set, with known RT loci and prior D1 contexts retained only as calibration/bridge rows. It does not rescore DNA or change candidate membership.\n\n','## Novelty locus classes\n']
    for k,v in sorted(novelty_counts.items()): notes.append(f'- {k}: {v}\n')
    notes.append('\nA complementary same-profile fragment chain means multiple same-strand local translation fragments collectively cover additional portions of one RVT HMM in transcription-consistent order. A single-stop bridge is a targeted sensitivity across one stop codon. Neither is evidence of actual translation through the stop or of biochemical RT activity.\n')
    (out/'STAGE05F2B_MANUSCRIPT_NOTES.md').write_text(''.join(notes))
    print(json.dumps(gate,indent=2))

if __name__=='__main__': main()
