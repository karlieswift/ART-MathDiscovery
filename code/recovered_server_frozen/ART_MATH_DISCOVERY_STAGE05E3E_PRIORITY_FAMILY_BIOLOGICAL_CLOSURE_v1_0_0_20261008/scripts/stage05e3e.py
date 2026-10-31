#!/usr/bin/env python3
import argparse,csv,hashlib,json,math,os,re,shutil,statistics,subprocess,tarfile,tempfile
from pathlib import Path

STAGE='Stage05E3E'; VERSION='v1.0.0'
PASS_E3D='PASS_E3D_FAMILY_DEFINING_PREDICTED_WITNESS_SEQUENCES_HASH_MATERIALIZED'


def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def norm_aa(s): return ''.join(str(s).split()).upper().rstrip('*')
def revcomp(s):
    return s.translate(str.maketrans('ACGTNacgtn','TGCANtgcan'))[::-1].upper()
def truthy(v): return str(v).strip().lower() in {'1','true','t','yes','y'}
def safe_extract(src,dst):
    dst=Path(dst).resolve(); dst.mkdir(parents=True,exist_ok=True)
    with tarfile.open(src,'r:*') as tf:
        for m in tf.getmembers():
            t=(dst/m.name).resolve()
            if not str(t).startswith(str(dst)+os.sep) and t!=dst: raise RuntimeError(f'unsafe archive member: {m.name}')
        tf.extractall(dst)
def resolve_input(src,root):
    p=Path(src)
    if p.exists(): return p
    hits=[x for x in Path(root).rglob(p.name) if x.is_file() or x.is_dir()]
    if len(hits)==1:
        print(f'[RESOLVE] {p} -> {hits[0]}'); return hits[0]
    if not hits: raise FileNotFoundError(f'{p}; exact-basename search under {root} found 0')
    # Permit byte-identical file copies only.
    files=[x for x in hits if x.is_file()]
    if len(files)==len(hits) and files:
        hs={sha256_file(x) for x in files}
        if len(hs)==1:
            pref=sorted(files,key=lambda x:(0 if '/outputs/' in str(x) else 1,len(str(x)),str(x)))[0]
            print(f'[RESOLVE] {len(files)} byte-identical copies -> {pref}'); return pref
    raise RuntimeError(f'ambiguous exact-basename resolution for {p.name}: {hits[:20]}')
def materialize(src,dst,root):
    src=resolve_input(src,root); dst=Path(dst)
    if dst.exists(): shutil.rmtree(dst)
    if src.is_dir(): shutil.copytree(src,dst)
    else: safe_extract(src,dst)
    return dst

def find_one(root,name):
    hits=list(Path(root).rglob(name))
    if len(hits)==1: return hits[0]
    if hits:
        files=[x for x in hits if x.is_file()]
        if len(files)==len(hits):
            hs={sha256_file(x) for x in files}
            if len(hs)==1: return sorted(files,key=lambda x:(0 if '/outputs/' in str(x) else 1,str(x)))[0]
    raise RuntimeError(f'expected one/byte-identical {name}; found {len(hits)}: {hits[:12]}')
def read_tsv(p):
    with open(p,'r',encoding='utf-8',errors='replace',newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(p,rows,fields):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,'w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore',lineterminator='\n'); w.writeheader(); w.writerows(rows)
def read_fasta(p):
    out=[]; h=None; seq=[]
    with open(p,'r',encoding='utf-8',errors='replace') as f:
        for line in f:
            line=line.rstrip('\n\r')
            if line.startswith('>'):
                if h is not None: out.append((h,''.join(seq).replace(' ','').upper()))
                h=line[1:].strip(); seq=[]
            elif h is not None: seq.append(line.strip())
    if h is not None: out.append((h,''.join(seq).replace(' ','').upper()))
    return out
def write_fasta(p,recs):
    with open(p,'w') as f:
        for h,s in recs:
            f.write('>'+h+'\n')
            for i in range(0,len(s),80): f.write(s[i:i+80]+'\n')
def first_token(h): return h.split()[0]
def parse_tblout(p):
    rows=[]
    with open(p,'r',encoding='utf-8',errors='replace') as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            parts=line.rstrip('\n').split(maxsplit=18)
            if len(parts)<18: continue
            rows.append({
                'target_name':parts[0], 'target_accession':parts[1], 'query_name':parts[2], 'query_accession':parts[3],
                'full_evalue':float(parts[4]), 'full_score':float(parts[5]), 'full_bias':float(parts[6]),
                'best_domain_evalue':float(parts[7]), 'best_domain_score':float(parts[8]), 'best_domain_bias':float(parts[9]),
                'description':parts[18] if len(parts)>18 else ''
            })
    return rows

def run_phmmer(query,db,outprefix,cpu,evalue):
    exe=shutil.which('phmmer')
    if not exe: raise RuntimeError('phmmer not found; activate the project art_rt_evidence environment (HMMER 3.4)')
    tbl=Path(str(outprefix)+'.tblout'); dom=Path(str(outprefix)+'.domtblout'); stdout=Path(str(outprefix)+'.stdout.txt')
    cmd=[exe,'--cpu',str(cpu),'--noali','-E',str(evalue),'--domE',str(evalue),'--incE',str(evalue),'--incdomE',str(evalue),'--tblout',str(tbl),'--domtblout',str(dom),str(query),str(db)]
    with open(stdout,'w') as fo:
        cp=subprocess.run(cmd,stdout=fo,stderr=subprocess.STDOUT,text=True)
    if cp.returncode!=0: raise RuntimeError(f'phmmer failed rc={cp.returncode}; inspect {stdout}')
    return tbl,dom,stdout,cmd

def parse_json(p): return json.loads(Path(p).read_text())
def get_policy(pkgroot): return json.loads((pkgroot/'config/POLICY.json').read_text())
def load_priority(pkgroot): return read_tsv(pkgroot/'config/PRIORITY_FAMILIES.tsv')

def resolve_proteome(e3d,project_root):
    manifest=read_tsv(find_one(e3d,'STAGE05E3D_SEQUENCE_MANIFEST.tsv'))
    candidates=[]
    for r in manifest:
        if r.get('sequence_type')=='family_defining_predicted_witness_aa' and r.get('source_file'):
            candidates.append(Path(r['source_file']))
    for p in candidates:
        if p.exists(): return p
    names={p.name for p in candidates}
    if not names: names={'FULL919_PREDICTED_PROTEINS.faa'}
    hits=[]
    for n in names: hits.extend(Path(project_root).rglob(n))
    hits=[x for x in hits if x.is_file()]
    if not hits: raise FileNotFoundError(f'could not resolve full predicted proteome from E3D source paths / {names}')
    hs={sha256_file(x) for x in hits}
    if len(hs)>1: raise RuntimeError(f'multiple non-identical proteome candidates: {hits[:10]}')
    return sorted(hits,key=lambda x:(0 if 'STAGE05C1' in str(x) else 1,len(str(x)),str(x)))[0]

def load_fasta_index(p):
    return {first_token(h):s for h,s in read_fasta(p)}

def accession_from_record(record):
    return record.split('__',1)[0] if '__' in record else record.split('|',1)[0]

def aliases(header): return [x.strip() for x in header]

def discover_coordinate_registry(stage_dir,target_ids):
    id_alias=['record_key','protein_key','orf_id','orf_key','protein_id','predicted_protein_id','target_name']
    acc_alias=['sequence_id','accession','acc','genome_id']
    start_alias=['start0','orf_start0','cds_start0','gene_start0','start']
    end_alias=['end0','orf_end0','cds_end0','gene_end0','end']
    strand_alias=['strand','orf_strand','cds_strand','gene_strand']
    candidates=[]
    for p in Path(stage_dir).rglob('*.tsv'):
        try:
            with open(p,'r',encoding='utf-8',errors='replace') as f: hdr=f.readline().rstrip('\n\r').split('\t')
        except Exception: continue
        def choose(opts):
            for x in opts:
                if x in hdr: return x
            return None
        ic,ac,sc,ec,tc=choose(id_alias),choose(acc_alias),choose(start_alias),choose(end_alias),choose(strand_alias)
        if ic and ac and sc and ec and tc:
            candidates.append((p,ic,ac,sc,ec,tc))
    scored=[]
    sample=set(list(target_ids)[:1000])
    for p,ic,ac,sc,ec,tc in candidates:
        n=0; total=0
        try:
            with open(p,'r',encoding='utf-8',errors='replace',newline='') as f:
                rr=csv.DictReader(f,delimiter='\t')
                for r in rr:
                    total+=1
                    if r.get(ic,'') in sample: n+=1
        except Exception: continue
        if n: scored.append((n,total,p,ic,ac,sc,ec,tc))
    if not scored: return None,None
    scored.sort(key=lambda x:(-x[0],-x[1],str(x[2])))
    best=scored[0]
    if len(scored)>1 and scored[1][0]==best[0] and scored[1][1]==best[1] and sha256_file(scored[1][2])!=sha256_file(best[2]):
        # same coverage and size but non-identical: still use the one in outputs if unique preference; record diagnostics.
        pass
    p,ic,ac,sc,ec,tc=best[2:]
    idx={}
    with open(p,'r',encoding='utf-8',errors='replace',newline='') as f:
        rr=csv.DictReader(f,delimiter='\t')
        for r in rr:
            rid=r.get(ic,'')
            if not rid: continue
            try: start=int(float(r.get(sc,''))); end=int(float(r.get(ec,'')))
            except: continue
            idx[rid]={'sequence_id':r.get(ac,'') or accession_from_record(rid),'start0':start,'end0':end,'strand':r.get(tc,''),'registry_path':str(p),'id_column':ic}
    return idx,{'path':str(p),'sha256':sha256_file(p),'id_col':ic,'acc_col':ac,'start_col':sc,'end_col':ec,'strand_col':tc,'matched_target_sample_n':best[0],'rows_n':best[1]}

def resolve_genome_fasta(project_root,acc):
    names=[f'{acc}.frozen.fasta',f'{acc}.fasta',f'{acc}.fa',f'{acc}.fna']
    hits=[]
    for n in names: hits.extend(Path(project_root).rglob(n))
    hits=[x for x in hits if x.is_file() and 'E3E' not in str(x)]
    if not hits: return None
    # prefer Stage05A frozen input.
    hits=sorted(hits,key=lambda x:(0 if ('STAGE05A' in str(x) and '.frozen.fasta' in x.name) else 1, len(str(x)), str(x)))
    # If alternatives disagree, use frozen Stage05A if present; otherwise fail to avoid accidental different record.
    if not ('STAGE05A' in str(hits[0]) and '.frozen.fasta' in hits[0].name):
        hs={sha256_file(x) for x in hits}
        if len(hs)>1: return None
    return hits[0]

def read_single_genome(p):
    recs=read_fasta(p)
    if len(recs)!=1: raise RuntimeError(f'expected single-record genome FASTA: {p}, found {len(recs)}')
    return recs[0][1].upper()
def scan_motif(seq,motif):
    motif=motif.upper(); rc=revcomp(motif); hits=[]
    for label,m in [('+',motif),('-',rc)]:
        start=0
        while True:
            i=seq.find(m,start)
            if i<0: break
            hits.append((i,i+len(m),label)); start=i+1
    # dedupe palindromic duplicates
    seen=set(); out=[]
    for h in sorted(hits):
        k=(h[0],h[1],h[2] if motif!=rc else '.')
        if k not in seen: seen.add(k); out.append(h)
    return out

def topology_and_gap(gs,ge,strand,rs,re):
    if re<=gs:
        side='UPSTREAM' if strand=='+' else 'DOWNSTREAM'; gap=gs-re
    elif rs>=ge:
        side='DOWNSTREAM' if strand=='+' else 'UPSTREAM'; gap=rs-ge
    else: side='OVERLAP'; gap=0
    return side,gap

def global_identity(a,b):
    # Needleman-Wunsch with linear gap; descriptive only.
    a=norm_aa(a); b=norm_aa(b); n=len(a); m=len(b)
    match=2; mismatch=-1; gap=-2
    prev=[j*gap for j in range(m+1)]; ptr=[]
    # Keep traceback rows as bytearray directions 0 diag,1 up,2 left. Sizes are small for phmmer hits.
    dirs=[bytearray(m+1)]
    for i in range(1,n+1):
        cur=[i*gap]+[0]*m; drow=bytearray(m+1); drow[0]=1
        ai=a[i-1]
        for j in range(1,m+1):
            vals=(prev[j-1]+(match if ai==b[j-1] else mismatch), prev[j]+gap, cur[j-1]+gap)
            k=max(range(3),key=lambda x:vals[x]); cur[j]=vals[k]; drow[j]=k
        prev=cur; dirs.append(drow)
    i,j=n,m; matches=0; aln=0; q_res=0; t_res=0
    while i>0 or j>0:
        if i>0 and j>0 and dirs[i][j]==0:
            aln+=1; q_res+=1; t_res+=1; matches+= (a[i-1]==b[j-1]); i-=1; j-=1
        elif i>0 and (j==0 or dirs[i][j]==1): aln+=1; q_res+=1; i-=1
        else: aln+=1; t_res+=1; j-=1
    ident=matches/aln if aln else 0.0
    return ident, (q_res/n if n else 0.0), (t_res/m if m else 0.0), aln

def resolve_rt_loci(project_root):
    hits=list(Path(project_root).rglob('RT_LOCI_ALL919.tsv'))
    if not hits: return None,None
    if len(hits)>1:
        hs={sha256_file(x) for x in hits}
        if len(hs)>1: return None,{'status':'AMBIGUOUS_NONIDENTICAL','paths':[str(x) for x in hits[:20]]}
    p=sorted(hits,key=lambda x:(0 if '/outputs/' in str(x) else 1,str(x)))[0]
    rows=read_tsv(p); idx={}
    acc_alias=['sequence_id','accession','acc']
    start_alias=['orf_start0','start0','cds_start0','rt_start0','start']
    end_alias=['orf_end0','end0','cds_end0','rt_end0','end']
    def col(opts):
        if not rows: return None
        for x in opts:
            if x in rows[0]: return x
        return None
    ac,sc,ec=col(acc_alias),col(start_alias),col(end_alias)
    if not ac or not sc or not ec: return None,{'status':'SCHEMA_UNRESOLVED','path':str(p),'columns':list(rows[0]) if rows else []}
    for r in rows:
        try: s=int(float(r[sc])); e=int(float(r[ec]))
        except: continue
        idx.setdefault(r[ac],[]).append((s,e,r))
    return idx,{'status':'RESOLVED','path':str(p),'sha256':sha256_file(p),'rows_n':len(rows),'acc_col':ac,'start_col':sc,'end_col':ec}
def interval_distance(a0,a1,b0,b1):
    if a1<=b0: return b0-a1
    if b1<=a0: return a0-b1
    return 0

def median(xs): return statistics.median(xs) if xs else ''

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--e3d',required=True); ap.add_argument('--project-root',required=True); ap.add_argument('--outdir',default='outputs'); ap.add_argument('--workdir',default='work'); ap.add_argument('--cpu',type=int,default=16)
    a=ap.parse_args(); project=Path(a.project_root); pkgroot=Path(__file__).resolve().parents[1]
    out=Path(a.outdir); work=Path(a.workdir)
    for p in (out,work):
        if p.exists(): shutil.rmtree(p)
        p.mkdir(parents=True)
    policy=get_policy(pkgroot); pri=load_priority(pkgroot); priority=[r['family_id'] for r in pri]; pmap={r['family_id']:r for r in pri}
    e3d=materialize(a.e3d,work/'e3d',project)
    gate=parse_json(find_one(e3d,'STAGE05E3D_GATE.json'))
    if not gate.get('gate_pass') or gate.get('gate')!=PASS_E3D: raise RuntimeError(f'E3D gate not exact PASS: {gate}')
    famsum=read_tsv(find_one(e3d,'STAGE05E3D_FAMILY_SEQUENCE_SUMMARY.tsv'))
    faudit=read_tsv(find_one(e3d,'STAGE05E3D_PREDICTED_WITNESS_AUDIT.tsv'))
    famsum_map={r['family_id']:r for r in famsum}
    locks=[]; queries={}; motifs={}; witness_ids={}
    qdir=out/'queries'; qdir.mkdir()
    for fam in priority:
        if fam not in famsum_map: raise RuntimeError(f'missing priority family {fam} in E3D summary')
        s=famsum_map[fam]; exp=pmap[fam]
        rows=[r for r in faudit if r.get('family_id')==fam]
        clusters=sorted({r.get('signature_FAM50_cluster','') for r in rows if r.get('signature_FAM50_cluster','')})
        if clusters != [exp['fam50_cluster']]: raise RuntimeError(f'{fam} FAM50 mismatch: observed={clusters}, expected={exp["fam50_cluster"]}')
        faa=find_one(e3d,f'{fam}.predicted_witness.faa'); fna=find_one(e3d,f'{fam}.repeat_core.fna')
        arecs=read_fasta(faa); rrecs=read_fasta(fna)
        unique_aa={norm_aa(x[1]) for x in arecs}; unique_rep={x[1].upper() for x in rrecs}
        obs_w=len(arecs); obs_u=len(unique_aa); obs_r=len(unique_rep)
        if obs_w!=int(exp['expected_official_supported_witness_n']) or obs_u!=1 or obs_r!=1:
            raise RuntimeError(f'{fam} lock mismatch: witness={obs_w}, uniqueAA={obs_u}, repeat={obs_r}')
        qseq=next(iter(unique_aa)); motif=next(iter(unique_rep)); queries[fam]=qseq; motifs[fam]=motif; witness_ids[fam]={r['record_key'] for r in rows}
        write_fasta(qdir/f'{fam}.representative.faa',[(f'{fam}|{exp["fam50_cluster"]}|E3D_REP',qseq)])
        (qdir/f'{fam}.repeat_core.fna').write_text(f'>{fam}|E3D_REPEAT_CORE\n{motif}\n')
        locks.append({'family_id':fam,'fam50_cluster':exp['fam50_cluster'],'witness_n':obs_w,'unique_witness_aa_n':obs_u,'representative_aa_length':len(qseq),'representative_aa_sha256':hashlib.sha256(qseq.encode()).hexdigest(),'repeat_core_unique_n':obs_r,'repeat_core_length':len(motif),'repeat_core_sha256':hashlib.sha256(motif.encode()).hexdigest(),'frozen_e3b_genomes':exp['frozen_e3b_genomes'],'frozen_e3b_taxonomy_tail_groups':exp['frozen_e3b_taxonomy_tail_groups'],'frozen_e3b_modal_side':exp['frozen_e3b_modal_side'],'frozen_e3b_median_gap_bp':exp['frozen_e3b_median_gap_bp']})
    write_tsv(out/'STAGE05E3E_PRIORITY_FAMILY_LOCK.tsv',locks,list(locks[0]))

    proteome=resolve_proteome(e3d,project); proteome_sha=sha256_file(proteome); pidx=load_fasta_index(proteome)
    stage_dir=proteome.parent.parent if proteome.parent.name=='outputs' else proteome.parent
    # Run phmmer separately to preserve per-family output.
    hitrows=[]; commands=[]
    for fam in priority:
        tbl,dom,stdout,cmd=run_phmmer(qdir/f'{fam}.representative.faa',proteome,work/f'{fam}.phmmer',a.cpu,policy['phmmer_report_evalue'])
        commands.append({'family_id':fam,'command':' '.join(cmd),'tblout':str(tbl),'domtblout':str(dom),'stdout':str(stdout)})
        for h in parse_tblout(tbl):
            target=h['target_name']; tseq=pidx.get(target,'')
            if not tseq: continue
            ident,qcov,tcov,alen=global_identity(queries[fam],tseq)
            hitrows.append({'family_id':fam,'fam50_cluster':pmap[fam]['fam50_cluster'],'target_name':target,'sequence_id':accession_from_record(target),'full_evalue':h['full_evalue'],'full_score':h['full_score'],'best_domain_evalue':h['best_domain_evalue'],'best_domain_score':h['best_domain_score'],'target_aa_length':len(tseq),'global_alignment_identity':ident,'query_coverage_global':qcov,'target_coverage_global':tcov,'global_alignment_length':alen,'exact_query_sequence':tseq==queries[fam],'e3d_official_witness':target in witness_ids[fam],'target_aa_sha256':hashlib.sha256(tseq.encode()).hexdigest()})
    write_tsv(out/'STAGE05E3E_PHMMER_HITS.tsv',hitrows,['family_id','fam50_cluster','target_name','sequence_id','full_evalue','full_score','best_domain_evalue','best_domain_score','target_aa_length','global_alignment_identity','query_coverage_global','target_coverage_global','global_alignment_length','exact_query_sequence','e3d_official_witness','target_aa_sha256'])
    write_tsv(out/'STAGE05E3E_PHMMER_COMMANDS.tsv',commands,['family_id','command','tblout','domtblout','stdout'])

    # Resolve coordinate registry based on the actual phmmer target IDs.
    regidx,regmeta=discover_coordinate_registry(stage_dir,{r['target_name'] for r in hitrows})
    if regidx is None:
        raise RuntimeError(f'could not resolve Stage05C1 predicted ORF coordinate registry under {stage_dir}')

    # Repeat colocalization around each homolog.
    genome_cache={}; genome_paths={}; coloc=[]
    for h in hitrows:
        rid=h['target_name']; fam=h['family_id']; rec=regidx.get(rid)
        if rec is None:
            coloc.append({**{k:h[k] for k in ('family_id','target_name','sequence_id','full_evalue','e3d_official_witness')},'coordinate_resolved':False,'genome_resolved':False,'nearest_repeat_found':False,'exact_repeat_within_500bp':False,'exact_repeat_within_2000bp':False})
            continue
        acc=rec['sequence_id'] or h['sequence_id']; gp=genome_paths.get(acc)
        if acc not in genome_paths:
            gp=resolve_genome_fasta(project,acc); genome_paths[acc]=gp
        if gp is None:
            coloc.append({**{k:h[k] for k in ('family_id','target_name','sequence_id','full_evalue','e3d_official_witness')},'coordinate_resolved':True,'genome_resolved':False,'orf_start0':rec['start0'],'orf_end0':rec['end0'],'orf_strand':rec['strand'],'nearest_repeat_found':False,'exact_repeat_within_500bp':False,'exact_repeat_within_2000bp':False})
            continue
        if acc not in genome_cache: genome_cache[acc]=read_single_genome(gp)
        g=genome_cache[acc]; rad=int(policy['exact_repeat_search_radius_bp']); lo=max(0,rec['start0']-rad); hi=min(len(g),rec['end0']+rad)
        local=g[lo:hi]; mh=scan_motif(local,motifs[fam]); events=[]
        for rs,re_,ori in mh:
            rs+=lo; re_+=lo; side,gap=topology_and_gap(rec['start0'],rec['end0'],rec['strand'],rs,re_); events.append((gap,rs,re_,ori,side))
        events.sort(key=lambda x:(x[0],x[1],x[2],x[3]))
        if events:
            gap,rs,re_,ori,side=events[0]
            found=True
        else:
            gap=rs=re_=''; ori=side=''; found=False
        coloc.append({'family_id':fam,'target_name':rid,'sequence_id':acc,'full_evalue':h['full_evalue'],'e3d_official_witness':h['e3d_official_witness'],'coordinate_resolved':True,'genome_resolved':True,'orf_start0':rec['start0'],'orf_end0':rec['end0'],'orf_strand':rec['strand'],'repeat_core':motifs[fam],'repeat_core_length':len(motifs[fam]),'repeat_match_n_within_2kb_window':len(events),'nearest_repeat_found':found,'nearest_repeat_start0':rs,'nearest_repeat_end0':re_,'nearest_repeat_genomic_orientation':ori,'nearest_repeat_relative_side':side,'nearest_repeat_gap_bp':gap,'exact_repeat_within_500bp':bool(found and gap<=int(policy['exact_repeat_close_radius_bp'])),'exact_repeat_within_2000bp':bool(found),'genome_fasta':str(gp),'genome_sha256':sha256_file(gp)})
    write_tsv(out/'STAGE05E3E_REPEAT_COLOCALIZATION.tsv',coloc,['family_id','target_name','sequence_id','full_evalue','e3d_official_witness','coordinate_resolved','genome_resolved','orf_start0','orf_end0','orf_strand','repeat_core','repeat_core_length','repeat_match_n_within_2kb_window','nearest_repeat_found','nearest_repeat_start0','nearest_repeat_end0','nearest_repeat_genomic_orientation','nearest_repeat_relative_side','nearest_repeat_gap_bp','exact_repeat_within_500bp','exact_repeat_within_2000bp','genome_fasta','genome_sha256'])

    # RT context is optional and only uses pre-existing Stage05C1 RT loci.
    rtidx,rtmeta=resolve_rt_loci(project); rtrows=[]
    if rtidx:
        for c in coloc:
            if not truthy(c.get('coordinate_resolved')): continue
            acc=c['sequence_id']; arr=rtidx.get(acc,[]); best=None
            for s,e,r in arr:
                d=interval_distance(int(c['orf_start0']),int(c['orf_end0']),s,e)
                if best is None or d<best[0]: best=(d,s,e,r)
            rtrows.append({'family_id':c['family_id'],'target_name':c['target_name'],'sequence_id':acc,'nearest_rt_found':best is not None,'nearest_rt_distance_bp':best[0] if best else '','nearest_rt_start0':best[1] if best else '','nearest_rt_end0':best[2] if best else '','within_10kb_of_stage05c1_rt':bool(best and best[0]<=int(policy['rt_context_radius_bp']))})
    write_tsv(out/'STAGE05E3E_RT_CONTEXT.tsv',rtrows,['family_id','target_name','sequence_id','nearest_rt_found','nearest_rt_distance_bp','nearest_rt_start0','nearest_rt_end0','within_10kb_of_stage05c1_rt'])

    # Family summaries.
    fsum=[]; colsum=[]; rtsum=[]
    for fam in priority:
        hh=[x for x in hitrows if x['family_id']==fam]; cc=[x for x in coloc if x['family_id']==fam]
        strong=[x for x in hh if float(x['full_evalue'])<=float(policy['phmmer_strong_evalue'])]
        exact=[x for x in hh if truthy(x['exact_query_sequence'])]
        close=[x for x in cc if truthy(x.get('exact_repeat_within_500bp'))]
        wcc=[x for x in cc if truthy(x.get('e3d_official_witness'))]
        witness_close=[x for x in wcc if truthy(x.get('exact_repeat_within_500bp'))]
        sides=[x.get('nearest_repeat_relative_side') for x in close if x.get('nearest_repeat_relative_side')]
        side_mode=max(set(sides),key=sides.count) if sides else ''
        gaps=[int(x['nearest_repeat_gap_bp']) for x in close if str(x.get('nearest_repeat_gap_bp','')).isdigit()]
        fsum.append({'family_id':fam,'fam50_cluster':pmap[fam]['fam50_cluster'],'e3d_witness_n':len(witness_ids[fam]),'representative_aa_length':len(queries[fam]),'phmmer_report_hits_n':len(hh),'phmmer_report_genomes_n':len({x['sequence_id'] for x in hh}),'phmmer_strong_hits_n':len(strong),'phmmer_strong_genomes_n':len({x['sequence_id'] for x in strong}),'exact_sequence_hits_n':len(exact),'exact_sequence_genomes_n':len({x['sequence_id'] for x in exact}),'homologs_with_exact_repeat_within_500bp_n':len(close),'homolog_genomes_with_exact_repeat_within_500bp_n':len({x['sequence_id'] for x in close}),'e3d_witnesses_with_exact_repeat_within_500bp_n':len(witness_close),'e3d_witness_repeat_close_fraction':(len(witness_close)/len(wcc) if wcc else ''),'colocalized_repeat_modal_side':side_mode,'colocalized_repeat_median_gap_bp':median(gaps),'frozen_e3b_modal_side':pmap[fam]['frozen_e3b_modal_side'],'frozen_e3b_median_gap_bp':pmap[fam]['frozen_e3b_median_gap_bp']})
        colsum.append({'family_id':fam,'phmmer_homolog_n':len(hh),'coordinate_resolved_n':sum(truthy(x.get('coordinate_resolved')) for x in cc),'genome_resolved_n':sum(truthy(x.get('genome_resolved')) for x in cc),'exact_repeat_within_500bp_n':len(close),'exact_repeat_within_2000bp_n':sum(truthy(x.get('exact_repeat_within_2000bp')) for x in cc),'modal_side_within_500bp':side_mode,'median_gap_within_500bp':median(gaps)})
        rr=[x for x in rtrows if x['family_id']==fam]
        rtsum.append({'family_id':fam,'rt_context_available':bool(rtidx),'homologs_with_rt_context_n':len(rr),'within_10kb_of_stage05c1_rt_n':sum(truthy(x.get('within_10kb_of_stage05c1_rt')) for x in rr),'within_10kb_of_stage05c1_rt_genomes_n':len({x['sequence_id'] for x in rr if truthy(x.get('within_10kb_of_stage05c1_rt'))})})
    write_tsv(out/'STAGE05E3E_FAMILY_CLOSURE_SUMMARY.tsv',fsum,list(fsum[0]))
    write_tsv(out/'STAGE05E3E_EXACT_REPEAT_COLOCALIZATION_SUMMARY.tsv',colsum,list(colsum[0]))
    write_tsv(out/'STAGE05E3E_RT_CONTEXT_SUMMARY.tsv',rtsum,list(rtsum[0]))

    # Core integrity gate: every priority family must return its E3D witnesses in phmmer and resolve/scan their genomes.
    core_checks={}
    for fam in priority:
        hids={x['target_name'] for x in hitrows if x['family_id']==fam}
        cc=[x for x in coloc if x['family_id']==fam and truthy(x.get('e3d_official_witness'))]
        core_checks[f'{fam}_all_e3d_witnesses_recovered_by_phmmer']=witness_ids[fam].issubset(hids)
        core_checks[f'{fam}_all_e3d_witness_coordinates_resolved']=len(cc)==len(witness_ids[fam]) and all(truthy(x.get('coordinate_resolved')) for x in cc)
        core_checks[f'{fam}_all_e3d_witness_genomes_resolved']=len(cc)==len(witness_ids[fam]) and all(truthy(x.get('genome_resolved')) for x in cc)
    gate_pass=all(core_checks.values())
    resource={'phmmer_path':shutil.which('phmmer'),'phmmer_version':subprocess.run([shutil.which('phmmer'),'-h'],capture_output=True,text=True).stdout.splitlines()[1:3] if shutil.which('phmmer') else [],'full919_proteome_path':str(proteome),'full919_proteome_sha256':proteome_sha,'coordinate_registry':regmeta,'rt_loci':rtmeta,'resolved_genome_fasta_n':sum(1 for x in genome_paths.values() if x is not None),'requested_genome_fasta_n':len(genome_paths)}
    (out/'STAGE05E3E_RESOURCE_STATUS.json').write_text(json.dumps(resource,indent=2)+"\n")
    gateobj={'stage':STAGE,'version':VERSION,'gate':'PASS_E3E_PRIORITY_FAMILY_LOCAL_HOMOLOGY_AND_REPEAT_COLOCALIZATION_CLOSURE' if gate_pass else 'HOLD_E3E_PRIORITY_FAMILY_LOCAL_CLOSURE','gate_pass':gate_pass,'priority_families':priority,'core_checks':core_checks,'family_membership_changed':False,'candidate_filtering_performed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,'retuning_authorized':False,'external_global_novelty_audit_performed':False,'structure_search_performed':False,'claim_boundary':policy['claim_boundary']}
    (out/'STAGE05E3E_GATE.json').write_text(json.dumps(gateobj,indent=2)+"\n")
    summary={'stage':STAGE,'version':VERSION,'gate':gateobj['gate'],'gate_pass':gate_pass,'priority_family_n':len(priority),'phmmer_total_hits_n':len(hitrows),'phmmer_unique_target_n':len({x['target_name'] for x in hitrows}),'families':{r['family_id']:r for r in fsum},'rt_context_status':rtmeta,'next_required_after_e3e':'Interpret local expansion; then perform separately frozen external/public-database remote-homology/structure/novelty work before any novelty or function claim.'}
    (out/'STAGE05E3E_SUMMARY.json').write_text(json.dumps(summary,indent=2)+"\n")

    notes=['# Stage05E3E manuscript-facing notes','', 'Stage05E3E is a local 919-genome biological-closure layer. It does not change the E2A/E3D family membership and does not establish global novelty or function.','']
    for r in fsum:
        notes.append(f"## {r['family_id']}")
        notes.append(f"- Frozen family: {r['fam50_cluster']}; E3D official-supported witnesses: {r['e3d_witness_n']}; representative length: {r['representative_aa_length']} aa.")
        notes.append(f"- Local phmmer expansion: {r['phmmer_report_hits_n']} proteins in {r['phmmer_report_genomes_n']} genomes at E<={policy['phmmer_report_evalue']}; {r['phmmer_strong_hits_n']} proteins in {r['phmmer_strong_genomes_n']} genomes at E<={policy['phmmer_strong_evalue']}.")
        notes.append(f"- Exact frozen repeat core within 500 bp of a local homolog: {r['homologs_with_exact_repeat_within_500bp_n']} homologs in {r['homolog_genomes_with_exact_repeat_within_500bp_n']} genomes; modal side={r['colocalized_repeat_modal_side'] or 'NA'}, median gap={r['colocalized_repeat_median_gap_bp'] or 'NA'} bp.")
        notes.append(f"- E3D witness exact-repeat close recovery: {r['e3d_witnesses_with_exact_repeat_within_500bp_n']}/{r['e3d_witness_n']}. Frozen E3B reference topology was {r['frozen_e3b_modal_side']} with median gap {r['frozen_e3b_median_gap_bp']} bp.")
        notes.append('')
    notes += ['## Claim boundary','Positive local homology and repeat colocalization support a replicated architecture; they do not by themselves prove RNA expression, translation, interaction, reverse-transcription activity, ART identity, or natural-world uniqueness. Negative phmmer/RT context is not proof of novelty.','', '## PB50 manuscript lock','The separate PB50 manuscript lock in `docs/PB50_MANUSCRIPT_LOCK.md` remains mandatory for final paper integration.']
    (out/'STAGE05E3E_MANUSCRIPT_NOTES.md').write_text('\n'.join(notes)+'\n')
    print(json.dumps(gateobj,indent=2))
    if not gate_pass: raise SystemExit(2)

if __name__=='__main__': main()
