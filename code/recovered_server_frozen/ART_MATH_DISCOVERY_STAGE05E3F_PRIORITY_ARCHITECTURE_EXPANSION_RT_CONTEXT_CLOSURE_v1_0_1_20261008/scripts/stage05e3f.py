#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, hashlib, io, json, os, re, shutil, statistics, tarfile

FAMS=['ARCH250_0010','ARCH250_0013','ARCH250_0033']

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()

def read_tsv(p):
    with open(p,newline='') as f: return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(p,rows,cols=None):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    if cols is None: cols=list(rows[0]) if rows else []
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n',extrasaction='ignore')
        if cols: w.writeheader(); w.writerows(rows)

def truthy(x): return str(x).lower() in ('1','true','yes','y')

def materialize(src, project, dest):
    src=Path(src)
    if not src.exists():
        hits=list(Path(project).rglob(src.name))
        hits=[x for x in hits if x.is_file()]
        if len(hits)!=1: raise RuntimeError(f'cannot uniquely resolve input {src.name}; hits={len(hits)}')
        print(f'[RESOLVE] {src} -> {hits[0]}')
        src=hits[0]
    if dest.exists(): shutil.rmtree(dest)
    dest.mkdir(parents=True)
    if src.is_dir():
        for x in src.rglob('*'):
            if x.is_file():
                q=dest/x.relative_to(src); q.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(x,q)
    else:
        with tarfile.open(src,'r:*') as t:
            for m in t.getmembers():
                if not m.isfile(): continue
                pp=Path(m.name)
                if pp.is_absolute() or '..' in pp.parts: raise RuntimeError('unsafe archive member')
            t.extractall(dest)
    return src

def one_by_name(root,name):
    hits=list(Path(root).rglob(name))
    if len(hits)!=1: raise RuntimeError(f'expected one {name}; found {len(hits)}')
    return hits[0]

def parse_fasta(p):
    d={}; name=None; seq=[]
    for line in open(p):
        line=line.strip()
        if not line: continue
        if line.startswith('>'):
            if name is not None: d[name]=''.join(seq)
            name=line[1:].split()[0]; seq=[]
        else: seq.append(line)
    if name is not None: d[name]=''.join(seq)
    return d

def find_project_file(project,name, preferred_fragment=None):
    hits=[]
    for x in Path(project).rglob(name):
        if not x.is_file(): continue
        s=str(x)
        if 'STAGE05E3F_PRIORITY_ARCHITECTURE' in s: continue
        hits.append(x)
    if preferred_fragment:
        pref=[x for x in hits if preferred_fragment in str(x)]
        if len(pref)==1: return pref[0], {'route':'direct_preferred','path':str(pref[0]),'sha256':sha256_file(pref[0])}
    if len(hits)==1: return hits[0], {'route':'direct_unique','path':str(hits[0]),'sha256':sha256_file(hits[0])}
    if len(hits)>1:
        hs={sha256_file(x) for x in hits}
        if len(hs)==1:
            p=sorted(hits,key=lambda x:(0 if '/outputs/' in str(x) else 1,len(str(x)),str(x)))[0]
            return p, {'route':'direct_identical_duplicate','path':str(p),'sha256':sha256_file(p),'duplicate_n':len(hits)}
    return None, {'route':'direct_not_resolved','hits':[str(x) for x in hits[:20]]}

def extract_member_from_stage05c1_tar(project, work):
    candidates=[]
    exact='ART_STAGE05C1_NEW919_BIOLOGY_RESULTS_20261006_092047.tar.gz'
    for x in Path(project).rglob('*.tar.gz'):
        if not x.is_file(): continue
        n=x.name
        if n==exact or ('STAGE05C1' in n and 'RESULT' in n): candidates.append(x)
    found=[]
    for tar in sorted(candidates):
        try:
            with tarfile.open(tar,'r:*') as t:
                ms=[m for m in t.getmembers() if m.isfile() and Path(m.name).name=='RT_LOCI_ALL919.tsv']
                for m in ms:
                    b=t.extractfile(m).read(); found.append((tar,m.name,b,hashlib.sha256(b).hexdigest()))
        except tarfile.TarError:
            pass
    if not found: return None, {'route':'stage05c1_tar_not_resolved','candidate_tar_n':len(candidates)}
    hashes={x[3] for x in found}
    if len(hashes)!=1:
        return None, {'route':'stage05c1_tar_ambiguous_nonidentical','found':[(str(a),b,d) for a,b,c,d in found]}
    found.sort(key=lambda x:(0 if x[0].name==exact else 1,str(x[0]),x[1]))
    tar,member,b,h=found[0]
    out=work/'RT_LOCI_ALL919.from_tar.tsv'; out.write_bytes(b)
    return out, {'route':'stage05c1_result_tar','tar':str(tar),'member':member,'sha256':h,'matching_copies_n':len(found)}

def resolve_rt_loci(project, work):
    p,meta=find_project_file(project,'RT_LOCI_ALL919.tsv','ART_MATH_DISCOVERY_STAGE05C1_NEW919_RT_ORF_REPEAT_ARCHITECTURE_v1_0_0_20261006')
    if p: return p,meta
    return extract_member_from_stage05c1_tar(project,work)

def load_rt_accession_lock():
    p=Path(__file__).resolve().parents[1]/'config'/'STAGE05C1_RT_CORE_ACCESSION_LOCK.tsv'
    rows=read_tsv(p)
    acc=[r.get('accession','') for r in rows if r.get('accession','')]
    if len(acc)!=10 or len(set(acc))!=10:
        raise RuntimeError(f'RT accession lock must contain exactly 10 unique accessions; observed={len(acc)} unique={len(set(acc))}')
    return p,rows,set(acc)

def choose_col(rows,aliases):
    if not rows: return None
    cols=rows[0].keys()
    for x in aliases:
        if x in cols: return x
    return None

def parse_rt_index(rt_path, gene_registry):
    rows=read_tsv(rt_path); greg=read_tsv(gene_registry)
    gm={r.get('record_key',''):r for r in greg if r.get('record_key')}
    ac=choose_col(rows,['sequence_id','accession','acc','genome_id'])
    sc=choose_col(rows,['orf_start0','start0','cds_start0','rt_start0','start'])
    ec=choose_col(rows,['orf_end0','end0','cds_end0','rt_end0','end'])
    rc=choose_col(rows,['record_key','protein_id','target_name','orf_id'])
    idx={}; parsed=0; via_registry=0
    for r in rows:
        acc=r.get(ac,'') if ac else ''
        s=r.get(sc,'') if sc else ''; e=r.get(ec,'') if ec else ''
        if (not acc or s=='' or e=='') and rc and r.get(rc,'') in gm:
            z=gm[r[rc]]; acc=z.get('sequence_id',''); s=z.get('start0',''); e=z.get('end0',''); via_registry+=1
        try: s=int(float(s)); e=int(float(e))
        except: continue
        if not acc and rc and r.get(rc): acc=r[rc].split('__PRODIGAL_')[0]
        if not acc: continue
        idx.setdefault(acc,[]).append((s,e,r)); parsed+=1
    meta={'rows_n':len(rows),'parsed_loci_n':parsed,'via_gene_registry_n':via_registry,'columns':list(rows[0]) if rows else [],'accession_col':ac,'start_col':sc,'end_col':ec,'record_col':rc}
    if parsed==0: raise RuntimeError('RT_LOCI_ALL919.tsv resolved but zero loci could be parsed')
    return idx,meta

def interval_distance(a0,a1,b0,b1):
    if a1<=b0: return b0-a1
    if b1<=a0: return a0-b1
    return 0

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--e3e',required=True); ap.add_argument('--project-root',required=True); ap.add_argument('--outdir',default='outputs'); ap.add_argument('--workdir',default='work')
    a=ap.parse_args(); project=Path(a.project_root); out=Path(a.outdir); work=Path(a.workdir)
    for p in [out,work]:
        if p.exists(): shutil.rmtree(p)
        p.mkdir(parents=True)
    e3e_src=materialize(a.e3e,project,work/'e3e')
    gate=json.load(open(one_by_name(work/'e3e','STAGE05E3E_GATE.json')))
    if not gate.get('gate_pass') or gate.get('gate')!='PASS_E3E_PRIORITY_FAMILY_LOCAL_HOMOLOGY_AND_REPEAT_COLOCALIZATION_CLOSURE': raise RuntimeError('E3E gate not accepted')
    hits=read_tsv(one_by_name(work/'e3e','STAGE05E3E_PHMMER_HITS.tsv'))
    coloc=read_tsv(one_by_name(work/'e3e','STAGE05E3E_REPEAT_COLOCALIZATION.tsv'))
    lock=read_tsv(one_by_name(work/'e3e','STAGE05E3E_PRIORITY_FAMILY_LOCK.tsv'))
    rstat=json.load(open(one_by_name(work/'e3e','STAGE05E3E_RESOURCE_STATUS.json')))
    if set(x['family_id'] for x in lock)!=set(FAMS): raise RuntimeError('priority family lock mismatch')
    lmap={x['family_id']:x for x in lock}; cmap={(x['family_id'],x['target_name']):x for x in coloc}
    proteome=Path(rstat['full919_proteome_path']); greg=Path(rstat['coordinate_registry']['path'])
    if not proteome.exists() or not greg.exists(): raise FileNotFoundError('Stage05C1 proteome/gene registry not present at frozen paths')
    if sha256_file(proteome)!=rstat['full919_proteome_sha256']: raise RuntimeError('full919 proteome SHA drift')
    if sha256_file(greg)!=rstat['coordinate_registry']['sha256']: raise RuntimeError('gene registry SHA drift')
    seqs=parse_fasta(proteome)
    missing=[h['target_name'] for h in hits if h['target_name'] not in seqs]
    if missing: raise RuntimeError(f'missing {len(missing)} E3E hit sequences from frozen proteome')
    # Export family hit FASTAs and descriptive expansion rows.
    exp=[]; fsum=[]
    for fam in FAMS:
        hh=[h for h in hits if h['family_id']==fam]
        fasta=out/f'{fam}.all_local_homologs.faa'
        with fasta.open('w') as f:
            for h in hh:
                rid=h['target_name']; f.write(f'>{rid} sequence_id={h["sequence_id"]} evalue={h["full_evalue"]} witness={h["e3d_official_witness"]}\n')
                s=seqs[rid]
                for i in range(0,len(s),80): f.write(s[i:i+80]+'\n')
        non=[h for h in hh if not truthy(h['e3d_official_witness'])]
        nonstrong=[h for h in non if float(h['full_evalue'])<=1e-5]
        non500=[]; non2k=[]; same500=[]
        for h in non:
            c=cmap[(fam,h['target_name'])]
            r={'family_id':fam,'target_name':h['target_name'],'sequence_id':h['sequence_id'],'full_evalue':h['full_evalue'],'target_aa_length':h['target_aa_length'],'global_alignment_identity':h['global_alignment_identity'],'strong_e1e5':float(h['full_evalue'])<=1e-5,'exact_repeat_within_500bp':truthy(c['exact_repeat_within_500bp']),'exact_repeat_within_2000bp':truthy(c['exact_repeat_within_2000bp']),'repeat_side':c.get('nearest_repeat_relative_side',''),'repeat_gap_bp':c.get('nearest_repeat_gap_bp',''),'same_as_frozen_modal_side':c.get('nearest_repeat_relative_side','')==lmap[fam]['frozen_e3b_modal_side'] if c.get('nearest_repeat_relative_side','') else False}
            exp.append(r)
            if r['exact_repeat_within_500bp']: non500.append(r)
            if r['exact_repeat_within_2000bp']: non2k.append(r)
            if r['exact_repeat_within_500bp'] and r['same_as_frozen_modal_side']: same500.append(r)
        if same500:
            cls='EXPANDED_SAME_ARCHITECTURE'
        elif non2k:
            cls='WEAK_2KB_EXTENSION_ONLY'
        elif nonstrong:
            cls='PROTEIN_FAMILY_BROAD_REPEAT_RESTRICTED'
        else:
            cls='WITNESS_CORE_ONLY'
        fsum.append({'family_id':fam,'fam50_cluster':lmap[fam]['fam50_cluster'],'e3d_witness_n':lmap[fam]['witness_n'],'local_homolog_n':len(hh),'nonwitness_local_homolog_n':len(non),'nonwitness_strong_homolog_n':len(nonstrong),'nonwitness_repeat_within_500bp_n':len(non500),'nonwitness_repeat_within_2000bp_n':len(non2k),'nonwitness_same_side_repeat_within_500bp_n':len(same500),'frozen_modal_side':lmap[fam]['frozen_e3b_modal_side'],'frozen_median_gap_bp':lmap[fam]['frozen_e3b_median_gap_bp'],'e3f_architecture_expansion_class':cls})
    write_tsv(out/'STAGE05E3F_NONWITNESS_ARCHITECTURE_EXPANSION.tsv',exp)
    # RT context resolution. Prefer the exact Stage05C1 coordinate table. If that
    # compact resource is absent, a frozen ten-accession lock can still prove
    # zero RT proximity when NONE of the E3E homolog genomes belongs to the
    # Stage05C1 RT-core genome set. This is accession-level logic, not a rescore.
    rt_path,rt_source=resolve_rt_loci(project,work)
    hit_accessions={h['sequence_id'] for h in hits}
    lock_path,lock_rows,rt_accession_lock=load_rt_accession_lock()
    lock_intersection=sorted(hit_accessions & rt_accession_lock)
    rtrows=[]
    if rt_path is not None:
        rtidx,rtparse=parse_rt_index(rt_path,greg); rt_source['parse']=rtparse
        rt_source['resolution_level']='EXACT_COORDINATE_TABLE'
        rt_source['frozen_rt_accession_lock_sha256']=sha256_file(lock_path)
        rt_source['e3e_homolog_accession_intersection_with_rt_lock']=lock_intersection
        for h in hits:
            c=cmap[(h['family_id'],h['target_name'])]
            acc=h['sequence_id']; a0=int(c['orf_start0']); a1=int(c['orf_end0']); arr=rtidx.get(acc,[]); best=None
            for s,e,r in arr:
                d=interval_distance(a0,a1,s,e)
                if best is None or d<best[0]: best=(d,s,e,r)
            rtrows.append({'family_id':h['family_id'],'target_name':h['target_name'],'sequence_id':acc,'e3d_official_witness':h['e3d_official_witness'],'rt_context_resolution_level':'EXACT_COORDINATE_TABLE','same_genome_stage05c1_rt':best is not None,'nearest_rt_found':best is not None,'nearest_rt_distance_bp':best[0] if best else '','nearest_rt_start0':best[1] if best else '','nearest_rt_end0':best[2] if best else '','within_10kb_of_stage05c1_rt':bool(best and best[0]<=10000)})
    else:
        if lock_intersection:
            raise RuntimeError('exact RT_LOCI_ALL919.tsv is unavailable and E3E homolog genomes intersect the frozen 10-genome RT-core accession lock; exact coordinates are required for those accessions: '+','.join(lock_intersection))
        rt_source={'route':'frozen_stage05c1_rt_accession_lock_zero_intersection',
                   'resolution_level':'ACCESSION_LEVEL_ZERO_INTERSECTION',
                   'exact_rt_coordinate_table_resolved':False,
                   'frozen_rt_accession_lock_path':str(lock_path),
                   'frozen_rt_accession_lock_sha256':sha256_file(lock_path),
                   'frozen_rt_accession_n':len(rt_accession_lock),
                   'e3e_homolog_accession_n':len(hit_accessions),
                   'e3e_homolog_accession_intersection_with_rt_lock':lock_intersection,
                   'logic':'No E3E homolog genome is one of the ten frozen Stage05C1 RT-core genomes; therefore no homolog can be within 10 kb of a Stage05C1 RT locus. Coordinates are unnecessary for this zero-intersection conclusion.'}
        for h in hits:
            rtrows.append({'family_id':h['family_id'],'target_name':h['target_name'],'sequence_id':h['sequence_id'],'e3d_official_witness':h['e3d_official_witness'],'rt_context_resolution_level':'ACCESSION_LEVEL_ZERO_INTERSECTION','same_genome_stage05c1_rt':False,'nearest_rt_found':False,'nearest_rt_distance_bp':'','nearest_rt_start0':'','nearest_rt_end0':'','within_10kb_of_stage05c1_rt':False})
    write_tsv(out/'STAGE05E3F_RT_CONTEXT.tsv',rtrows)
    rtsum=[]
    for fam in FAMS:
        rr=[r for r in rtrows if r['family_id']==fam]; near=[r for r in rr if truthy(r['within_10kb_of_stage05c1_rt'])]
        levels=sorted({r['rt_context_resolution_level'] for r in rr})
        rtsum.append({'family_id':fam,'rt_context_available':True,'rt_context_resolution_level':';'.join(levels),'homologs_audited_n':len(rr),'homologs_with_any_same_genome_stage05c1_rt_n':sum(truthy(r['same_genome_stage05c1_rt']) for r in rr),'within_10kb_of_stage05c1_rt_n':len(near),'within_10kb_of_stage05c1_rt_genomes_n':len({r['sequence_id'] for r in near})})
    write_tsv(out/'STAGE05E3F_RT_CONTEXT_SUMMARY.tsv',rtsum)
    # Add RT summary to family summary.
    rtmap={r['family_id']:r for r in rtsum}
    for r in fsum:
        z=rtmap[r['family_id']]; r.update({'within_10kb_of_stage05c1_rt_n':z['within_10kb_of_stage05c1_rt_n'],'within_10kb_of_stage05c1_rt_genomes_n':z['within_10kb_of_stage05c1_rt_genomes_n']})
    write_tsv(out/'STAGE05E3F_FAMILY_EXPANSION_SUMMARY.tsv',fsum)
    # Manuscript notes and gate.
    notes=['# Stage05E3F manuscript-facing notes','', 'E3F separates protein-family expansion from exact repeat-linked architecture expansion and repairs the previously unavailable Stage05C1 RT-context join. When the exact RT coordinate table is absent, accession-level zero-intersection against the frozen ten-genome RT-core lock is allowed only if none of the E3E homolog genomes belongs to that lock. It does not alter frozen family membership or detector scores.','']
    for r in fsum:
        notes += [f'## {r["family_id"]}',f'- Descriptive class: **{r["e3f_architecture_expansion_class"]}**.',f'- E3D witnesses: {r["e3d_witness_n"]}; local homologs: {r["local_homolog_n"]}; non-witness strong homologs: {r["nonwitness_strong_homolog_n"]}.',f'- Non-witness exact frozen repeat core within 500 bp: {r["nonwitness_repeat_within_500bp_n"]}; within 2 kb: {r["nonwitness_repeat_within_2000bp_n"]}; same frozen modal side within 500 bp: {r["nonwitness_same_side_repeat_within_500bp_n"]}.',f'- Stage05C1 RT-core loci within 10 kb of E3E local homologs: {r["within_10kb_of_stage05c1_rt_n"]} in {r["within_10kb_of_stage05c1_rt_genomes_n"]} genomes.','']
    notes += ['## Claim boundary','These are within-project replication and context results. They do not prove global novelty, function, RNA expression, translation, interaction, ART identity, or biochemical activity. Negative RT proximity is not proof that a protein is unrelated to reverse-transcription biology.','', '## PB50 lock','The mandatory PB50 manuscript lock remains in `docs/PB50_MANUSCRIPT_LOCK.md` and must be integrated unchanged in scientific meaning.']
    (out/'STAGE05E3F_MANUSCRIPT_NOTES.md').write_text('\n'.join(notes)+'\n')
    (out/'STAGE05E3F_RT_RESOURCE_STATUS.json').write_text(json.dumps(rt_source,indent=2)+'\n')
    gate={'stage':'Stage05E3F','version':'v1.0.1','gate':'PASS_E3F_PRIORITY_ARCHITECTURE_EXPANSION_AND_RT_CONTEXT_CLOSURE','gate_pass':True,'priority_families':FAMS,'e3e_input_sha256':sha256_file(e3e_src) if Path(e3e_src).is_file() else None,'rt_context_resolved':True,'rt_context_resolution_level':rt_source.get('resolution_level','UNKNOWN'),'rt_coordinate_table_resolved':rt_path is not None,'frozen_rt_accession_lock_intersection_n':len(lock_intersection),'rt_loci_source':rt_source,'all_e3e_hit_sequences_recovered':True,'family_membership_changed':False,'candidate_filtering_performed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,'retuning_authorized':False,'external_global_novelty_audit_performed':False,'structure_search_performed':False,'claim_boundary':'Within-project architecture expansion and pre-existing RT-context closure only; no global novelty/function/ART/wet-lab claim.'}
    (out/'STAGE05E3F_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate,indent=2))
    print('[PASS] Stage05E3F architecture expansion + RT context closure')

if __name__=='__main__': main()
