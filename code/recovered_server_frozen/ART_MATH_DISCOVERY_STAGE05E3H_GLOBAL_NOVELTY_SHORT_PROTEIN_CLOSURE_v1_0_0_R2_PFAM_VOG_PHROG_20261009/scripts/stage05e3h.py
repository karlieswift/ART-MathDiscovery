#!/usr/bin/env python3
import argparse, csv, hashlib, json, math, os, re, shutil, subprocess, tarfile, tempfile
from pathlib import Path
from collections import Counter, defaultdict

FAMS = ["ARCH250_0010","ARCH250_0013","ARCH250_0033"]
KD = {"I":4.5,"V":4.2,"L":3.8,"F":2.8,"C":2.5,"M":1.9,"A":1.8,"G":-0.4,"T":-0.7,"S":-0.8,"W":-0.9,"Y":-1.3,"P":-1.6,"H":-3.2,"E":-3.5,"Q":-3.5,"D":-3.5,"N":-3.5,"K":-3.9,"R":-4.5}
HYDRO=set("AVILMFWYC")
HELIX_PROXY=set("ALEKQRHM")
AA20=set("ACDEFGHIKLMNPQRSTVWY")


def sha256_text(s): return hashlib.sha256(s.encode()).hexdigest()
def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20), b''): h.update(b)
    return h.hexdigest()

def resolve(path, root):
    p=Path(path)
    if p.exists(): return p.resolve()
    name=p.name
    hits=[]
    for q in Path(root).rglob(name):
        if q.is_file(): hits.append(q.resolve())
    if len(hits)==1:
        print(f"[RESOLVE] {path} -> {hits[0]}")
        return hits[0]
    raise RuntimeError(f"could not uniquely resolve {path}; hits={len(hits)}")

def extract_tar(tar_path, dest):
    dest=Path(dest); dest.mkdir(parents=True,exist_ok=True)
    with tarfile.open(tar_path,'r:gz') as tf: tf.extractall(dest)
    return dest

def read_tsv(path):
    with open(path,newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(path, rows, fields):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(rows)

def fasta(path):
    recs=[]; h=None; seq=[]
    with open(path) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if h is not None: recs.append((h,''.join(seq).upper()))
                h=line[1:]; seq=[]
            else: seq.append(line)
        if h is not None: recs.append((h,''.join(seq).upper()))
    return recs

def entropy(seq):
    c=Counter(seq); n=len(seq)
    if not n: return 0.0
    return -sum((v/n)*math.log2(v/n) for v in c.values())

def max_window(seq, w, fn):
    if len(seq)<=w: return fn(seq) if seq else 0.0
    return max(fn(seq[i:i+w]) for i in range(len(seq)-w+1))

def longest_run(seq, charset):
    best=cur=0
    for a in seq:
        if a in charset: cur+=1; best=max(best,cur)
        else: cur=0
    return best

def seq_features(fam, header, seq, witness):
    n=len(seq); cnt=Counter(seq)
    unknown=sorted(set(seq)-AA20)
    if unknown: raise RuntimeError(f"noncanonical AA {unknown} in {fam} {header}")
    kd_mean=sum(KD[a] for a in seq)/n
    kd19=max_window(seq,19,lambda s: sum(KD[a] for a in s)/len(s))
    hyd19=max_window(seq,19,lambda s: sum(a in HYDRO for a in s)/len(s))
    helix19=max_window(seq,19,lambda s: sum(a in HELIX_PROXY for a in s)/len(s))
    nterm=seq[:min(25,n)]
    nterm_kd=sum(KD[a] for a in nterm)/len(nterm)
    charge=(cnt['K']+cnt['R']+0.1*cnt['H'])-(cnt['D']+cnt['E'])
    tmh = bool(n>=19 and kd19>=1.6 and hyd19>=0.68)
    return {
      'family_id':fam,'target_name':header.split()[0],'witness':str(bool(witness)),'aa_length':n,'aa_sha256':sha256_text(seq),
      'mean_kd_hydropathy':round(kd_mean,6),'max_19aa_kd_hydropathy':round(kd19,6),'max_19aa_hydrophobic_fraction':round(hyd19,6),
      'nterm25_mean_kd_hydropathy':round(nterm_kd,6),'longest_hydrophobic_run':longest_run(seq,HYDRO),
      'net_charge_proxy_pH7':round(charge,3),'acidic_fraction':round((cnt['D']+cnt['E'])/n,6),'basic_fraction':round((cnt['K']+cnt['R']+cnt['H'])/n,6),
      'alanine_fraction':round(cnt['A']/n,6),'glutamate_fraction':round(cnt['E']/n,6),'proline_fraction':round(cnt['P']/n,6),'glycine_fraction':round(cnt['G']/n,6),
      'max_19aa_helix_propensity_proxy':round(helix19,6),'sequence_entropy_bits':round(entropy(seq),6),
      'strong_hydrophobic_segment_heuristic':str(tmh)
    }

def find_one(root, basename):
    hits=[p for p in Path(root).rglob(basename) if p.is_file()]
    return sorted(hits)[0] if hits else None

def evidence_index_from_e3g(e3gdir):
    p=Path(e3gdir)/'outputs/STAGE05E3G_EXISTING_ORTHOGONAL_EVIDENCE.tsv'
    rows=read_tsv(p)
    by=defaultdict(list)
    for r in rows: by[r['family_id']].append(r)
    return by

def classify_existing(lines, layer):
    text=' || '.join(lines)
    if not text: return 'NO_MATCH_IN_RECOVERED_EXISTING_ARTIFACTS'
    if layer=='D4_R2_SWISSPROT':
        if 'NO_SWISSPROT_HIT_NOT_NEGATIVE' in text or re.search(r'\t0\t', '\t'+text+'\t'): return 'NO_SWISSPROT_HIT_NOT_NEGATIVE'
        return 'EXISTING_SWISSPROT_RECORD_REVIEW'
    if layer=='D4_R3_STRUCTURE_TIER':
        if 'TIER_S_SHORT_LT180' in text: return 'TIER_S_SHORT_LT180'
        if 'TIER_A_COMPACT_STRUCTURE_READY' in text: return 'TIER_A_COMPACT_STRUCTURE_READY'
        return 'EXISTING_STRUCTURE_TIER_RECORD_REVIEW'
    if layer=='D4_R4_PROSTT5_FOLDSEEK':
        if 'NO_AFSP_STRUCTURAL_HIT_NOT_NEGATIVE' in text: return 'NO_AFSP_STRUCTURAL_HIT_NOT_NEGATIVE'
        if 'TIER_S_SHORT_LT180' in text: return 'NOT_SEARCHED_OLD_D4_SHORT_LT180'
        return 'EXISTING_R4_RECORD_REVIEW'
    if layer=='D4_R5_FULL3D':
        if 'NO_FULL3D_AFSP_HIT_NOT_NEGATIVE' in text or 'NO_HIT' in text: return 'NO_FULL3D_CURATED_OR_RT_REFERENCE_HIT_NOT_NEGATIVE'
        return 'EXISTING_R5_RECORD_REVIEW'
    return 'MATCHED_EXISTING_ARTIFACT'

def run_cmd(cmd, stdout=None, timeout=3600):
    try:
        p=subprocess.run(cmd,stdout=stdout or subprocess.PIPE,stderr=subprocess.PIPE,text=(stdout is None),timeout=timeout,check=False)
        return p.returncode, (p.stdout if stdout is None else ''), p.stderr
    except Exception as e: return 99,'',str(e)


def recover_all_homolog_existing_evidence(root, feat_rows, e3g_existing):
    specs=[
      ('D4_R2_SWISSPROT','STAGE05D4R2_REMOTE_HOMOLOGY_TRIAGE.tsv'),
      ('D4_R3_STRUCTURE_TIER','STAGE05D4R3_NR90_STRUCTURE_PRIORITY.tsv'),
      ('D4_R4_PROSTT5_FOLDSEEK','STAGE05D4R4_TIER_A_STRUCTURAL_TRIAGE.tsv'),
      ('D4_R5_FULL3D','STAGE05D4R5_P1_FULL3D_TRIAGE.tsv')]
    files={layer:find_one(root,base) for layer,base in specs}
    file_lines={}
    for layer,p in files.items():
        if p:
            try: file_lines[layer]=p.read_text(errors='replace').splitlines()
            except Exception: file_lines[layer]=[]
        else: file_lines[layer]=[]
    out=[]
    for fr in feat_rows:
        fam=fr['family_id']; sha=fr['aa_sha256']; target=fr['target_name']
        for layer,_ in specs:
            lines=[]
            for line in file_lines[layer]:
                if sha in line or ('AA_'+sha) in line: lines.append(line[:4000])
            route='canonical_project_output' if lines else ''
            # Representative fallback to E3G exact evidence if canonical file not visible.
            if not lines:
                es=[r for r in e3g_existing.get(fam,[]) if r['stage_layer']==layer and (sha in r.get('matched_tokens','') or sha in r.get('matched_line_excerpt',''))]
                if es:
                    lines=[r.get('matched_line_excerpt','') for r in es]; route='e3g_recovered_exact_evidence'
            out.append({'family_id':fam,'target_name':target,'aa_sha256':sha,'stage_layer':layer,'semantic_status':classify_existing(lines,layer),'matched_line_n':len(lines),'resolution_route':route,'source_file':str(files[layer] or '')})
    return out

def optional_pfam(query_faa, outdir, resource_rows, threads):
    hmm=os.environ.get('E3H_PFAM_HMM','').strip(); hmmscan=shutil.which('hmmscan')
    if not hmm:
        resource_rows.append({'resource':'Pfam-A','status':'NOT_PROVIDED','path':'','action':'set E3H_PFAM_HMM=/path/to/Pfam-A.hmm'})
        return False
    if not Path(hmm).exists() or not hmmscan:
        resource_rows.append({'resource':'Pfam-A','status':'UNAVAILABLE','path':hmm,'action':'provide valid HMM and hmmscan'})
        return False
    # hmmscan needs pressed DB; press in-place only if writable, otherwise fail closed.
    pressed=all(Path(hmm+s).exists() for s in ['.h3f','.h3i','.h3m','.h3p'])
    if not pressed:
        hp=shutil.which('hmmpress')
        if hp and os.access(str(Path(hmm).parent),os.W_OK):
            rc,_,err=run_cmd([hp,'-f',hmm],timeout=3600)
            pressed=(rc==0)
        if not pressed:
            resource_rows.append({'resource':'Pfam-A','status':'UNPRESSED_OR_READONLY','path':hmm,'action':'hmmpress Pfam-A.hmm or provide pressed database'})
            return False
    dom=Path(outdir)/'STAGE05E3H_PFAM.domtblout'; txt=Path(outdir)/'STAGE05E3H_PFAM.txt'
    with open(txt,'w') as fo:
        p=subprocess.run([hmmscan,'--cpu',str(threads),'--cut_ga','--domtblout',str(dom),hmm,str(query_faa)],stdout=fo,stderr=subprocess.PIPE,text=True,check=False)
    resource_rows.append({'resource':'Pfam-A','status':'SEARCHED_PASS' if p.returncode==0 else 'SEARCH_FAILED','path':hmm,'action':'' if p.returncode==0 else p.stderr[:300]})
    return p.returncode==0

def optional_mmseqs(query_faa, envvar, label, outdir, resource_rows, threads):
    db=os.environ.get(envvar,'').strip(); exe=shutil.which('mmseqs')
    if not db:
        resource_rows.append({'resource':label,'status':'NOT_PROVIDED','path':'','action':f'set {envvar}=/path/to/{label}_proteins.faa'})
        return False
    if not Path(db).exists() or not exe:
        resource_rows.append({'resource':label,'status':'UNAVAILABLE','path':db,'action':'provide FASTA and mmseqs'})
        return False
    out=Path(outdir)/f'STAGE05E3H_{label}.m8'; tmp=Path(outdir)/f'tmp_{label.lower()}'
    cmd=[exe,'easy-search',str(query_faa),db,str(out),str(tmp),'--threads',str(threads),'-s','7.5','-e','1e-3','--max-seqs','25','--format-output','query,target,fident,alnlen,qcov,tcov,evalue,bits']
    rc,_,err=run_cmd(cmd,timeout=7200)
    resource_rows.append({'resource':label,'status':'SEARCHED_PASS' if rc==0 else 'SEARCH_FAILED','path':db,'action':'' if rc==0 else err[:300]})
    return rc==0

def optional_remote_blast(reps, outdir, resource_rows):
    enabled=os.environ.get('E3H_REMOTE_BLAST','0')=='1'; exe=shutil.which('blastp')
    if not enabled:
        resource_rows.append({'resource':'NCBI_nr_remote_BLASTP','status':'DISABLED_DEFAULT','path':'','action':'set E3H_REMOTE_BLAST=1 to attempt; network failure is nonfatal'})
        return False
    if not exe:
        resource_rows.append({'resource':'NCBI_nr_remote_BLASTP','status':'BLASTP_NOT_FOUND','path':'','action':'install NCBI BLAST+ or leave external audit pending'})
        return False
    ok=True
    for fam,(hdr,seq) in reps.items():
        q=Path(outdir)/f'{fam}.rep.faa'; q.write_text(f'>{fam}|{hdr}\n{seq}\n')
        o=Path(outdir)/f'{fam}.remote_nr.tsv'
        cmd=[exe,'-remote','-db','nr','-query',str(q),'-evalue','1e-3','-max_target_seqs','20','-outfmt','6 qseqid sacc pident length qcovs evalue bitscore stitle','-out',str(o)]
        rc,_,err=run_cmd(cmd,timeout=1800); ok &= (rc==0)
        if rc!=0: (Path(outdir)/f'{fam}.remote_nr.ERROR.txt').write_text(err)
    resource_rows.append({'resource':'NCBI_nr_remote_BLASTP','status':'SEARCHED_PASS' if ok else 'SEARCH_FAILED_NONFATAL','path':'remote nr','action':'' if ok else 'inspect *.ERROR.txt; global novelty remains unresolved'})
    return ok

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--e3g',required=True); ap.add_argument('--e3f',required=True); ap.add_argument('--e3d',required=True); ap.add_argument('--project-root',required=True); ap.add_argument('--threads',type=int,default=16); ap.add_argument('--outdir',default='outputs')
    a=ap.parse_args(); root=Path(a.project_root).resolve(); out=Path(a.outdir).resolve(); out.mkdir(parents=True,exist_ok=True)
    work=out.parent/'work'/'resolved_inputs'; shutil.rmtree(work,ignore_errors=True); work.mkdir(parents=True)
    e3g=resolve(a.e3g,root); e3f=resolve(a.e3f,root); e3d=resolve(a.e3d,root)
    gdir=extract_tar(e3g,work/'e3g'); fdir=extract_tar(e3f,work/'e3f'); ddir=extract_tar(e3d,work/'e3d')
    policy=read_tsv(Path(__file__).resolve().parents[1]/'config/PRIORITY_FAMILIES.tsv')
    pol={r['family_id']:r for r in policy}
    reps_rows=read_tsv(gdir/'outputs/STAGE05E3G_FAMILY_REPRESENTATIVES.tsv'); reps_by={r['family_id']:r for r in reps_rows}
    if set(reps_by)!=set(FAMS): raise RuntimeError(f'E3G family set mismatch: {sorted(reps_by)}')
    existing=evidence_index_from_e3g(gdir)
    feat=[]; famsum=[]; all_unique={}; reps={}; motif_rows=[]; existing_rows=[]
    for fam in FAMS:
        pr=pol[fam]; rr=reps_by[fam]
        if rr['target_name']!=pr['representative_target'] or rr['aa_sha256']!=pr['expected_aa_sha256'] or int(rr['aa_length'])!=int(pr['expected_aa_length']):
            raise RuntimeError(f'{fam}: representative freeze mismatch')
        hfile=fdir/f'outputs/{fam}.all_local_homologs.faa'; hrec=fasta(hfile)
        if not hrec: raise RuntimeError(f'{fam}: no E3F homolog sequences')
        rseq=None
        for h,s in hrec:
            witness='witness=True' in h
            row=seq_features(fam,h,s,witness); feat.append(row); all_unique.setdefault(row['aa_sha256'],(fam,h,s))
            if h.split()[0]==rr['target_name']: rseq=s
        if rseq is None or sha256_text(rseq)!=rr['aa_sha256']: raise RuntimeError(f'{fam}: representative sequence not recovered')
        reps[fam]=(rr['target_name'],rseq)
        mfile=ddir/f'outputs/family_sequences/{fam}.repeat_core.fna'; mrec=fasta(mfile)
        motifs=sorted(set(s for _,s in mrec))
        if motifs!=[pr['exact_repeat_core']]: raise RuntimeError(f'{fam}: repeat core mismatch {motifs}')
        motif_rows.append({'family_id':fam,'exact_repeat_core':motifs[0],'repeat_length_nt':len(motifs[0]),'frozen_topology':pr['frozen_topology'],'source_fasta':str(mfile)})
        fr=[r for r in feat if r['family_id']==fam]
        def med(key):
            vals=sorted(float(x[key]) for x in fr); n=len(vals); return vals[n//2] if n%2 else (vals[n//2-1]+vals[n//2])/2
        famsum.append({'family_id':fam,'fam50_cluster':pr['fam50_cluster'],'priority_role':pr['priority_role'],'local_homolog_n':len(fr),'unique_aa_n':len(set(r['aa_sha256'] for r in fr)),'aa_length_min':min(int(r['aa_length']) for r in fr),'aa_length_max':max(int(r['aa_length']) for r in fr),'median_mean_kd_hydropathy':round(med('mean_kd_hydropathy'),6),'median_net_charge_proxy_pH7':round(med('net_charge_proxy_pH7'),3),'any_strong_hydrophobic_segment_heuristic':str(any(r['strong_hydrophobic_segment_heuristic']=='True' for r in fr)),'exact_repeat_core':pr['exact_repeat_core'],'repeat_length_nt':pr['repeat_length_nt'],'frozen_topology':pr['frozen_topology']})
        # representative-level existing evidence from E3G, now semantically classified.
        for layer in ['D1_BROAD_RVT','D4_R2_SWISSPROT','D4_R3_STRUCTURE_TIER','D4_R4_PROSTT5_FOLDSEEK','D4_R5_FULL3D']:
            es=[r for r in existing[fam] if r['stage_layer']==layer]
            lines=[r.get('matched_line_excerpt','') for r in es]
            existing_rows.append({'family_id':fam,'aa_sha256':rr['aa_sha256'],'stage_layer':layer,'semantic_status':classify_existing(lines,layer) if layer!='D1_BROAD_RVT' else ('MATCHED_D1_CONTEXT' if lines else 'NO_D1_CONTEXT_MATCH'),'matching_artifact_n':len(es),'source_files':' ; '.join(sorted(set(r.get('source_file','') for r in es)))})
    # query bundle
    qdir=out/'query_bundle'; qdir.mkdir(exist_ok=True)
    with open(qdir/'E3H_ALL_UNIQUE_LOCAL_HOMOLOGS.faa','w') as f:
        for sha,(fam,h,s) in sorted(all_unique.items()): f.write(f'>{fam}|AA_{sha}|{h.split()[0]}\n{s}\n')
    with open(qdir/'E3H_THREE_REPRESENTATIVES.faa','w') as f:
        for fam,(h,s) in reps.items(): f.write(f'>{fam}|{h}\n{s}\n')
    with open(qdir/'E3H_EXACT_REPEAT_CORES.fna','w') as f:
        for m in motif_rows: f.write(f">{m['family_id']}|{m['frozen_topology']}\n{m['exact_repeat_core']}\n")
    write_tsv(out/'STAGE05E3H_SEQUENCE_FEATURES.tsv',feat,list(feat[0]))
    write_tsv(out/'STAGE05E3H_FAMILY_INTRINSIC_SUMMARY.tsv',famsum,list(famsum[0]))
    write_tsv(out/'STAGE05E3H_REPEAT_CORE_LOCK.tsv',motif_rows,list(motif_rows[0]))
    write_tsv(out/'STAGE05E3H_EXISTING_EVIDENCE_SEMANTIC_SUMMARY.tsv',existing_rows,list(existing_rows[0]))
    all_homolog_evidence=recover_all_homolog_existing_evidence(root,feat,existing)
    write_tsv(out/'STAGE05E3H_ALL_HOMOLOG_EXISTING_EVIDENCE.tsv',all_homolog_evidence,list(all_homolog_evidence[0]))
    # optional resources
    resources=[]
    pf=optional_pfam(qdir/'E3H_ALL_UNIQUE_LOCAL_HOMOLOGS.faa',out,resources,a.threads)
    ph=optional_mmseqs(qdir/'E3H_ALL_UNIQUE_LOCAL_HOMOLOGS.faa','E3H_PHROG_FASTA','PHROG',out,resources,a.threads)
    vg=optional_mmseqs(qdir/'E3H_ALL_UNIQUE_LOCAL_HOMOLOGS.faa','E3H_VOG_FASTA','VOG',out,resources,a.threads)
    rb=optional_remote_blast(reps,out,resources)
    write_tsv(out/'STAGE05E3H_EXTERNAL_RESOURCE_STATUS.tsv',resources,['resource','status','path','action'])
    # literature/search query bundle; not executed here.
    qrows=[]
    for fam in FAMS:
        pr=pol[fam]
        qrows += [
          {'family_id':fam,'query_type':'exact_motif','query':f'"{pr["exact_repeat_core"]}" phage'},
          {'family_id':fam,'query_type':'representative_accession_architecture','query':f'"{reps_by[fam]["sequence_id"]}" repeat protein phage'},
          {'family_id':fam,'query_type':'fam50_internal_name','query':f'"{pr["fam50_cluster"]}"'}]
    write_tsv(out/'STAGE05E3H_PUBLIC_LITERATURE_SEARCH_QUERIES.tsv',qrows,['family_id','query_type','query'])
    global_complete = bool(pf and ph and vg and rb)
    # recommendations
    rec=[]
    for fam in FAMS:
        if fam=='ARCH250_0013': nxt='GLOBAL_CATALOG_AUDIT_THEN_SHORT_PROTEIN_STRUCTURE_AND_EVOLUTIONARY_COCONSERVATION'
        elif fam=='ARCH250_0010': nxt='GLOBAL_CATALOG_AUDIT_THEN_SHORT_PROTEIN_STRUCTURE_AND_LINEAGE_DEREPLICATION'
        else: nxt='GLOBAL_CATALOG_AUDIT_THEN_REPEAT_POSITIVE_VS_NEGATIVE_SUBCLADE_EVOLUTIONARY_COMPARISON'
        rec.append({'family_id':fam,'current_claim':'UNCHARACTERIZED_REPEAT_ASSOCIATED_PROTEIN_ARCHITECTURE','global_novelty_established':'False','recommended_next':nxt,'wetlab_priority':'HIGH' if fam=='ARCH250_0013' else ('MEDIUM_HIGH' if fam=='ARCH250_0010' else 'CONTROL')})
    write_tsv(out/'STAGE05E3H_DECISION_TABLE.tsv',rec,list(rec[0]))
    gate={
      'stage':'Stage05E3H','version':'v1.0.0','gate':'PASS_E3H_ACCELERATED_LOCAL_NOVELTY_FUNCTIONAL_CLOSURE_AND_QUERY_BUNDLE','gate_pass':True,
      'e3d_input_sha256':sha256_file(e3d),'e3f_input_sha256':sha256_file(e3f),'e3g_input_sha256':sha256_file(e3g),
      'priority_family_n':3,'local_unique_protein_sequence_n':len(all_unique),'repeat_core_lock_pass':True,'representative_freeze_pass':True,
      'pfam_search_performed':pf,'phrog_search_performed':ph,'vog_search_performed':vg,'remote_nr_blast_performed':rb,
      'global_novelty_complete':global_complete,
      'global_novelty_status':'COMPLETE_REQUIRED_RESOURCE_SET_SEARCHED_BUT_STILL_NOT_ABSOLUTE_PROOF' if global_complete else 'PENDING_EXTERNAL_DATABASE_CATALOG_AUDIT',
      'family_membership_changed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,'new_structure_prediction_performed':False,
      'claim_boundary':'Deterministic local/available-database closure. External database absence is never treated as proof of novelty or function.'}
    (out/'STAGE05E3H_GATE.json').write_text(json.dumps(gate,indent=2)+"\n")
    notes=['# Stage05E3H manuscript-facing notes','', 'E3H does not establish global novelty unless the prespecified external resource set is explicitly searched; even then, database no-hit is evidence of non-annotation, not proof of natural-world uniqueness.','']
    for r in famsum:
        notes += [f"## {r['family_id']}",f"- Exact frozen repeat core: `{r['exact_repeat_core']}` ({r['repeat_length_nt']} nt), topology={r['frozen_topology']}.",f"- Local homologs={r['local_homolog_n']}; unique AA={r['unique_aa_n']}; length range={r['aa_length_min']}-{r['aa_length_max']} aa.",f"- Strong-hydrophobic-segment heuristic in any local homolog: {r['any_strong_hydrophobic_segment_heuristic']} (descriptive heuristic only).",'']
    notes += ['## Claim boundary','Use `previously uncharacterized repeat-associated protein architecture` until an explicit global catalog/literature audit closes the stronger novelty wording.','', '## PB50 lock','The PB50 manuscript lock remains mandatory and unchanged.']
    (out/'STAGE05E3H_MANUSCRIPT_NOTES.md').write_text('\n'.join(notes)+'\n')
    print(json.dumps(gate,indent=2))

if __name__=='__main__': main()
