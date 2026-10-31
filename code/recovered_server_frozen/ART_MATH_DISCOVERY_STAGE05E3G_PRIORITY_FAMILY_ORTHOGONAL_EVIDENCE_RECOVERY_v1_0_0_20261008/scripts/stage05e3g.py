#!/usr/bin/env python3
from pathlib import Path
import argparse, tarfile, tempfile, shutil, hashlib, csv, json, re, os

FAMS=['ARCH250_0010','ARCH250_0013','ARCH250_0033']
STAGES=[
 ('D1_BROAD_RVT',['STAGE05D1','BROAD_RVT','RVT_ALL','RVT-']),
 ('D4_R1_FAMILY_COMPRESSION',['STAGE05D4','SEQUENCE_FAMILY_COMPRESSION','FAM50','NR90']),
 ('D4_R2_SWISSPROT',['SWISSPROT','SWISS-PROT','SENSITIVE_HOMOLOGY','NR90_SWISSPROT']),
 ('D4_R3_STRUCTURE_TIER',['STRUCTURE_QUERY_TIER','COMPACT_STRUCTURE','TIER_A','TIER_B']),
 ('D4_R4_PROSTT5_FOLDSEEK',['PROSTT5','FOLDSEEK','REMOTE_STRUCTURE']),
 ('D4_R5_FULL3D',['BOLTZ','FULL3D','FULL_3D','R5C','R5D'])]
TEXT_EXT={'.tsv','.csv','.json','.txt','.md','.log','.out','.yaml','.yml','.faa','.fasta','.fa'}
PATH_MARKERS=('STAGE05D1','STAGE05D4','DIVERGENT_RT','BOLTZ','FULL3D','FULL_3D','SWISSPROT','PROSTT5','FOLDSEEK','SEQUENCE_FAMILY')

def sha(s): return hashlib.sha256(s.encode()).hexdigest()

def materialize(src:Path,dst:Path):
    if src.is_dir():
        shutil.copytree(src,dst,dirs_exist_ok=True); return dst
    if not src.exists(): raise FileNotFoundError(src)
    dst.mkdir(parents=True,exist_ok=True)
    with tarfile.open(src,'r:*') as t:
        for m in t.getmembers():
            p=(dst/m.name).resolve()
            if not str(p).startswith(str(dst.resolve())): raise RuntimeError('unsafe tar path')
        t.extractall(dst)
    return dst

def parse_fasta(path):
    rec=[]; h=None; seq=[]
    for line in path.read_text(errors='replace').splitlines():
        if line.startswith('>'):
            if h is not None: rec.append((h,''.join(seq)))
            h=line[1:]; seq=[]
        else: seq.append(line.strip())
    if h is not None: rec.append((h,''.join(seq)))
    return rec

def pick_reps(mat):
    out=[]
    for fam in FAMS:
        hits=list(mat.rglob(f'{fam}.all_local_homologs.faa'))
        if len(hits)!=1: raise RuntimeError(f'{fam}: expected one E3F all_local_homologs.faa, found {len(hits)}')
        rec=parse_fasta(hits[0])
        wit=[x for x in rec if 'witness=True' in x[0]]
        if not wit: raise RuntimeError(f'{fam}: no witness=True record')
        h,s=wit[0]; target=h.split()[0]; acc=target.split('__PRODIGAL_')[0]
        out.append(dict(family_id=fam,target_name=target,sequence_id=acc,aa_length=len(s),aa_sha256=sha(s),aa_sequence=s,local_homolog_n=len(rec),witness_homolog_n=len(wit)))
    return out

def relevant_files(root):
    out=[]
    for dp,dns,fns in os.walk(root):
        p=Path(dp)
        # avoid scanning E3E/E3F/E3G self-generated trees and massive generic trees
        ps=str(p)
        if any(x in ps for x in ('STAGE05E3E','STAGE05E3F','STAGE05E3G')):
            dns[:] = [] if p!=root else dns
            continue
        if not any(m.lower() in ps.lower() for m in PATH_MARKERS):
            # still descend, because marker may appear below
            pass
        for fn in fns:
            fp=p/fn
            if fp.suffix.lower() not in TEXT_EXT: continue
            if not any(m.lower() in str(fp).lower() for m in PATH_MARKERS): continue
            try:
                if fp.stat().st_size>100_000_000: continue
            except OSError: continue
            out.append(fp)
    return out

def classify_stage(path, text):
    s=(str(path)+'\n'+text[:5000]).upper()
    scores=[]
    for stage,marks in STAGES:
        n=sum(m.upper() in s for m in marks)
        if n: scores.append((n,stage))
    return max(scores)[1] if scores else 'OTHER_D1_D4_RELATED'

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--e3f',required=True); ap.add_argument('--project-root',required=True); ap.add_argument('--out',required=True); ap.add_argument('--work',required=True)
    a=ap.parse_args(); out=Path(a.out); work=Path(a.work); out.mkdir(parents=True,exist_ok=True); work.mkdir(parents=True,exist_ok=True)
    e3f=Path(a.e3f)
    if not e3f.exists():
        hits=list(Path(a.project_root).rglob(e3f.name))
        if len(hits)==1:
            e3f=hits[0]; print(f'[RESOLVE] {Path(a.e3f)} -> {e3f}')
        else: raise FileNotFoundError(a.e3f)
    mat=materialize(e3f,work/'e3f')
    reps=pick_reps(mat)
    # E3F family summary for architecture class
    sums=list(mat.rglob('STAGE05E3F_FAMILY_EXPANSION_SUMMARY.tsv'))
    arch={}
    if len(sums)==1:
        with sums[0].open() as f:
            for r in csv.DictReader(f,delimiter='\t'): arch[r['family_id']]=r
    for r in reps:
        z=arch.get(r['family_id'],{})
        r['e3f_architecture_expansion_class']=z.get('e3f_architecture_expansion_class','')
    rf=out/'STAGE05E3G_FAMILY_REPRESENTATIVES.tsv'
    cols=['family_id','target_name','sequence_id','aa_length','aa_sha256','local_homolog_n','witness_homolog_n','e3f_architecture_expansion_class']
    with rf.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t'); w.writeheader(); w.writerows([{k:r[k] for k in cols} for r in reps])
    files=relevant_files(Path(a.project_root))
    evidence=[]
    for fp in files:
        try: txt=fp.read_text(errors='replace')
        except Exception: continue
        st=classify_stage(fp,txt)
        for r in reps:
            toks=[r['aa_sha256'],'AA_'+r['aa_sha256'],r['target_name']]
            found=[t for t in toks if t and t in txt]
            if not found: continue
            # retain up to 5 matching lines for compact audit
            lines=[]
            for ln in txt.splitlines():
                if any(t in ln for t in found):
                    lines.append(ln[:2000])
                    if len(lines)>=5: break
            evidence.append(dict(family_id=r['family_id'],stage_layer=st,source_file=str(fp),matched_tokens=';'.join(found),matched_line_excerpt=' || '.join(lines)))
    ef=out/'STAGE05E3G_EXISTING_ORTHOGONAL_EVIDENCE.tsv'
    ecols=['family_id','stage_layer','source_file','matched_tokens','matched_line_excerpt']
    with ef.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=ecols,delimiter='\t'); w.writeheader(); w.writerows(evidence)
    # summarize per stage and family
    ss=[]
    for r in reps:
        for stage,_ in STAGES:
            n=sum(1 for x in evidence if x['family_id']==r['family_id'] and x['stage_layer']==stage)
            status='MATCHED_EXISTING_ARTIFACT' if n else 'NOT_RECOVERED_FROM_EXISTING_ARTIFACTS'
            if stage in ('D4_R4_PROSTT5_FOLDSEEK','D4_R5_FULL3D') and r['aa_length']<180:
                status='OLD_D4_LENGTH_PRIORITY_EXCLUDED_LT180' if n==0 else status
            ss.append(dict(family_id=r['family_id'],aa_length=r['aa_length'],stage_layer=stage,status=status,matching_artifact_n=n))
    sf=out/'STAGE05E3G_STAGE_EVIDENCE_SUMMARY.tsv'
    scols=['family_id','aa_length','stage_layer','status','matching_artifact_n']
    with sf.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=scols,delimiter='\t'); w.writeheader(); w.writerows(ss)
    # short protein plan and decision support
    plan=[]
    for r in reps:
        if r['aa_length']<180:
            plan.append(dict(family_id=r['family_id'],aa_length=r['aa_length'],reason='below_old_D4_structure_priority_length_180aa',recommended_next='targeted_short_protein_global_sequence_annotation_then_structure_only_if_sequence_remains_unresolved',priority='HIGH' if r['family_id']=='ARCH250_0013' else 'MEDIUM'))
    pf=out/'STAGE05E3G_SHORT_PROTEIN_EXTENSION_PLAN.tsv'
    pcols=['family_id','aa_length','reason','recommended_next','priority']
    with pf.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=pcols,delimiter='\t'); w.writeheader(); w.writerows(plan)
    # Gate: reps must be exact and at least one D1/D4 artifact collection found project-wide. Evidence match per family not mandatory.
    stage_files_n=len(files)
    gate_pass=(len(reps)==3 and len({r['aa_sha256'] for r in reps})==3 and stage_files_n>0)
    gate={
      'stage':'Stage05E3G','version':'v1.0.0','gate':'PASS_E3G_PRIORITY_FAMILY_ORTHOGONAL_EVIDENCE_RECOVERY' if gate_pass else 'HOLD_E3G_EVIDENCE_RECOVERY_INCOMPLETE','gate_pass':gate_pass,
      'e3f_input':str(e3f),'e3f_input_sha256':hashlib.sha256(e3f.read_bytes()).hexdigest(),'priority_family_n':len(reps),'existing_d1_d4_text_artifacts_scanned_n':stage_files_n,'exact_evidence_rows_n':len(evidence),
      'new_external_search_performed':False,'new_structure_prediction_performed':False,'family_membership_changed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,
      'claim_boundary':'Existing downstream protein-evidence provenance recovery only; missing match is unresolved, not negative or novelty proof.'}
    (out/'STAGE05E3G_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    # manuscript notes
    byfam={r['family_id']:r for r in reps}
    lines=['# Stage05E3G manuscript-facing notes','', 'E3G is an evidence-recovery layer. It does not add a novelty claim.','']
    for fam in FAMS:
        r=byfam[fam]; lines += [f'## {fam}',f'- Representative: {r["target_name"]}; {r["aa_length"]} aa; SHA256={r["aa_sha256"]}.',f'- E3F architecture class: {r["e3f_architecture_expansion_class"]}.']
        st=[x for x in ss if x['family_id']==fam]
        for x in st: lines.append(f'- {x["stage_layer"]}: {x["status"]} ({x["matching_artifact_n"]} matching artifacts).')
        lines.append('')
    lines += ['## Interpretation boundary','Existing-stage absence does not establish biological absence. Short proteins <180 aa were not admitted to the old D4 coordinate-structure priority tier by design. Any targeted short-protein extension must be labeled post-E3F and must not feed back into DNA ranks or family membership.','','## PB50 lock','The PB50 manuscript lock remains mandatory and unchanged in scientific meaning.']
    (out/'STAGE05E3G_MANUSCRIPT_NOTES.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(gate,indent=2))

if __name__=='__main__': main()
