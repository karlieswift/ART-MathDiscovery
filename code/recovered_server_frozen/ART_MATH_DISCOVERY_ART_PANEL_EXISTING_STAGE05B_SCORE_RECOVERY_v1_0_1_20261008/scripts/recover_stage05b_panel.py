#!/usr/bin/env python3
import argparse, csv, hashlib, json, os, re, shutil, tarfile
from pathlib import Path

TARGETS=['OR770614.1','OR836606.1','LC680885.1','MN091626.1','MZ779063.1']
EXPECTED_STAGE05B_SHA='4fc4f1f09ed92f4e0be9817b14c0c14fe7f9eafd56cc0918cad5e6a73a054ebd'


def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def safe_extract(src,dst):
    dst=Path(dst).resolve(); dst.mkdir(parents=True,exist_ok=True)
    with tarfile.open(src,'r:*') as tf:
        for m in tf.getmembers():
            t=(dst/m.name).resolve()
            if not str(t).startswith(str(dst)+os.sep) and t!=dst:
                raise RuntimeError(f'unsafe archive member: {m.name}')
        tf.extractall(dst)

def resolve(src, root, exact_basename=None):
    if src:
        p=Path(src)
        if p.exists(): return p
        name=p.name
    else:
        name=exact_basename
    hits=[p for p in Path(root).rglob(name) if p.is_file() or p.is_dir()]
    if len(hits)==1:
        print(f'[RESOLVE] {name} -> {hits[0]}')
        return hits[0]
    if not hits: raise FileNotFoundError(f'{name}: 0 exact-basename matches under {root}')
    raise RuntimeError(f'{name}: {len(hits)} exact-basename matches under {root}: {hits[:12]}')

def materialize(src,dst):
    src=Path(src); dst=Path(dst)
    if dst.exists(): shutil.rmtree(dst)
    if src.is_dir(): shutil.copytree(src,dst)
    else: safe_extract(src,dst)
    return dst

def read_tsv(p):
    with open(p,'r',encoding='utf-8',errors='replace',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(p,rows,fields):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,'w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore',lineterminator='\n')
        w.writeheader(); w.writerows(rows)

def find_prov(root):
    hits=[p for p in Path(root).rglob('ART_EXACT_SOURCE_FASTA_PROVENANCE.tsv') if p.is_file()]
    if not hits:
        raise RuntimeError('expected provenance table; found 0')
    # Collector packages intentionally contain byte-identical copies in inputs/ and
    # outputs/. Accept duplication only when the bytes are identical; otherwise HOLD.
    by={}
    for p in hits:
        by.setdefault(sha256_file(p),[]).append(p)
    if len(by)!=1:
        raise RuntimeError('multiple non-identical provenance tables: '+str([str(x) for x in hits]))
    group=next(iter(by.values()))
    # Prefer outputs/ as the frozen generated evidence, then the shortest path.
    group=sorted(group,key=lambda p:(0 if 'outputs' in p.parts else 1,len(str(p)),str(p)))
    if len(group)>1:
        print(f'[RESOLVE] {len(group)} byte-identical provenance copies -> {group[0]}')
    return group[0]

def find_fasta(root,acc):
    hits=[p for p in Path(root).rglob(f'{acc}.frozen.fasta') if p.is_file()]
    if len(hits)!=1: raise RuntimeError(f'{acc}: expected one frozen fasta; found {len(hits)}')
    return hits[0]

def acc_from_rows(rows):
    if not rows: return None
    cols=rows[0].keys()
    candidates=[c for c in cols if c.lower() in ('sequence_id','accession','seq_id','genome_id','record_id')]
    for c in candidates:
        vals={str(r.get(c,'')).strip() for r in rows if str(r.get(c,'')).strip()}
        if len(vals)==1:
            return next(iter(vals))
    return None

def file_matches_acc(p,acc):
    if acc in str(p): return True
    try:
        if p.suffix=='.tsv':
            rows=read_tsv(p)
            a=acc_from_rows(rows)
            if a==acc: return True
        elif p.suffix=='.json':
            txt=p.read_text(errors='replace')
            if acc in txt: return True
    except Exception:
        pass
    return False

def unique_identical(paths):
    if not paths: return None
    by={}
    for p in paths: by.setdefault(sha256_file(p),[]).append(p)
    if len(by)>1:
        raise RuntimeError('multiple non-identical candidate files: '+str([str(x) for x in paths[:12]]))
    return sorted(paths,key=lambda x:(len(str(x)),str(x)))[0]

def find_evidence(root,acc,kind):
    root=Path(root)
    if kind=='union': patterns=['STAGE03A_UNION_REGIONS.tsv','*UNION_REGIONS.tsv']
    elif kind=='proposal': patterns=['STAGE03A_SELECTED_PROPOSAL_WINDOWS.tsv','*SELECTED*PROPOSAL*.tsv']
    elif kind=='completion': patterns=['*COMPLETION*.json','*FINAL_GATE.json','*GATE*.json']
    else: raise ValueError(kind)
    paths=[]
    seen=set()
    for pat in patterns:
        for p in root.rglob(pat):
            if p.is_file() and p not in seen and file_matches_acc(p,acc):
                paths.append(p); seen.add(p)
    return unique_identical(paths)

def famtext(row):
    for c in row:
        if c.lower() in ('support_families','families','support_family','family_support'):
            return str(row.get(c,''))
    vals=[]
    for c,v in row.items():
        cl=c.lower()
        if ('family' in cl or 'support' in cl) and v: vals.append(str(v))
    return ' '.join(vals)

def is_dual(row):
    s=famtext(row).upper()
    return ('ACNV' in s) and ('RECURRENCE' in s)

def rank_value(row):
    for c in row:
        if c.lower() in ('region_rank','rank','priority_rank','region_priority','priority'):
            try: return float(row[c])
            except Exception: pass
    return None

def proposal_family(row):
    for c in row:
        if c.lower() in ('family','score_family','proposal_family','support_family'):
            return str(row.get(c,'')).upper()
    return ' '.join(str(v) for v in row.values()).upper()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--project-root',required=True)
    ap.add_argument('--exact-inputs',required=True)
    ap.add_argument('--stage05b',default='')
    ap.add_argument('--outdir',default='outputs')
    ap.add_argument('--workdir',default='work')
    args=ap.parse_args()
    out=Path(args.outdir); work=Path(args.workdir)
    if out.exists(): shutil.rmtree(out)
    if work.exists(): shutil.rmtree(work)
    out.mkdir(parents=True); work.mkdir(parents=True)

    inp=resolve(args.exact_inputs,args.project_root)
    inp_mat=materialize(inp,work/'exact_inputs')
    prov=read_tsv(find_prov(inp_mat))
    pmap={r['accession']:r for r in prov}
    if sorted(pmap)!=sorted(TARGETS):
        raise RuntimeError(f'exact-source panel accessions mismatch: {sorted(pmap)}')
    source_rows=[]
    for acc in TARGETS:
        f=find_fasta(inp_mat,acc); h=sha256_file(f)
        if h!=pmap[acc]['sha256']:
            raise RuntimeError(f'{acc}: frozen FASTA hash mismatch {h} != {pmap[acc]["sha256"]}')
        source_rows.append({'accession':acc,'fasta_sha256':h,'length_bp':pmap[acc].get('length_bp',''),'source_provenance_verified':True})

    s5=resolve(args.stage05b,args.project_root,'ART_STAGE05B_NEW919_SCORE_RESULTS_20261006_081746.tar.gz')
    if s5.is_file():
        h=sha256_file(s5)
        if s5.name=='ART_STAGE05B_NEW919_SCORE_RESULTS_20261006_081746.tar.gz' and h!=EXPECTED_STAGE05B_SHA:
            raise RuntimeError(f'Stage05B archive SHA256 mismatch: {h}')
    else: h='DIRECTORY_INPUT'
    s5mat=materialize(s5,work/'stage05b')

    panel=out/'per_genome'; panel.mkdir()
    sums=[]; manifest=[]
    for acc in TARGETS:
        u=find_evidence(s5mat,acc,'union')
        p=find_evidence(s5mat,acc,'proposal')
        g=find_evidence(s5mat,acc,'completion')
        if not u or not p or not g:
            raise RuntimeError(f'{acc}: incomplete Stage05B evidence union={u} proposal={p} completion={g}')
        ad=panel/acc; ad.mkdir()
        for label,src in [('UNION_REGIONS.tsv',u),('SELECTED_PROPOSALS.tsv',p),('COMPLETION_OR_GATE.json',g)]:
            dst=ad/label; shutil.copy2(src,dst)
            manifest.append({'accession':acc,'evidence_type':label,'source_path':str(src),'copied_path':str(dst),'sha256':sha256_file(dst)})
        ur=read_tsv(u); pr=read_tsv(p)
        dual=[r for r in ur if is_dual(r)]
        ranks=[rank_value(r) for r in dual]; ranks=[x for x in ranks if x is not None]
        acnv=sum('ACNV' in proposal_family(r) for r in pr)
        rec=sum('RECURRENCE' in proposal_family(r) for r in pr)
        sums.append({'accession':acc,'all_union_regions_n':len(ur),'dual_arm_regions_n':len(dual),'best_dual_region_rank':min(ranks) if ranks else '', 'selected_proposals_n':len(pr),'selected_acnv_proposals_n':acnv,'selected_recurrence_proposals_n':rec,'stage05b_score_evidence_recovered':True})

    write_tsv(out/'ART_PANEL_SOURCE_PROVENANCE_VERIFIED.tsv',source_rows,['accession','fasta_sha256','length_bp','source_provenance_verified'])
    write_tsv(out/'ART_PANEL_STAGE05B_SCORE_RECOVERY_SUMMARY.tsv',sums,['accession','all_union_regions_n','dual_arm_regions_n','best_dual_region_rank','selected_proposals_n','selected_acnv_proposals_n','selected_recurrence_proposals_n','stage05b_score_evidence_recovered'])
    write_tsv(out/'ART_PANEL_RECOVERED_FILE_MANIFEST.tsv',manifest,['accession','evidence_type','source_path','copied_path','sha256'])
    gate={
      'stage':'ART_PANEL_STAGEA1_EXISTING_SCORE_RECOVERY','version':'v1.0.1','gate':'PASS_ART_PANEL_FIVE_EXISTING_STAGE05B_FROZEN_SCORE_BOARDS_RECOVERED','gate_pass':True,
      'target_accessions':TARGETS,'n_target':len(TARGETS),'n_recovered':len(sums),'stage05b_input':str(s5),'stage05b_input_sha256':h,
      'score_recomputed':False,'raw_features_recomputed':False,'dna_rank_recomputed':False,'retuning_authorized':False,'truth_coordinates_used':False,'repeat_motif_used':False,
      'scientific_boundary':'These five genomes were already part of the frozen 919-genome Stage05B scoring universe. This stage recovers their existing frozen proposal/union evidence and does not constitute a new blind validation run or a rescore.'
    }
    (out/'ART_PANEL_STAGEA1_GATE.json').write_text(json.dumps(gate,indent=2,sort_keys=True)+'\n')
    print(json.dumps(gate,indent=2,sort_keys=True))

if __name__=='__main__': main()
