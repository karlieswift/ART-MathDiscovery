#!/usr/bin/env python3
import argparse, csv, hashlib, json, os, re, shutil, statistics, tarfile
from pathlib import Path

STAGE='Stage05E3D'; VERSION='v1.0.1'
FAM_ALIASES=['family_id','family','arch_family_id','e2a_family_id']
CAND_ALIASES=['candidate_id','candidate','candidate_key','locus_id','region_id']
SEQID_ALIASES=['sequence_id','seq_id','accession','record_id','genome_id']
AA_NAME_RX=re.compile(r'(official.*aa|aa.*sequence|protein.*sequence|translation|predicted.*aa|orf.*aa)',re.I)
REP_NAME_RX=re.compile(r'(repeat.*(sequence|seq|core|consensus)|signature.*(sequence|seq)|core.*sequence)',re.I)
PROT_ID_RX=re.compile(r'(protein.*id|official.*protein|protein_accession|protein_acc)',re.I)

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
            if not str(t).startswith(str(dst)+os.sep) and t!=dst: raise RuntimeError(f'unsafe archive member {m.name}')
        tf.extractall(dst)

def resolve_input(src, project_root=None):
    src=Path(src)
    if src.exists():
        return src
    if project_root is None:
        raise FileNotFoundError(src)
    root=Path(project_root)
    if not root.exists():
        raise FileNotFoundError(f'project root not found: {root}')
    # Collector packages are often emitted inside their stage directory while the
    # command line historically referenced the project root. Recover only by exact
    # basename and only when the match is unique; never guess by prefix/version.
    hits=[p for p in root.rglob(src.name) if p.is_file() or p.is_dir()]
    if len(hits)==1:
        print(f'[RESOLVE] {src} -> {hits[0]}')
        return hits[0]
    if not hits:
        raise FileNotFoundError(f'{src}; recursive exact-basename search under {root} found 0 matches')
    raise RuntimeError(f'{src}; recursive exact-basename search under {root} found {len(hits)} matches: {hits[:12]}')

def materialize(src,dst,project_root=None):
    src=resolve_input(src, project_root); dst=Path(dst)
    if dst.exists(): shutil.rmtree(dst)
    if src.is_dir(): shutil.copytree(src,dst)
    else: safe_extract(src,dst)
    return dst

def find1(root,name):
    hits=list(Path(root).rglob(name))
    if len(hits)!=1: raise RuntimeError(f'expected exactly one {name}; found {len(hits)}: {hits[:8]}')
    return hits[0]

def read_tsv(p):
    with open(p,'r',encoding='utf-8',errors='replace',newline='') as f: return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(p,rows,fields=None):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    if fields is None:
        fields=[]
        for r in rows:
            for k in r:
                if k not in fields: fields.append(k)
    with open(p,'w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore',lineterminator='\n'); w.writeheader(); w.writerows(rows)

def first_col(row,aliases):
    lower={k.lower():k for k in row}
    for a in aliases:
        if a in lower and str(row.get(lower[a],'')).strip(): return lower[a]
    return None

def norm_dna(s): return re.sub('[^ACGTN]','',str(s).upper())
def norm_aa(s): return re.sub('[^A-Z]','',str(s).upper().replace('*',''))
def looks_dna(s):
    x=norm_dna(s); raw=re.sub(r'\s+','',str(s).upper())
    return len(x)>=6 and len(x)>=0.9*max(1,len(raw)) and re.fullmatch('[ACGTN]+',x) is not None

def looks_aa(s):
    x=norm_aa(s); raw=re.sub(r'\s+','',str(s).upper().replace('*',''))
    if len(x)<10 or len(x)<0.9*max(1,len(raw)): return False
    return re.fullmatch('[ACDEFGHIKLMNPQRSTVWYXBZUOJ]+',x) is not None and not re.fullmatch('[ACGTN]+',x)

def parse_fasta(p):
    recs=[]; h=None; seq=[]
    with open(p,'r',encoding='utf-8',errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                if h is not None: recs.append((h,''.join(seq).replace(' ','').upper()))
                h=line[1:].strip(); seq=[]
            else: seq.append(line.strip())
        if h is not None: recs.append((h,''.join(seq).replace(' ','').upper()))
    return recs

def fasta_id_aliases(header):
    vals=set()
    first=header.split()[0]
    vals.add(first); vals.add(first.split('|')[-1]); vals.add(first.split('|')[0])
    for tok in re.split(r'[|;,\s]+',header):
        tok=tok.strip()
        if tok: vals.add(tok)
    more=set()
    for v in vals:
        more.add(v.split('.')[0])
    return vals|more

def write_fasta(path,recs):
    with open(path,'w') as f:
        for rid,seq,meta in recs:
            f.write(f'>{rid} {meta}\n')
            for i in range(0,len(seq),80): f.write(seq[i:i+80]+'\n')

def med(vals): return statistics.median(vals) if vals else ''

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--e3c',required=True); ap.add_argument('--e3a',required=True); ap.add_argument('--project-root',required=True)
    ap.add_argument('--outdir',default='outputs'); ap.add_argument('--workdir',default='work'); ap.add_argument('--expected-families',type=int,default=12)
    args=ap.parse_args()
    out=Path(args.outdir); work=Path(args.workdir)
    if out.exists(): shutil.rmtree(out)
    if work.exists(): shutil.rmtree(work)
    out.mkdir(parents=True); work.mkdir(parents=True)
    e3c=materialize(args.e3c,work/'e3c',args.project_root); e3a=materialize(args.e3a,work/'e3a',args.project_root)
    q=read_tsv(find1(e3c,'STAGE05E3C_REVIEW_QUEUE.tsv'))
    if len(q)!=args.expected_families: raise RuntimeError(f'E3C family count {len(q)} != {args.expected_families}')
    fams=[r['family_id'] for r in q]; tiers={r['family_id']:r.get('e3c_review_tier','') for r in q}
    priority=[f for f in fams if tiers[f]=='E3C1_PRIORITY_BIOLOGICAL_CLOSURE']
    if len(priority)!=3: raise RuntimeError(f'expected 3 E3C1 families, got {priority}')

    candidate_files=[]
    for name in ['STAGE05E3A_E2A_CANDIDATES.tsv','STAGE05E3A_E2A_WITNESS_PAIRS.tsv']:
        try: candidate_files.append(find1(e3a,name))
        except Exception: pass
    if not candidate_files: raise RuntimeError('no frozen E2A candidate/witness tables found in E3A input')

    cand_to_fam={}; seqid_to_fams={}; tabrows=[]
    for p in candidate_files:
        rows=read_tsv(p); tabrows.append((p,rows))
        for r in rows:
            fc=first_col(r,FAM_ALIASES); cc=first_col(r,CAND_ALIASES); sc=first_col(r,SEQID_ALIASES)
            fam=str(r.get(fc,'')).strip() if fc else ''
            if fam not in fams: continue
            if cc: cand_to_fam[str(r.get(cc,'')).strip()]=fam
            if sc: seqid_to_fams.setdefault(str(r.get(sc,'')).strip(),set()).add(fam)

    repeats={f:[] for f in fams}; proteins={f:[] for f in fams}; diagnostics=[]
    for p,rows in tabrows:
        for i,r in enumerate(rows,2):
            fc=first_col(r,FAM_ALIASES); cc=first_col(r,CAND_ALIASES)
            fam=str(r.get(fc,'')).strip() if fc else ''
            if fam not in fams and cc: fam=cand_to_fam.get(str(r.get(cc,'')).strip(),'')
            if fam not in fams: continue
            rid=(str(r.get(cc,'')).strip() if cc else '') or f'{p.name}:{i}'
            for c,v in r.items():
                if v and REP_NAME_RX.search(c) and looks_dna(v):
                    seq=norm_dna(v); repeats[fam].append((rid,seq,f'{p.name}:{i}:{c}'))
                if v and AA_NAME_RX.search(c) and looks_aa(v):
                    seq=norm_aa(v); proteins[fam].append((rid,seq,f'{p.name}:{i}:{c}:DIRECT'))

    # official proteins and flexible candidate->protein mapping
    faa=find1(e3a,'STAGE05E3A_NEIGHBORHOOD_OFFICIAL_PROTEINS.faa')
    fasta_recs=parse_fasta(faa)
    alias_to_rec={}
    for h,s in fasta_recs:
        for a in fasta_id_aliases(h): alias_to_rec.setdefault(a,(h,s))

    context_files=[]
    for name in ['STAGE05E3A_CANDIDATE_OFFICIAL_CONTEXT.tsv','STAGE05E3A_NEIGHBORHOOD_OFFICIAL_CDS.tsv']:
        try: context_files.append(find1(e3a,name))
        except Exception: pass
    for p in context_files:
        rows=read_tsv(p)
        for i,r in enumerate(rows,2):
            fc=first_col(r,FAM_ALIASES); cc=first_col(r,CAND_ALIASES); sc=first_col(r,SEQID_ALIASES)
            fam=str(r.get(fc,'')).strip() if fc else ''
            if fam not in fams and cc: fam=cand_to_fam.get(str(r.get(cc,'')).strip(),'')
            if fam not in fams and sc:
                fs=seqid_to_fams.get(str(r.get(sc,'')).strip(),set())
                if len(fs)==1: fam=next(iter(fs))
            if fam not in fams: continue
            rid=(str(r.get(cc,'')).strip() if cc else '') or (str(r.get(sc,'')).strip() if sc else '') or f'{p.name}:{i}'
            # direct AA columns
            for c,v in r.items():
                if v and AA_NAME_RX.search(c) and looks_aa(v): proteins[fam].append((rid,norm_aa(v),f'{p.name}:{i}:{c}:DIRECT'))
            # explicit protein-id columns first
            tokens=[]
            for c,v in r.items():
                if v and PROT_ID_RX.search(c): tokens.extend(re.split(r'[|;,\s]+',str(v)))
            # then all tokens, to tolerate schema drift
            if not tokens:
                tokens=re.split(r'[|;,\s]+',' '.join(str(v) for v in r.values() if v))
            seen=set()
            for tok in tokens:
                tok=tok.strip()
                if not tok or tok in seen: continue
                seen.add(tok)
                rec=alias_to_rec.get(tok) or alias_to_rec.get(tok.split('.')[0])
                if rec:
                    h,s=rec
                    if looks_aa(s): proteins[fam].append((rid,norm_aa(s),f'{p.name}:{i}:FASTA:{h.split()[0]}'))

    # deduplicate exact sequences within family
    for d in (repeats,proteins):
        for fam in fams:
            seen=set(); u=[]
            for rec in d[fam]:
                if rec[1] in seen: continue
                seen.add(rec[1]); u.append(rec)
            d[fam]=u

    seqdir=out/'family_sequences'; seqdir.mkdir()
    man=[]; family_summary=[]
    for fam in fams:
        preps=[]; prots=[]
        for j,(rid,seq,src) in enumerate(repeats[fam],1):
            sid=f'{fam}__REP{j:03d}__{re.sub("[^A-Za-z0-9_.:-]","_",rid)[:80]}'
            preps.append((sid,seq,f'source={src}'))
            man.append({'family_id':fam,'sequence_type':'repeat_dna','record_id':rid,'materialized_id':sid,'length':len(seq),'sha256':hashlib.sha256(seq.encode()).hexdigest(),'source':src})
        for j,(rid,seq,src) in enumerate(proteins[fam],1):
            sid=f'{fam}__PROT{j:03d}__{re.sub("[^A-Za-z0-9_.:-]","_",rid)[:80]}'
            prots.append((sid,seq,f'source={src}'))
            man.append({'family_id':fam,'sequence_type':'protein_aa','record_id':rid,'materialized_id':sid,'length':len(seq),'sha256':hashlib.sha256(seq.encode()).hexdigest(),'source':src,'yxdd':bool(re.search(r'Y.DD',seq))})
        write_fasta(seqdir/f'{fam}.repeat.fna',preps)
        write_fasta(seqdir/f'{fam}.protein.faa',prots)
        family_summary.append({'family_id':fam,'e3c_review_tier':tiers[fam],'repeat_sequence_n':len(preps),'protein_sequence_n':len(prots),'repeat_median_nt':med([len(x[1]) for x in preps]),'protein_median_aa':med([len(x[1]) for x in prots]),'priority_sequence_complete':(len(preps)>0 and len(prots)>0) if fam in priority else ''})

    write_tsv(out/'STAGE05E3D_SEQUENCE_MANIFEST.tsv',man)
    write_tsv(out/'STAGE05E3D_FAMILY_SEQUENCE_SUMMARY.tsv',family_summary)
    complete_priority=[r['family_id'] for r in family_summary if r['family_id'] in priority and r['repeat_sequence_n']>0 and r['protein_sequence_n']>0]
    missing_priority=sorted(set(priority)-set(complete_priority))
    gate_pass=(len(missing_priority)==0)
    gate='PASS_E3D_PRIORITY_FAMILY_EXACT_SEQUENCES_MATERIALIZED' if gate_pass else 'HOLD_E3D_PRIORITY_SEQUENCE_MATERIALIZATION_INCOMPLETE'
    summary={'stage':STAGE,'version':VERSION,'n_families':len(fams),'priority_families':priority,'priority_complete_families':complete_priority,'priority_missing_families':missing_priority,'n_families_with_repeat_sequence':sum(r['repeat_sequence_n']>0 for r in family_summary),'n_families_with_protein_sequence':sum(r['protein_sequence_n']>0 for r in family_summary),'total_repeat_sequences':sum(r['repeat_sequence_n'] for r in family_summary),'total_protein_sequences':sum(r['protein_sequence_n'] for r in family_summary),'family_membership_changed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,'gate':gate,'gate_pass':gate_pass,'claim_boundary':'Exact sequence materialization only; no ART identity, novelty, function, or wet-lab claim.'}
    json.dump(summary,open(out/'STAGE05E3D_SUMMARY.json','w'),indent=2,sort_keys=True)
    json.dump({'gate':gate,'gate_pass':gate_pass,'expected_family_n':args.expected_families,'observed_family_n':len(fams),'priority_missing_families':missing_priority,'family_membership_changed':False,'dna_rank_recomputed':False},open(out/'STAGE05E3D_GATE.json','w'),indent=2,sort_keys=True)
    # schema diagnostics for any HOLD
    diags=[]
    for p in candidate_files+context_files:
        rows=read_tsv(p); diags.append({'file':str(p),'columns':' | '.join(rows[0].keys()) if rows else ''})
    write_tsv(out/'STAGE05E3D_SCHEMA_DIAGNOSTICS.tsv',diags,['file','columns'])
    print(json.dumps(summary,indent=2,sort_keys=True))
    raise SystemExit(0 if gate_pass else 2)

if __name__=='__main__': main()
