#!/usr/bin/env python3
import argparse,csv,hashlib,json,os,re,shutil,statistics,tarfile
from pathlib import Path
STAGE='Stage05E3D'; VERSION='v1.0.3'

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

def resolve_input(src,root):
    p=Path(src)
    if p.exists(): return p
    hits=[x for x in Path(root).rglob(p.name) if x.is_file() or x.is_dir()]
    if len(hits)==1:
        print(f'[RESOLVE] {p} -> {hits[0]}'); return hits[0]
    if not hits: raise FileNotFoundError(f'{p}; exact-basename search under {root} found 0')
    raise RuntimeError(f'{p}; exact-basename search under {root} found {len(hits)}: {hits[:12]}')

def materialize(src,dst,root):
    src=resolve_input(src,root); dst=Path(dst)
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

def write_tsv(p,rows,fields):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,'w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore',lineterminator='\n'); w.writeheader(); w.writerows(rows)

def parse_fasta(p):
    out=[]; h=None; seq=[]
    with open(p,encoding='utf-8',errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                if h is not None: out.append((h,''.join(seq).replace(' ','').upper().replace('*','')))
                h=line[1:].strip(); seq=[]
            else: seq.append(line.strip())
        if h is not None: out.append((h,''.join(seq).replace(' ','').upper().replace('*','')))
    return out

def aliases(h):
    vals=set(re.split(r'[|;,\s]+',h)); vals.add(h.split()[0])
    more=set()
    for v in vals:
        if v: more.add(v.split('.')[0])
    return {v for v in vals|more if v}

def write_fasta(path,recs):
    with open(path,'w') as f:
        for rid,seq,meta in recs:
            f.write(f'>{rid} {meta}\n')
            for i in range(0,len(seq),80): f.write(seq[i:i+80]+'\n')

def truthy(v): return str(v).strip().lower() in {'1','true','t','yes','y'}
def num(v,default=0):
    try: return float(v)
    except: return default

def median(xs): return statistics.median(xs) if xs else ''

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--e3c',required=True); ap.add_argument('--e3a',required=True); ap.add_argument('--project-root',required=True)
    ap.add_argument('--outdir',default='outputs'); ap.add_argument('--workdir',default='work'); ap.add_argument('--expected-families',type=int,default=12)
    a=ap.parse_args(); out=Path(a.outdir); work=Path(a.workdir)
    for p in (out,work):
        if p.exists(): shutil.rmtree(p)
        p.mkdir(parents=True)
    e3c=materialize(a.e3c,work/'e3c',a.project_root); e3a=materialize(a.e3a,work/'e3a',a.project_root)
    q=read_tsv(find1(e3c,'STAGE05E3C_REVIEW_QUEUE.tsv'))
    if len(q)!=a.expected_families: raise RuntimeError(f'E3C family count {len(q)} != {a.expected_families}')
    fams=[r['family_id'] for r in q]
    tiers={r['family_id']:r.get('e3c_review_tier','') for r in q}
    expected_off={r['family_id']:int(num(r.get('official_supported_witness_n',''),0)) for r in q}
    expected_med={r['family_id']:r.get('official_witness_median_aa','') for r in q}
    priority=[f for f in fams if tiers[f]=='E3C1_PRIORITY_BIOLOGICAL_CLOSURE']
    if len(priority)!=3: raise RuntimeError(f'expected 3 priority families; got {priority}')

    witnesses=read_tsv(find1(e3a,'STAGE05E3A_E2A_WITNESS_PAIRS.tsv'))
    ctx=read_tsv(find1(e3a,'STAGE05E3A_CANDIDATE_OFFICIAL_CONTEXT.tsv'))
    faa=find1(e3a,'STAGE05E3A_NEIGHBORHOOD_OFFICIAL_PROTEINS.faa')
    recs=parse_fasta(faa)
    amap={}
    for h,s in recs:
        for x in aliases(h):
            if x not in amap:
                amap[x]=(h,s)
            else:
                prev=amap[x]
                if prev is None:
                    # Alias was already proven ambiguous by an earlier conflicting record.
                    # Keep it ambiguous; do not dereference None on later records.
                    continue
                if prev[1] != s:
                    amap[x]=None

    repeat_by={f:[] for f in fams}; pred_hash_rows=[]
    seen_rep={f:set() for f in fams}; seen_ph=set()
    for r in witnesses:
        fam=r.get('family_id',''); core=re.sub('[^ACGTN]','',r.get('signature_repeat_core','').upper())
        if fam not in repeat_by: continue
        if core and core not in seen_rep[fam]:
            seen_rep[fam].add(core); repeat_by[fam].append((r.get('candidate_key',''),core,r.get('repeat_group_id',''),r.get('repeat_assay','')))
        aa=r.get('aa_sha256','').strip(); cand=r.get('candidate_key','').strip()
        if aa and (fam,cand,aa) not in seen_ph:
            seen_ph.add((fam,cand,aa)); pred_hash_rows.append({'family_id':fam,'candidate_key':cand,'sequence_id':r.get('sequence_id',''),'record_key':r.get('record_key',''),'aa_sha256':aa,'signature_FAM50_cluster':r.get('signature_FAM50_cluster','')})

    off_by={f:[] for f in fams}; seen_candidate={}
    audit=[]
    for r in ctx:
        fam=r.get('family_id',''); cand=r.get('candidate_key','').strip()
        if fam not in off_by or not truthy(r.get('has_exact_primary_joint_witness','')): continue
        pid=r.get('witness_official_protein_id','').strip()
        mclass=r.get('witness_official_match_class','').strip()
        if not pid: continue
        if pid in amap:
            rec=amap[pid]
        else:
            rec=amap.get(pid.split('.')[0])
        if rec is None: raise RuntimeError(f'{fam} {cand}: official witness protein id {pid} absent/ambiguous in frozen official FASTA')
        h,s=rec
        key=(fam,cand)
        if key in seen_candidate and seen_candidate[key][0]!=pid:
            raise RuntimeError(f'{fam} {cand}: conflicting official witness proteins {seen_candidate[key][0]} vs {pid}')
        seen_candidate[key]=(pid,s)
        off_by[fam].append((cand,pid,s,mclass,h))
        audit.append({'family_id':fam,'candidate_key':cand,'sequence_id':r.get('sequence_id',''),'witness_official_match_class':mclass,'witness_official_protein_id':pid,'witness_official_product':r.get('witness_official_product',''),'aa_length':len(s),'aa_sha256':hashlib.sha256(s.encode()).hexdigest()})

    # one official witness record per candidate; duplicates in context are not silently multiplied
    for fam in fams:
        d={}
        for x in off_by[fam]:
            if x[0] in d and d[x[0]][2]!=x[2]: raise RuntimeError(f'{fam} {x[0]}: nonidentical duplicated witness sequence')
            d[x[0]]=x
        off_by[fam]=list(d.values())

    seqdir=out/'family_sequences'; seqdir.mkdir()
    sums=[]; mismatches=[]; manifest=[]
    for fam in fams:
        rrecs=[]
        for i,(cand,seq,rg,assay) in enumerate(repeat_by[fam],1):
            rid=f'{fam}__REP{i:03d}__{re.sub("[^A-Za-z0-9_.:-]","_",cand)[:80]}'
            rrecs.append((rid,seq,f'candidate={cand} repeat_group={rg} assay={assay}'))
            manifest.append({'family_id':fam,'sequence_type':'repeat_core_dna','candidate_key':cand,'source_id':rg,'materialized_id':rid,'length':len(seq),'sha256':hashlib.sha256(seq.encode()).hexdigest()})
        precs=[]
        for i,(cand,pid,seq,mclass,h) in enumerate(off_by[fam],1):
            rid=f'{fam}__OFFWIT{i:03d}__{re.sub("[^A-Za-z0-9_.:-]","_",cand)[:80]}'
            precs.append((rid,seq,f'candidate={cand} protein_id={pid} match_class={mclass}'))
            manifest.append({'family_id':fam,'sequence_type':'official_witness_protein_aa','candidate_key':cand,'source_id':pid,'materialized_id':rid,'length':len(seq),'sha256':hashlib.sha256(seq.encode()).hexdigest()})
        write_fasta(seqdir/f'{fam}.repeat_core.fna',rrecs)
        write_fasta(seqdir/f'{fam}.official_witness.faa',precs)
        obs=len(precs); exp=expected_off[fam]; med=median([len(x[1]) for x in precs])
        medexp=expected_med[fam]
        count_ok=(obs==exp)
        med_ok=True
        if exp>0 and str(medexp).strip(): med_ok=abs(float(med)-float(medexp))<1e-9
        if not count_ok or not med_ok: mismatches.append({'family_id':fam,'expected_official_witness_n':exp,'observed_official_witness_n':obs,'expected_median_aa':medexp,'observed_median_aa':med})
        sums.append({'family_id':fam,'e3c_review_tier':tiers[fam],'repeat_core_unique_n':len(rrecs),'expected_official_witness_n':exp,'official_witness_record_n':obs,'official_witness_unique_aa_n':len({x[2] for x in off_by[fam]}),'official_witness_median_aa':med,'e3c_expected_median_aa':medexp,'official_count_matches_e3c':count_ok,'median_matches_e3c':med_ok,'priority_sequence_complete':(len(rrecs)>0 and obs==exp and obs>0) if fam in priority else ''})

    write_tsv(out/'STAGE05E3D_FAMILY_SEQUENCE_SUMMARY.tsv',sums,list(sums[0]))
    write_tsv(out/'STAGE05E3D_SEQUENCE_MANIFEST.tsv',manifest,['family_id','sequence_type','candidate_key','source_id','materialized_id','length','sha256'])
    write_tsv(out/'STAGE05E3D_OFFICIAL_WITNESS_AUDIT.tsv',audit,['family_id','candidate_key','sequence_id','witness_official_match_class','witness_official_protein_id','witness_official_product','aa_length','aa_sha256'])
    write_tsv(out/'STAGE05E3D_PREDICTED_WITNESS_AA_HASHES.tsv',pred_hash_rows,['family_id','candidate_key','sequence_id','record_key','aa_sha256','signature_FAM50_cluster'])
    write_tsv(out/'STAGE05E3D_E3C_CONCORDANCE_MISMATCHES.tsv',mismatches,['family_id','expected_official_witness_n','observed_official_witness_n','expected_median_aa','observed_median_aa'])
    priority_missing=[r['family_id'] for r in sums if r['family_id'] in priority and not truthy(r['priority_sequence_complete'])]
    total_expected=sum(expected_off.values()); total_obs=sum(r['official_witness_record_n'] for r in sums)
    gate_pass=(not mismatches and not priority_missing and total_obs==total_expected)
    gate='PASS_E3D_EXACT_WITNESS_SEQUENCES_CONCORDANT_WITH_E3C' if gate_pass else 'HOLD_E3D_EXACT_WITNESS_SEQUENCE_CONCORDANCE'
    summary={'stage':STAGE,'version':VERSION,'n_families':len(fams),'priority_families':priority,'priority_missing_families':priority_missing,'total_repeat_core_unique_n':sum(r['repeat_core_unique_n'] for r in sums),'total_official_witness_expected_n':total_expected,'total_official_witness_materialized_n':total_obs,'total_predicted_witness_hash_records':len(pred_hash_rows),'concordance_mismatch_n':len(mismatches),'neighborhood_proteins_materialized':False,'v1_0_1_neighborhood_protein_output_superseded':True,'family_membership_changed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,'gate':gate,'gate_pass':gate_pass,'claim_boundary':'Only exact E2A repeat cores and E3A official proteins explicitly linked to exact-primary witness ORFs are materialized. Whole-neighborhood proteins are excluded. No ART identity, novelty, function, or wet-lab claim.'}
    (out/'STAGE05E3D_SUMMARY.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    (out/'STAGE05E3D_GATE.json').write_text(json.dumps({'gate':gate,'gate_pass':gate_pass,'expected_family_n':a.expected_families,'observed_family_n':len(fams),'total_official_witness_expected_n':total_expected,'total_official_witness_materialized_n':total_obs,'priority_missing_families':priority_missing,'concordance_mismatch_n':len(mismatches),'v1_0_1_superseded':True},indent=2,sort_keys=True)+'\n')
    print(json.dumps(summary,indent=2,sort_keys=True))
    raise SystemExit(0 if gate_pass else 2)
if __name__=='__main__': main()
