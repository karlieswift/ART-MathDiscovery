#!/usr/bin/env python3
import argparse,csv,hashlib,json,os,re,shutil,statistics,tarfile,gzip
from pathlib import Path
STAGE='Stage05E3D'; VERSION='v1.0.5'

AA_ALLOWED=set('ACDEFGHIKLMNPQRSTVWYBXZJUO*')

def safe_extract(src,dst):
    dst=Path(dst).resolve(); dst.mkdir(parents=True,exist_ok=True)
    with tarfile.open(src,'r:*') as tf:
        for m in tf.getmembers():
            t=(dst/m.name).resolve()
            if not str(t).startswith(str(dst)+os.sep) and t!=dst:
                raise RuntimeError(f'unsafe archive member {m.name}')
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
    with open(p,'r',encoding='utf-8',errors='replace',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(p,rows,fields):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,'w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore',lineterminator='\n')
        w.writeheader(); w.writerows(rows)

def truthy(v): return str(v).strip().lower() in {'1','true','t','yes','y'}
def num(v,default=0):
    try: return float(v)
    except: return default

def median(xs): return statistics.median(xs) if xs else ''

def norm_aa(s):
    s=''.join(str(s).split()).upper()
    if s.endswith('*'): s=s[:-1]
    return s

def aa_hash(s): return hashlib.sha256(norm_aa(s).encode()).hexdigest()

def parse_fasta_stream(handle):
    h=None; seq=[]
    for raw in handle:
        if isinstance(raw,bytes): raw=raw.decode('utf-8','replace')
        line=raw.rstrip('\n\r')
        if line.startswith('>'):
            if h is not None: yield h,norm_aa(''.join(seq))
            h=line[1:].strip(); seq=[]
        else:
            if h is not None: seq.append(line.strip())
    if h is not None: yield h,norm_aa(''.join(seq))

def candidate_protein_fasta_files(root):
    root=Path(root)
    strong_dirs=[]; other=[]
    keyparts=('STAGE05D','STAGE05E1','STAGE05E2','STAGE05C1','STAGE05C2')
    skip_parts=('ART_MATH_DISCOVERY_STAGE05E3D_','ART95_PHASE1','inputs/genomes','work/logan')
    for dp,dirs,files in os.walk(root):
        rel=str(Path(dp).relative_to(root)) if Path(dp)!=root else ''
        if any(x in rel for x in skip_parts):
            dirs[:] = []
            continue
        for fn in files:
            low=fn.lower()
            if not (low.endswith('.faa') or low.endswith('.fa') or low.endswith('.fasta') or low.endswith('.faa.gz') or low.endswith('.fa.gz') or low.endswith('.fasta.gz')):
                continue
            p=Path(dp)/fn
            try:
                if p.stat().st_size > 400_000_000: continue
            except OSError: continue
            # avoid obvious nucleotide/genome files when filename says so
            if any(z in low for z in ('.fna','genome','contig','assembly','nucleotide')):
                continue
            target=strong_dirs if any(k in rel for k in keyparts) else other
            target.append(p)
    # stable order: likely stage sources first, then filename hints, then path
    def score(p):
        n=p.name.lower()
        hints=sum(h in n for h in ('protein','orf','aa','unique','query','pred'))
        return (-hints,str(p))
    return sorted(strong_dirs,key=score)+sorted(other,key=score)

def recover_by_hash(project_root,wanted_hashes):
    wanted=set(wanted_hashes); found={}; scanned=[]
    for p in candidate_protein_fasta_files(project_root):
        if not wanted: break
        try:
            op=gzip.open if p.name.lower().endswith('.gz') else open
            with op(p,'rt',encoding='utf-8',errors='replace') as f:
                first_seen=False
                nrec=0
                for h,s in parse_fasta_stream(f):
                    nrec+=1
                    if not s or len(s)>10000: continue
                    # protein-safe check; skip obvious DNA-only records quickly
                    chars=set(s)
                    if not chars.issubset(AA_ALLOWED): continue
                    hv=hashlib.sha256(s.encode()).hexdigest()
                    if hv in wanted:
                        found[hv]={'sequence':s,'source_file':str(p),'source_header':h,'source_mode':'fasta_sha256_exact'}
                        wanted.remove(hv)
                if nrec: scanned.append({'path':str(p),'records_seen':nrec})
        except Exception:
            continue
    return found,wanted,scanned

def choose_predicted_witness(ctxrow,cand_rows):
    # 1) explicit AA hash in context, if present and valid for E2A witness rows
    hashes={r.get('aa_sha256','').strip() for r in cand_rows if r.get('aa_sha256','').strip()}
    for k,v in ctxrow.items():
        if 'aa_sha' in k.lower() and str(v).strip() in hashes:
            hv=str(v).strip(); rr=[r for r in cand_rows if r.get('aa_sha256','').strip()==hv]
            return rr[0],f'context_hash:{k}'
    # 2) explicit record/orf key in context matching E2A record_key
    record_values=[]
    for k,v in ctxrow.items():
        lk=k.lower(); sv=str(v).strip()
        if sv and (('record' in lk and 'key' in lk) or ('orf' in lk and ('id' in lk or 'key' in lk))):
            record_values.append((k,sv))
    for k,sv in record_values:
        rr=[r for r in cand_rows if r.get('record_key','').strip()==sv]
        if len(rr)==1: return rr[0],f'context_record:{k}'
    # 3) exactly one unique hash among all exact-primary witness rows for this candidate
    if len(hashes)==1:
        hv=next(iter(hashes)); rr=[r for r in cand_rows if r.get('aa_sha256','').strip()==hv]
        return rr[0],'candidate_unique_aa_sha256'
    return None,f'ambiguous_candidate_witness_hashes:{len(hashes)}'

def write_fasta(path,recs):
    with open(path,'w') as f:
        for rid,seq,meta in recs:
            f.write(f'>{rid} {meta}\n')
            for i in range(0,len(seq),80): f.write(seq[i:i+80]+'\n')

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
    expected_off_med={r['family_id']:r.get('official_witness_median_aa','') for r in q}
    priority=[f for f in fams if tiers[f]=='E3C1_PRIORITY_BIOLOGICAL_CLOSURE']
    if len(priority)!=3: raise RuntimeError(f'expected 3 priority families; got {priority}')

    witnesses=read_tsv(find1(e3a,'STAGE05E3A_E2A_WITNESS_PAIRS.tsv'))
    ctx=read_tsv(find1(e3a,'STAGE05E3A_CANDIDATE_OFFICIAL_CONTEXT.tsv'))

    # repeat cores: family-level exact E2A signatures, independent of official protein export
    repeat_by={f:[] for f in fams}; seen_rep={f:set() for f in fams}
    witness_by_candidate={}
    for r in witnesses:
        fam=r.get('family_id',''); cand=r.get('candidate_key','').strip()
        if fam not in repeat_by: continue
        core=re.sub('[^ACGTN]','',r.get('signature_repeat_core','').upper())
        if core and core not in seen_rep[fam]:
            seen_rep[fam].add(core); repeat_by[fam].append((cand,core,r.get('repeat_group_id',''),r.get('repeat_assay','')))
        if cand:
            witness_by_candidate.setdefault((fam,cand),[]).append(r)

    # E3A official support determines which E2A family-defining predicted witnesses are admitted.
    selected=[]; selection_diag=[]; seen_keys=set()
    for r in ctx:
        fam=r.get('family_id',''); cand=r.get('candidate_key','').strip()
        if fam not in expected_off or not cand or not truthy(r.get('has_exact_primary_joint_witness','')): continue
        pid=r.get('witness_official_protein_id','').strip()
        if not pid: continue
        key=(fam,cand)
        if key in seen_keys: continue
        seen_keys.add(key)
        cr=witness_by_candidate.get(key,[])
        chosen,route=choose_predicted_witness(r,cr)
        if chosen is None:
            selection_diag.append({'family_id':fam,'candidate_key':cand,'sequence_id':r.get('sequence_id',''),'official_protein_id':pid,'status':'UNRESOLVED','selection_route':route,'candidate_e2a_witness_rows_n':len(cr),'candidate_unique_aa_hashes_n':len({x.get('aa_sha256','').strip() for x in cr if x.get('aa_sha256','').strip()}),'record_key':'','aa_sha256':''})
            continue
        hv=chosen.get('aa_sha256','').strip()
        if not re.fullmatch(r'[0-9a-fA-F]{64}',hv):
            raise RuntimeError(f'{fam} {cand}: selected E2A witness has invalid aa_sha256={hv}')
        selected.append({'family_id':fam,'candidate_key':cand,'sequence_id':r.get('sequence_id',''),'official_protein_id':pid,'official_product':r.get('witness_official_product',''),'official_match_class':r.get('witness_official_match_class',''),'record_key':chosen.get('record_key',''),'aa_sha256':hv.lower(),'signature_FAM50_cluster':chosen.get('signature_FAM50_cluster',''),'selection_route':route})
        selection_diag.append({'family_id':fam,'candidate_key':cand,'sequence_id':r.get('sequence_id',''),'official_protein_id':pid,'status':'SELECTED','selection_route':route,'candidate_e2a_witness_rows_n':len(cr),'candidate_unique_aa_hashes_n':len({x.get('aa_sha256','').strip() for x in cr if x.get('aa_sha256','').strip()}),'record_key':chosen.get('record_key',''),'aa_sha256':hv.lower()})

    write_tsv(out/'STAGE05E3D_WITNESS_SELECTION_DIAGNOSTICS.tsv',selection_diag,['family_id','candidate_key','sequence_id','official_protein_id','status','selection_route','candidate_e2a_witness_rows_n','candidate_unique_aa_hashes_n','record_key','aa_sha256'])

    unresolved_selection=[x for x in selection_diag if x['status']!='SELECTED']
    wanted={x['aa_sha256'] for x in selected}
    recovered,missing_hashes,scanned=recover_by_hash(a.project_root,wanted)
    write_tsv(out/'STAGE05E3D_SEQUENCE_SOURCE_SCAN.tsv',scanned,['path','records_seen'])

    seqdir=out/'family_sequences'; seqdir.mkdir()
    byfam={f:[] for f in fams}; audit=[]; missing_rows=[]
    for x in selected:
        hv=x['aa_sha256']; rec=recovered.get(hv)
        if rec is None:
            missing_rows.append({**x,'reason':'aa_sha256_not_found_in_project_protein_fastas'})
            continue
        seq=rec['sequence']
        if hashlib.sha256(seq.encode()).hexdigest()!=hv:
            raise RuntimeError(f'hash mismatch after recovery for {hv}')
        y={**x,'aa_length':len(seq),'source_file':rec['source_file'],'source_header':rec['source_header'],'source_mode':rec['source_mode'],'sequence':seq}
        byfam[x['family_id']].append(y); audit.append(y)

    write_tsv(out/'STAGE05E3D_MISSING_PREDICTED_WITNESS_SEQUENCES.tsv',missing_rows,['family_id','candidate_key','sequence_id','official_protein_id','official_product','official_match_class','record_key','aa_sha256','signature_FAM50_cluster','selection_route','reason'])

    manifest=[]; sums=[]; mismatches=[]
    for fam in fams:
        rrecs=[]
        for i,(cand,seq,rg,assay) in enumerate(repeat_by[fam],1):
            rid=f'{fam}__REP{i:03d}__{re.sub("[^A-Za-z0-9_.:-]","_",cand)[:80]}'
            rrecs.append((rid,seq,f'candidate={cand} repeat_group={rg} assay={assay}'))
            manifest.append({'family_id':fam,'sequence_type':'repeat_core_dna','candidate_key':cand,'source_id':rg,'materialized_id':rid,'length':len(seq),'sha256':hashlib.sha256(seq.encode()).hexdigest(),'source_file':'E3A:E2A_WITNESS_PAIRS'})
        precs=[]
        for i,y in enumerate(sorted(byfam[fam],key=lambda z:z['candidate_key']),1):
            rid=f'{fam}__PREDWIT{i:03d}__{re.sub("[^A-Za-z0-9_.:-]","_",y["candidate_key"])[:80]}'
            precs.append((rid,y['sequence'],f'candidate={y["candidate_key"]} record_key={y["record_key"]} official_support={y["official_protein_id"]} FAM50={y["signature_FAM50_cluster"]}'))
            manifest.append({'family_id':fam,'sequence_type':'family_defining_predicted_witness_aa','candidate_key':y['candidate_key'],'source_id':y['record_key'],'materialized_id':rid,'length':len(y['sequence']),'sha256':y['aa_sha256'],'source_file':y['source_file']})
        write_fasta(seqdir/f'{fam}.repeat_core.fna',rrecs)
        write_fasta(seqdir/f'{fam}.predicted_witness.faa',precs)
        obs=len(precs); exp=expected_off[fam]; pred_med=median([len(y['sequence']) for y in byfam[fam]])
        if obs!=exp:
            mismatches.append({'family_id':fam,'expected_official_supported_witness_n':exp,'materialized_predicted_witness_n':obs,'e3c_official_witness_median_aa':expected_off_med[fam],'predicted_witness_median_aa':pred_med})
        sums.append({'family_id':fam,'e3c_review_tier':tiers[fam],'repeat_core_unique_n':len(rrecs),'expected_official_supported_witness_n':exp,'materialized_predicted_witness_n':obs,'materialized_unique_aa_n':len({y['aa_sha256'] for y in byfam[fam]}),'predicted_witness_median_aa':pred_med,'e3c_official_witness_median_aa':expected_off_med[fam],'count_matches_e3c':obs==exp,'priority_sequence_complete':(len(rrecs)>0 and obs==exp and obs>0) if fam in priority else ''})

    write_tsv(out/'STAGE05E3D_FAMILY_SEQUENCE_SUMMARY.tsv',sums,list(sums[0]))
    write_tsv(out/'STAGE05E3D_SEQUENCE_MANIFEST.tsv',manifest,['family_id','sequence_type','candidate_key','source_id','materialized_id','length','sha256','source_file'])
    aud_fields=['family_id','candidate_key','sequence_id','official_protein_id','official_product','official_match_class','record_key','aa_sha256','signature_FAM50_cluster','selection_route','aa_length','source_file','source_header','source_mode']
    write_tsv(out/'STAGE05E3D_PREDICTED_WITNESS_AUDIT.tsv',audit,aud_fields)
    write_tsv(out/'STAGE05E3D_E3C_CONCORDANCE_MISMATCHES.tsv',mismatches,['family_id','expected_official_supported_witness_n','materialized_predicted_witness_n','e3c_official_witness_median_aa','predicted_witness_median_aa'])

    priority_missing=[r['family_id'] for r in sums if r['family_id'] in priority and not truthy(r['priority_sequence_complete'])]
    total_expected=sum(expected_off.values()); total_obs=sum(r['materialized_predicted_witness_n'] for r in sums)
    gate_pass=(not unresolved_selection and not missing_rows and not mismatches and not priority_missing and total_obs==total_expected)
    gate='PASS_E3D_FAMILY_DEFINING_PREDICTED_WITNESS_SEQUENCES_HASH_MATERIALIZED' if gate_pass else 'HOLD_E3D_PREDICTED_WITNESS_SEQUENCE_RECOVERY'
    summary={
      'stage':STAGE,'version':VERSION,'n_families':len(fams),'priority_families':priority,'priority_missing_families':priority_missing,
      'total_repeat_core_unique_n':sum(r['repeat_core_unique_n'] for r in sums),'total_official_supported_witness_expected_n':total_expected,
      'total_predicted_witness_materialized_n':total_obs,'unresolved_witness_selection_n':len(unresolved_selection),'missing_aa_hash_n':len(missing_hashes),
      'missing_sequence_rows_n':len(missing_rows),'concordance_mismatch_n':len(mismatches),'protein_identity_source':'E2A family-defining predicted witness aa_sha256',
      'official_protein_fasta_required':False,'official_support_role':'E3A metadata/coordinate support only; not the family-defining sequence source',
      'v1_0_1_to_v1_0_4_official_fasta_materialization_superseded':True,'family_membership_changed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,
      'gate':gate,'gate_pass':gate_pass,'claim_boundary':'Materializes E2A family-defining predicted witness proteins by exact aa_sha256 for E3A official-supported candidates, plus exact E2A repeat cores. Official protein IDs are retained as support metadata but are not required to exist in the frozen neighborhood FASTA. No ART identity, novelty, function, or wet-lab claim.'}
    (out/'STAGE05E3D_SUMMARY.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    (out/'STAGE05E3D_GATE.json').write_text(json.dumps({'gate':gate,'gate_pass':gate_pass,'expected_family_n':a.expected_families,'observed_family_n':len(fams),'total_expected_n':total_expected,'total_materialized_n':total_obs,'unresolved_witness_selection_n':len(unresolved_selection),'missing_sequence_rows_n':len(missing_rows),'priority_missing_families':priority_missing,'concordance_mismatch_n':len(mismatches)},indent=2,sort_keys=True)+'\n')
    print(json.dumps(summary,indent=2,sort_keys=True))
    raise SystemExit(0 if gate_pass else 2)

if __name__=='__main__': main()
