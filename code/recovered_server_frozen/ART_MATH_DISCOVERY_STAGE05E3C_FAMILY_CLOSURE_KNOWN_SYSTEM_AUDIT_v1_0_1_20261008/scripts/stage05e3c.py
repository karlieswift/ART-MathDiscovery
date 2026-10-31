#!/usr/bin/env python3
import argparse, csv, gzip, hashlib, json, math, os, re, shutil, statistics, subprocess, sys, tarfile, tempfile
from pathlib import Path

VERSION = "v1.0.1"
STAGE = "Stage05E3C"

SEQ_COL_HINTS = ["official_aa_sequence", "aa_sequence", "protein_sequence", "translation", "protein_aa", "aa"]
REPEAT_COL_HINTS = ["repeat_sequence", "repeat_core", "repeat_consensus", "signature_sequence", "signature", "core_sequence"]
ID_COL_HINTS = ["candidate_id", "sequence_id", "seq_id", "record_id", "protein_id", "accession"]


def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()

def safe_extract(archive, dest):
    dest=Path(dest).resolve(); dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, 'r:*') as tf:
        for m in tf.getmembers():
            target=(dest/m.name).resolve()
            if not str(target).startswith(str(dest)+os.sep) and target != dest:
                raise RuntimeError(f"unsafe archive member: {m.name}")
        tf.extractall(dest)

def materialize_input(src, dest):
    src=Path(src)
    if not src.exists(): raise FileNotFoundError(src)
    dest=Path(dest)
    if src.is_dir():
        if dest.exists(): shutil.rmtree(dest)
        shutil.copytree(src, dest)
    else:
        if dest.exists(): shutil.rmtree(dest)
        safe_extract(src, dest)
    return dest

def find_basename(root, basename):
    hits=list(Path(root).rglob(basename))
    if len(hits)!=1:
        raise RuntimeError(f"expected exactly one {basename}, found {len(hits)}: {hits[:5]}")
    return hits[0]

def read_tsv(path):
    with open(path, 'r', encoding='utf-8', errors='replace', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(path, rows, fields=None):
    path=Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields=[]
        for r in rows:
            for k in r:
                if k not in fields: fields.append(k)
    with open(path,'w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore',lineterminator='\n')
        w.writeheader(); w.writerows(rows)

def num(x, default=0.0):
    try:
        if x is None or str(x).strip()=="": return default
        return float(x)
    except Exception: return default

def intval(x, default=0):
    try: return int(float(str(x)))
    except Exception: return default

def median_or_blank(vals):
    vals=[v for v in vals if v is not None]
    return statistics.median(vals) if vals else ""

def valid_aa(s):
    s=re.sub(r'\s+','',s.upper())
    return bool(s) and re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWYXBZUOJ*.-]+',s) is not None

def normalize_aa(s): return re.sub(r'[^A-Z]','',s.upper().replace('*',''))
def valid_dna(s):
    s=re.sub(r'\s+','',s.upper())
    return bool(s) and re.fullmatch(r'[ACGTN]+',s) is not None

def kmer_set(s,k):
    s=re.sub(r'[^A-Z]','',s.upper())
    return {s[i:i+k] for i in range(max(0,len(s)-k+1))}

def jaccard(a,b,k):
    A=kmer_set(a,k); B=kmer_set(b,k)
    if not A and not B: return 1.0
    if not A or not B: return 0.0
    return len(A&B)/len(A|B)

def pairwise_median_jaccard(seqs,k):
    vals=[]
    for i in range(len(seqs)):
        for j in range(i+1,len(seqs)):
            vals.append(jaccard(seqs[i],seqs[j],k))
    return median_or_blank(vals)

def load_terms(path):
    rows=read_tsv(path); out=[]
    for r in rows:
        try: rx=re.compile(r['pattern'], re.I)
        except re.error: rx=re.compile(re.escape(r['pattern']), re.I)
        out.append((r['term_id'],r['category'],r['pattern'],rx))
    return out

def scan_text(text, terms):
    hits=[]
    for tid,cat,pat,rx in terms:
        if rx.search(text or ''): hits.append((tid,cat,pat))
    return hits

def inventory(root):
    rows=[]
    for p in sorted(Path(root).rglob('*')):
        if p.is_file():
            rows.append({"relative_path":str(p.relative_to(root)),"size_bytes":p.stat().st_size,"sha256":sha256_file(p)})
    return rows

def collect_family_linked_records(roots, family_ids):
    aa={f:[] for f in family_ids}; rep={f:[] for f in family_ids}; text={f:[] for f in family_ids}; provenance=[]
    famset=set(family_ids)
    for root in roots:
        for p in sorted(Path(root).rglob('*.tsv')):
            try: rows=read_tsv(p)
            except Exception: continue
            if not rows: continue
            cols=set(rows[0])
            famcol=next((c for c in cols if c.lower() in ('family_id','family','arch_family_id','e2a_family_id')),None)
            if not famcol: continue
            for idx,r in enumerate(rows, start=2):
                fam=str(r.get(famcol,'')).strip()
                if fam not in famset: continue
                joined=' | '.join(str(v) for v in r.values() if v not in (None,''))
                text[fam].append(joined)
                idval=''
                for c in ID_COL_HINTS:
                    if c in r and str(r[c]).strip(): idval=str(r[c]).strip(); break
                for c in SEQ_COL_HINTS:
                    if c in r and valid_aa(str(r[c])):
                        seq=normalize_aa(str(r[c])); aa[fam].append((idval or f"{p.name}:{idx}",seq,str(p)))
                for c in REPEAT_COL_HINTS:
                    if c in r and valid_dna(str(r[c])):
                        seq=re.sub(r'[^ACGTN]','',str(r[c]).upper()); rep[fam].append((idval or f"{p.name}:{idx}",seq,str(p)))
                provenance.append({"family_id":fam,"source_file":str(p),"source_row":idx})
        # FASTA headers containing explicit family IDs
        for p in sorted(list(Path(root).rglob('*.faa'))+list(Path(root).rglob('*.fa'))+list(Path(root).rglob('*.fasta'))):
            try:
                cur=None; seq=[]
                def flush():
                    nonlocal cur,seq
                    if cur is None: return
                    s=''.join(seq).strip().upper()
                    fam=next((f for f in family_ids if f in cur),None)
                    if fam and valid_aa(s): aa[fam].append((cur.split()[0], normalize_aa(s), str(p)))
                with open(p,'r',encoding='utf-8',errors='replace') as f:
                    for line in f:
                        if line.startswith('>'):
                            flush(); cur=line[1:].strip(); seq=[]
                        else: seq.append(line.strip())
                    flush()
            except Exception: pass
    # deduplicate exact seq per family
    for fam in family_ids:
        seen=set(); uniq=[]
        for rec in aa[fam]:
            if rec[1] not in seen: seen.add(rec[1]); uniq.append(rec)
        aa[fam]=uniq
        seen=set(); uniq=[]
        for rec in rep[fam]:
            if rec[1] not in seen: seen.add(rec[1]); uniq.append(rec)
        rep[fam]=uniq
    return aa,rep,text,provenance

def parse_hmm_tbl(path):
    hits=[]
    with open(path,'r',encoding='utf-8',errors='replace') as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            parts=line.split()
            if len(parts) < 6: continue
            hits.append({"target":parts[0],"query":parts[2],"evalue":parts[4],"score":parts[5]})
    return hits

def optional_hmm_scan(hmm_dir, aa_by_family, outdir):
    results=[]; status={"requested":bool(hmm_dir),"hmmsearch_found":bool(shutil.which('hmmsearch')),"hmm_files_n":0,"queries_n":0,"hits_n":0}
    if not hmm_dir: return results,status
    hmm_dir=Path(hmm_dir)
    hmms=sorted(hmm_dir.rglob('*.hmm')) if hmm_dir.exists() else []
    status['hmm_files_n']=len(hmms)
    if not hmms or not shutil.which('hmmsearch'): return results,status
    qdir=Path(outdir)/'hmm_queries'; qdir.mkdir(parents=True,exist_ok=True)
    rawdir=Path(outdir)/'hmm_raw'; rawdir.mkdir(parents=True,exist_ok=True)
    for fam,recs in aa_by_family.items():
        if not recs: continue
        q=qdir/f'{fam}.faa'
        with open(q,'w') as f:
            for i,(rid,seq,src) in enumerate(recs,1): f.write(f'>{fam}__Q{i}__{rid}\n{seq}\n')
        status['queries_n'] += len(recs)
        for hmm in hmms:
            tbl=rawdir/f'{fam}__{hmm.stem}.tbl'
            cmd=['hmmsearch','--noali','--tblout',str(tbl),'-E','0.001',str(hmm),str(q)]
            cp=subprocess.run(cmd,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,text=True)
            if cp.returncode != 0:
                results.append({"family_id":fam,"hmm_file":str(hmm),"status":"ERROR","target":"","query":"","evalue":"","score":"","stderr":cp.stderr[:500]})
                continue
            for h in parse_hmm_tbl(tbl):
                results.append({"family_id":fam,"hmm_file":str(hmm),"status":"HIT",**h,"stderr":""})
    status['hits_n']=sum(1 for r in results if r.get('status')=='HIT')
    return results,status

def art95_audit(manifest, public_ref, outdir):
    expected={f'ART_{i:02d}' for i in range(1,96)}
    public_rows=read_tsv(public_ref)
    status={
      "author_reported_family_member_n":95,
      "author_reported_detectable_array_n":28,
      "public_reference_rows_in_release":len(public_rows),
      "exact_manifest_supplied":bool(manifest),
      "exact_manifest_ready":False,
      "exact_manifest_unique_art_id_n":0,
      "missing_art_ids":sorted(expected),
      "rt_aa_available_n":0,
      "yxdd_motif_n":0,
      "claim_boundary":"No 95/95 validation claim is allowed unless a complete exact 95-member catalogue is supplied and passes identity checks. Detectable-array sensitivity should be stratified; 95/95 dual-arm detection is not an expected endpoint because the author reports arrays for 28/95."
    }
    detail=[]
    if manifest:
        rows=read_tsv(manifest)
        ids=[r.get('art_id','').strip() for r in rows if r.get('art_id','').strip()]
        got=set(ids)
        status['exact_manifest_unique_art_id_n']=len(got)
        status['missing_art_ids']=sorted(expected-got)
        status['duplicate_art_id_n']=len(ids)-len(got)
        seqs=[normalize_aa(r.get('rt_aa','')) for r in rows if valid_aa(r.get('rt_aa',''))]
        status['rt_aa_available_n']=len(seqs)
        status['yxdd_motif_n']=sum(1 for s in seqs if re.search(r'Y.DD',s))
        status['exact_manifest_ready']=(got==expected and len(ids)==95 and len(seqs)==95)
        for r in rows:
            s=normalize_aa(r.get('rt_aa','')) if valid_aa(r.get('rt_aa','')) else ''
            detail.append({"art_id":r.get('art_id',''),"accession":r.get('accession',''),"aa_len":len(s) if s else '',"yxdd":bool(re.search(r'Y.DD',s)) if s else '',"has_coordinates":all(str(r.get(c,'')).strip() for c in ['rt_start','rt_end','strand'])})
    write_tsv(Path(outdir)/'STAGE05E3C_ART95_MEMBER_AUDIT.tsv', detail, ['art_id','accession','aa_len','yxdd','has_coordinates'])
    with open(Path(outdir)/'STAGE05E3C_ART95_AUDIT_STATUS.json','w') as f: json.dump(status,f,indent=2,sort_keys=True)
    return status

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--e3b',required=True); ap.add_argument('--e3a',required=True); ap.add_argument('--project-root',required=True)
    ap.add_argument('--outdir',default='outputs'); ap.add_argument('--workdir',default='work'); ap.add_argument('--expected-families',type=int,default=12)
    ap.add_argument('--terms',default=str(Path(__file__).resolve().parents[1]/'refs/KNOWN_SYSTEM_TERMS.tsv'))
    ap.add_argument('--public-art95',default=str(Path(__file__).resolve().parents[1]/'refs/ART95_PUBLIC_REFERENCE.tsv'))
    ap.add_argument('--art95-manifest',default=''); ap.add_argument('--hmm-dir',default='')
    args=ap.parse_args()
    out=Path(args.outdir); work=Path(args.workdir)
    if out.exists(): shutil.rmtree(out)
    if work.exists(): shutil.rmtree(work)
    out.mkdir(parents=True); work.mkdir(parents=True)
    e3b=materialize_input(args.e3b, work/'e3b'); e3a=materialize_input(args.e3a, work/'e3a')
    summary_path=find_basename(e3b,'STAGE05E3B_SUMMARY.json'); queue_path=find_basename(e3b,'STAGE05E3B_REVIEW_QUEUE.tsv')
    e3b_summary=json.load(open(summary_path)); qrows=read_tsv(queue_path)
    families=[r.get('family_id','').strip() for r in qrows if r.get('family_id','').strip()]
    if len(families)!=args.expected_families or len(set(families))!=args.expected_families:
        raise RuntimeError(f"family count mismatch: rows={len(families)} unique={len(set(families))} expected={args.expected_families}")
    if intval(e3b_summary.get('n_families'),-1)!=args.expected_families:
        raise RuntimeError('E3B summary n_families mismatch')
    inv=inventory(e3b)+inventory(e3a)
    for r in inv: r['source_group']='e3b_or_e3a'
    write_tsv(out/'STAGE05E3C_INPUT_INVENTORY.tsv',inv)
    aa,rep,text,prov=collect_family_linked_records([e3b,e3a],families)
    terms=load_terms(args.terms)
    evidence=[]; closure=[]
    for r in qrows:
        fam=r['family_id']
        text_hits=scan_text('\n'.join(text.get(fam,[])),terms)
        direct_known=intval(r.get('known_system_keyword_candidate_n',0))+intval(r.get('any_known_retro_system_keyword_candidate_n_10kb',0))+intval(r.get('strict_rt_keyword_candidate_n_10kb',0))
        cats=sorted(set(h[1] for h in text_hits))
        seqs=[x[1] for x in aa.get(fam,[])]; reps=[x[1] for x in rep.get(fam,[])]
        aa_lens=[len(s) for s in seqs]; rep_lens=[len(s) for s in reps]
        official=intval(r.get('official_supported_witness_n', r.get('official_same_strand_high_overlap_or_near_exact_witness_n',0)))
        witness=intval(r.get('signature_witness_candidate_n',0))
        stratum=r.get('stage05e3b_review_stratum','')
        if direct_known>0 or cats or 'CONFLICT' in stratum:
            tier='E3C5_KNOWN_SYSTEM_OR_SPECIAL_CONFLICT_REVIEW'
        elif official==0:
            tier='E3C4_UNRESOLVED_CODING_SUPPORT'
        elif stratum.startswith('E3B1') or stratum.startswith('E3B2'):
            tier='E3C1_PRIORITY_BIOLOGICAL_CLOSURE'
        elif stratum.startswith('E3B4'):
            tier='E3C2_NARROW_REPLICATED_CLOSURE'
        else:
            tier='E3C3_HETEROGENEOUS_RETAINED'
        closure.append({
          **r,
          'e3c_review_tier':tier,
          'family_linked_unique_aa_n':len(seqs),
          'family_linked_aa_median_len':median_or_blank(aa_lens),
          'family_linked_aa_pairwise_median_3mer_jaccard':pairwise_median_jaccard(seqs,3) if len(seqs)>=2 else '',
          'family_linked_repeat_sequence_n':len(reps),
          'family_linked_repeat_median_len':median_or_blank(rep_lens),
          'family_linked_repeat_pairwise_median_3mer_jaccard':pairwise_median_jaccard(reps,3) if len(reps)>=2 else '',
          'text_mapped_known_system_categories':';'.join(cats),
          'positive_known_system_evidence_present':bool(direct_known>0 or cats),
          'e3c_no_hit_is_not_novelty_proof':True
        })
        for tid,cat,pat in text_hits:
            evidence.append({'family_id':fam,'evidence_type':'TEXT_TERM','category':cat,'term_id':tid,'pattern':pat,'source':'family-linked E3A/E3B text','interpretation':'positive conflict evidence only'})
        if direct_known>0:
            evidence.append({'family_id':fam,'evidence_type':'E3B_CONTEXT','category':'KNOWN_SYSTEM_CONTEXT','term_id':'E3B','pattern':r.get('context_interpretation',''),'source':'E3B review queue','interpretation':'positive conflict evidence only'})
    hmm_results,hmm_status=optional_hmm_scan(args.hmm_dir,aa,out)
    if hmm_results:
        write_tsv(out/'STAGE05E3C_OPTIONAL_HMM_HITS.tsv',hmm_results)
        hitfams={r['family_id'] for r in hmm_results if r.get('status')=='HIT'}
        for c in closure:
            if c['family_id'] in hitfams:
                c['positive_hmm_conflict_evidence_present']=True
                c['e3c_review_tier']='E3C5_KNOWN_SYSTEM_OR_SPECIAL_CONFLICT_REVIEW'
            else: c['positive_hmm_conflict_evidence_present']=False
    else:
        for c in closure: c['positive_hmm_conflict_evidence_present']=False
        write_tsv(out/'STAGE05E3C_OPTIONAL_HMM_HITS.tsv',[],['family_id','hmm_file','status','target','query','evalue','score','stderr'])
    # AA and repeat audits
    aar=[]; rr=[]
    for fam in families:
        for rid,seq,src in aa[fam]: aar.append({'family_id':fam,'record_id':rid,'aa_len':len(seq),'yxdd':bool(re.search(r'Y.DD',seq)),'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),'source_file':src})
        for rid,seq,src in rep[fam]: rr.append({'family_id':fam,'record_id':rid,'nt_len':len(seq),'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),'source_file':src})
    write_tsv(out/'STAGE05E3C_FAMILY_PROTEIN_SEQUENCE_AUDIT.tsv',aar,['family_id','record_id','aa_len','yxdd','sequence_sha256','source_file'])
    write_tsv(out/'STAGE05E3C_REPEAT_SEQUENCE_AUDIT.tsv',rr,['family_id','record_id','nt_len','sequence_sha256','source_file'])
    write_tsv(out/'STAGE05E3C_KNOWN_SYSTEM_EVIDENCE.tsv',evidence,['family_id','evidence_type','category','term_id','pattern','source','interpretation'])
    # Stable priority order: tier, previous review order, family_id. This is a review queue, not a DNA rerank.
    tier_order={'E3C1_PRIORITY_BIOLOGICAL_CLOSURE':1,'E3C2_NARROW_REPLICATED_CLOSURE':2,'E3C3_HETEROGENEOUS_RETAINED':3,'E3C4_UNRESOLVED_CODING_SUPPORT':4,'E3C5_KNOWN_SYSTEM_OR_SPECIAL_CONFLICT_REVIEW':5}
    closure.sort(key=lambda x:(tier_order.get(x['e3c_review_tier'],9),intval(x.get('review_order',x.get('stage05e3b_review_order_group',9)),9),x['family_id']))
    for i,c in enumerate(closure,1): c['e3c_review_order']=i
    write_tsv(out/'STAGE05E3C_FAMILY_CLOSURE.tsv',closure)
    write_tsv(out/'STAGE05E3C_REVIEW_QUEUE.tsv',closure)
    # Family cards
    cards=out/'family_cards'; cards.mkdir()
    for c in closure:
        fam=c['family_id']
        lines=[f"# {fam} - Stage05E3C family card",'',f"E3C review tier: {c['e3c_review_tier']}",f"E3B stratum: {c.get('stage05e3b_review_stratum','')}",f"Candidates/genomes: {c.get('candidate_n','')} / {c.get('genome_n','')}",f"Lineage breadth: {c.get('lineage_breadth_class','')}",f"Official-supported witnesses: {c.get('official_supported_witness_n','')}",f"Median official AA: {c.get('official_witness_median_aa','')}",f"Topology: {c.get('modal_transcriptional_repeat_side_E2','')} (fraction {c.get('modal_side_fraction_E2','')})",f"Median repeat-ORF gap: {c.get('median_repeat_orf_gap_bp_E2','')}",f"Known-system context: {c.get('context_interpretation','')}",f"Family-linked AA sequences recovered in E3C input: {c.get('family_linked_unique_aa_n','0')}",f"Family-linked repeat sequences recovered in E3C input: {c.get('family_linked_repeat_sequence_n','0')}",'',"Claim boundary: retained architecture hypothesis; not ART identity or novelty proof."]
        (cards/f'{fam}.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    art95=art95_audit(args.art95_manifest or None,args.public_art95,out)
    counts={}
    for c in closure: counts[c['e3c_review_tier']]=counts.get(c['e3c_review_tier'],0)+1
    summary={
      'stage':STAGE,'version':VERSION,'n_families':len(closure),'n_candidates':intval(e3b_summary.get('n_candidates',0)),'n_sequence_ids':intval(e3b_summary.get('n_sequence_ids',0)),
      'e3b_gate':e3b_summary.get('gate',''),'e3b_gate_pass':bool(e3b_summary.get('gate_pass',False)),'review_tier_counts':counts,
      'n_families_with_positive_known_system_evidence':sum(1 for c in closure if c['positive_known_system_evidence_present'] or c['positive_hmm_conflict_evidence_present']),
      'n_families_with_family_linked_aa_sequences':sum(1 for c in closure if intval(c['family_linked_unique_aa_n'])>0),
      'n_families_with_family_linked_repeat_sequences':sum(1 for c in closure if intval(c['family_linked_repeat_sequence_n'])>0),
      'optional_hmm_status':hmm_status,'art95_exact_manifest_ready':art95['exact_manifest_ready'],
      'candidate_filtering_performed':False,'family_membership_changed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,'selection_uses_structure_results':False,
      'gate':'PASS_E3C_FAMILY_CLOSURE_AND_KNOWN_SYSTEM_EVIDENCE_AUDIT_FROZEN','gate_pass':True,
      'claim_boundary':'E3C freezes family-level closure/review evidence without changing E2A/E3B membership or DNA scores. Positive known-system evidence triggers conflict review. Negative evidence is not novelty proof; ART identity/function remain unresolved.'
    }
    with open(out/'STAGE05E3C_SUMMARY.json','w') as f: json.dump(summary,f,indent=2,sort_keys=True)
    with open(out/'STAGE05E3C_GATE.json','w') as f: json.dump({'gate':summary['gate'],'gate_pass':True,'expected_family_n':args.expected_families,'observed_family_n':len(closure),'family_membership_changed':False,'dna_rank_recomputed':False,'art95_exact_manifest_ready':art95['exact_manifest_ready'],'claim_boundary':summary['claim_boundary']},f,indent=2,sort_keys=True)
    print(json.dumps(summary,indent=2,sort_keys=True))

if __name__=='__main__': main()
