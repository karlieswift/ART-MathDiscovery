#!/usr/bin/env python3
import argparse, csv, hashlib, json, os, re, sys, time, urllib.parse, urllib.request
from pathlib import Path

try:
    from Bio import SeqIO
except Exception as e:
    raise SystemExit("[FAIL] Biopython is required in this stage: " + repr(e))

CAND_SHA = "73007864fdef7dba662c9b90a02e2b915359c1fbd93a6061edeea35afde07349"
MANIFEST_SHA = "24f69c4f463678413d772cc567332bd0b94ff18019d53ea7e454aa5e30fb5ea4"
BOARD_SHA = "272b59fd09d853bc7336e364a64e909c3cca7b082a684c69355d873707be9337"
EXPECTED_REGIONS = 383
EXPECTED_SEQS = 174
NB = 10000
RT_NEAR = 2000

RT_PATTERNS = [
    re.compile(r"reverse\s+transcriptase", re.I),
    re.compile(r"rna[- ]directed\s+dna\s+polymerase", re.I),
    re.compile(r"\bretron\b", re.I),
    re.compile(r"\bretroelement\b", re.I),
    re.compile(r"group\s*ii\s+intron.*(?:reverse\s+transcriptase|maturase)", re.I),
    re.compile(r"(?:reverse\s+transcriptase|maturase).*group\s*ii\s+intron", re.I),
]

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()

def sha256_seq(seq):
    s=''.join(str(seq).split()).upper().encode()
    return hashlib.sha256(s).hexdigest()

def read_tsv(p):
    with open(p, newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(p, rows, fields):
    with open(p,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=fields,extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

def feature_text(feat):
    vals=[]
    for k,v in feat.qualifiers.items():
        vals.extend([str(x) for x in v])
    return ' | '.join(vals)

def strict_rt(text):
    return any(p.search(text or '') for p in RT_PATTERNS)

def interval_distance(a0,a1,b0,b1):
    if a1 <= b0: return b0-a1
    if b1 <= a0: return a0-b1
    return 0

def wrap_extract(seq, start, end, circular):
    n=len(seq)
    if not circular:
        s=max(0,start); e=min(n,end)
        return str(seq[s:e]), s, e, int(start<0), int(end>n), 0
    width=end-start
    if width >= n:
        return str(seq), 0, n, 0, 0, 1
    s=start % n; e=end % n
    if s < e:
        x=str(seq[s:e])
    else:
        x=str(seq[s:]) + str(seq[:e])
    return x, s, e, 0, 0, 1

def _looks_like_genbank(data):
    if not data or len(data) < 1000:
        return False
    head=data[:10000].lstrip()
    return head.startswith(b"LOCUS") and b"ORIGIN" in data and data.rstrip().endswith(b"//")

def _payload_preview(data, n=240):
    if data is None:
        return ""
    x=data[:n].decode("utf-8", "replace")
    return re.sub(r"\s+", " ", x).strip()

def fetch_gbff(acc, out, api_key=None, retries=5):
    # Existing cache is trusted only if it is syntactically GenBank-like.
    if out.exists() and out.stat().st_size > 1000:
        data=out.read_bytes()
        if _looks_like_genbank(data):
            return
        out.unlink()

    tool=os.environ.get("NCBI_TOOL", "ART_MathDiscovery_Stage04B_R1")
    email=os.environ.get("NCBI_EMAIL", "")
    common={"db":"nuccore","id":acc,"retmode":"text"}
    if api_key:
        common["api_key"]=api_key
    if tool:
        common["tool"]=tool
    if email:
        common["email"]=email

    endpoint_specs=[]
    p1=dict(common); p1["rettype"]="gbwithparts"
    endpoint_specs.append(("eutils_gbwithparts", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"+urllib.parse.urlencode(p1)))
    p2=dict(common); p2["rettype"]="gb"
    endpoint_specs.append(("eutils_gb", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"+urllib.parse.urlencode(p2)))
    p3={"id":acc,"db":"nuccore","report":"genbank","retmode":"text"}
    endpoint_specs.append(("sviewer_genbank", "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?"+urllib.parse.urlencode(p3)))

    delay=0.15 if api_key else 0.45
    errors=[]
    for endpoint_name,url in endpoint_specs:
        last=None
        for i in range(retries):
            try:
                req=urllib.request.Request(url,headers={
                    "User-Agent":"ART-MathDiscovery-Stage04B-R1/1.0.2",
                    "Accept":"text/plain,*/*;q=0.8",
                    "Connection":"close",
                })
                with urllib.request.urlopen(req,timeout=120) as r:
                    data=r.read()
                if not _looks_like_genbank(data):
                    raise RuntimeError(
                        f"unexpected payload endpoint={endpoint_name} bytes={len(data)} preview={_payload_preview(data)!r}"
                    )
                out.write_bytes(data)
                time.sleep(delay)
                return
            except Exception as e:
                last=e
                time.sleep(min(30, 2**i))
        errors.append(f"{endpoint_name}:{last}")
    raise RuntimeError(f"fetch failed {acc}: " + " || ".join(errors))

def locate_fasta(stage04a_root, relpath, seq_id):
    cands=[stage04a_root/relpath, stage04a_root/'inputs'/f'{seq_id}.frozen.fasta']
    for p in cands:
        if p.exists(): return p
    raise FileNotFoundError(f"frozen FASTA not found for {seq_id}; checked: "+', '.join(map(str,cands)))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('stage04a_root')
    ap.add_argument('--skip-download',action='store_true')
    args=ap.parse_args()
    pkg=Path(__file__).resolve().parents[1]
    out=pkg/'outputs'; out.mkdir(exist_ok=True)
    work=pkg/'work'/'gbff'; work.mkdir(parents=True,exist_ok=True)
    stage=Path(args.stage04a_root).resolve()
    candp=pkg/'inputs'/'STAGE04A_R1_DUAL_ARM_REGIONS.tsv'
    manp=stage/'outputs'/'STAGE04A_UNIVERSE_MANIFEST.tsv'
    boardp=stage/'outputs'/'STAGE04A_ALL_FROZEN_UNION_REGIONS.tsv'
    errors=[]
    if sha256_file(candp)!=CAND_SHA: errors.append('candidate_registry_sha256_mismatch')
    if not manp.exists() or sha256_file(manp)!=MANIFEST_SHA: errors.append('universe_manifest_sha256_mismatch_or_missing')
    if not boardp.exists() or sha256_file(boardp)!=BOARD_SHA: errors.append('frozen_region_board_sha256_mismatch_or_missing')
    if errors:
        print('[FAIL] canonical source lock:',errors); return 2
    candidates=read_tsv(candp)
    if len(candidates)!=EXPECTED_REGIONS:
        print('[FAIL] expected 383 candidates, got',len(candidates)); return 2
    seq_ids=sorted(set(r['sequence_id'] for r in candidates))
    if len(seq_ids)!=EXPECTED_SEQS:
        print('[FAIL] expected 174 sequences, got',len(seq_ids)); return 2
    manifest={r['sequence_id']:r for r in read_tsv(manp)}
    provenance=[]; recs={}; frozen_seqs={}
    api_key=os.environ.get('NCBI_API_KEY')
    for i,sid in enumerate(seq_ids,1):
        if sid not in manifest: errors.append(f'{sid}:missing_manifest'); continue
        m=manifest[sid]
        try:
            fasta=locate_fasta(stage, m.get('frozen_fasta',''), sid)
            # Stage04A manifest fasta_sha256 is the SHA256 of the frozen FASTA FILE bytes,
            # not the normalized nucleotide sequence. Preserve that original byte-level lock,
            # then compute a separate normalized sequence digest for GenBank sequence identity.
            expected_file_sha=m['fasta_sha256']
            observed_file_sha=sha256_file(fasta)
            if observed_file_sha!=expected_file_sha:
                errors.append(f'{sid}:frozen_fasta_file_sha_mismatch')
                continue
            frozen_records=list(SeqIO.parse(str(fasta),'fasta'))
            if len(frozen_records)!=1:
                errors.append(f'{sid}:frozen_fasta_record_count_{len(frozen_records)}')
                continue
            frozen_rec=frozen_records[0]
            frozen_seq_sha=sha256_seq(frozen_rec.seq)
            gb=work/f'{sid}.gbff'
            if not args.skip_download: fetch_gbff(sid,gb,api_key)
            if not gb.exists():
                errors.append(f'{sid}:gbff_missing'); continue
            gbrecs=list(SeqIO.parse(str(gb),'genbank'))
            if len(gbrecs)!=1:
                errors.append(f'{sid}:gbff_record_count_{len(gbrecs)}'); continue
            gr=gbrecs[0]
            gsha=sha256_seq(gr.seq)
            match=(gsha==frozen_seq_sha)
            if not match:
                errors.append(f'{sid}:gbff_sequence_sha_mismatch')
            recs[sid]=gr; frozen_seqs[sid]=frozen_rec
            provenance.append({
                'sequence_id':sid,'frozen_fasta':str(fasta),
                'expected_frozen_fasta_file_sha256':expected_file_sha,
                'observed_frozen_fasta_file_sha256':observed_file_sha,
                'frozen_fasta_sequence_sha256':frozen_seq_sha,
                'gbff_path':str(gb),'gbff_file_sha256':sha256_file(gb),
                'gbff_sequence_sha256':gsha,'exact_sequence_match':str(match).lower(),
                'gbff_description':gr.description,'topology':gr.annotations.get('topology',''),
                'length_bp':len(gr.seq)
            })
            print(f'[{i}/{len(seq_ids)}] {sid} exact={match} features={len(gr.features)}')
        except Exception as e:
            errors.append(f'{sid}:{type(e).__name__}:{e}')
    write_tsv(out/'STAGE04B_R1_SEQUENCE_PROVENANCE.tsv', provenance,
              ['sequence_id','frozen_fasta','expected_frozen_fasta_file_sha256','observed_frozen_fasta_file_sha256','frozen_fasta_sequence_sha256','gbff_path','gbff_file_sha256','gbff_sequence_sha256','exact_sequence_match','gbff_description','topology','length_bp'])
    if errors:
        gate={'stage':'Stage04B-R1','gate':'HOLD_STAGE04B_R1_EXACT_SEQUENCE_PROVENANCE_INCOMPLETE','gate_pass':False,
              'n_expected_sequences':EXPECTED_SEQS,'n_provenance_rows':len(provenance),'errors':errors,
              'biological_interpretation_authorized':False,'mathematical_rescoring_authorized':False,'retuning_authorized':False}
        (out/'STAGE04B_R1_GATE.json').write_text(json.dumps(gate,indent=2),encoding='utf-8')
        print(json.dumps(gate,indent=2)); return 4

    feat_rows=[]; cand_rows=[]; bio_rows=[]; rt_rows=[]
    reg_fa=open(out/'STAGE04B_R1_CANDIDATE_REGIONS.fasta','w')
    nb_fa=open(out/'STAGE04B_R1_NEIGHBORHOODS.fasta','w')
    prot_fa=open(out/'STAGE04B_R1_NEIGHBORHOOD_PROTEINS.faa','w')
    for c in candidates:
        sid=c['sequence_id']; gr=recs[sid]; seq=frozen_seqs[sid].seq
        s=int(c['start0']); e=int(c['end0']); n=len(seq)
        topology=(gr.annotations.get('topology','') or '').lower(); circular=(topology=='circular')
        regseq=str(seq[s:e])
        nh_start=s-NB; nh_end=e+NB
        nhseq, xs, xe, lt, rt, wrapped=wrap_extract(seq,nh_start,nh_end,circular)
        key=f"{sid}|{c['region_id']}|rank={c['region_rank']}|{s}-{e}"
        reg_fa.write(f">{key}\n{regseq}\n")
        nb_fa.write(f">{key}|neighborhood=+/-{NB}|circular={int(circular)}\n{nhseq}\n")
        related=[]
        for fi,feat in enumerate(gr.features):
            if feat.type=='source': continue
            try:
                fs=int(feat.location.start); fe=int(feat.location.end); strand=feat.location.strand or 0
            except Exception:
                continue
            dist=interval_distance(s,e,fs,fe)
            if dist>NB: continue
            txt=feature_text(feat); isrt=strict_rt(txt)
            q=feat.qualifiers
            row={
                'sequence_id':sid,'region_id':c['region_id'],'region_rank':c['region_rank'],
                'candidate_start0':s,'candidate_end0':e,'feature_index':fi,'feature_type':feat.type,
                'feature_start0':fs,'feature_end0':fe,'strand':strand,'distance_to_candidate_bp':dist,
                'overlaps_candidate':str(dist==0).lower(),'rt_strict_keyword':str(isrt).lower(),
                'gene':';'.join(q.get('gene',[])),'locus_tag':';'.join(q.get('locus_tag',[])),
                'product':';'.join(q.get('product',[])),'protein_id':';'.join(q.get('protein_id',[])),
                'note':';'.join(q.get('note',[])),'db_xref':';'.join(q.get('db_xref',[]))
            }
            feat_rows.append(row); related.append(row)
            if feat.type=='CDS' and q.get('translation'):
                aa=''.join(q['translation'][0].split())
                pid=(q.get('protein_id') or q.get('locus_tag') or [f'{sid}_feature{fi}'])[0]
                prot_fa.write(f">{sid}|{c['region_id']}|{pid}|dist={dist}|rtkw={int(isrt)}\n{aa}\n")
        rt_feats=[r for r in related if r['feature_type']=='CDS' and r['rt_strict_keyword']=='true']
        nearest=min([int(r['distance_to_candidate_bp']) for r in rt_feats], default=None)
        if nearest is None: triage='NO_RT_KEYWORD_WITHIN_10KB'
        elif nearest<=RT_NEAR: triage='RT_KEYWORD_OVERLAP_OR_WITHIN_2KB'
        else: triage='RT_KEYWORD_WITHIN_10KB'
        cds=sum(r['feature_type']=='CDS' for r in related)
        rna=sum(r['feature_type'] in {'ncRNA','misc_RNA','tRNA','rRNA','tmRNA'} for r in related)
        reps=sum(r['feature_type']=='repeat_region' for r in related)
        cr={**c,'candidate_sequence_len':len(regseq),'neighborhood_radius_bp':NB,
            'neighborhood_extract_start':xs,'neighborhood_extract_end':xe,'neighborhood_wrapped':wrapped,
            'left_truncated':lt,'right_truncated':rt,'topology':topology}
        cand_rows.append(cr)
        br={
            'sequence_id':sid,'region_id':c['region_id'],'region_rank':c['region_rank'],'start0':s,'end0':e,
            'best_q_acnv_family':c['best_q_acnv_family'],'best_q_recurrence':c['best_q_recurrence'],
            'recurrence_best_k':c['recurrence_best_k'],'recurrence_best_lag_bp':c['recurrence_best_lag_bp'],
            'n_features_within_10kb':len(related),'n_cds_within_10kb':cds,'n_rna_features_within_10kb':rna,
            'n_repeat_region_features_within_10kb':reps,'n_rt_keyword_cds_within_10kb':len(rt_feats),
            'nearest_rt_keyword_distance_bp':'' if nearest is None else nearest,'rt_keyword_triage_class':triage,
            'annotation_only_no_candidate_filtering':'true'
        }
        bio_rows.append(br)
        if rt_feats:
            for r in rt_feats:
                rt_rows.append({**br,**{k:r[k] for k in ['feature_start0','feature_end0','strand','distance_to_candidate_bp','gene','locus_tag','product','protein_id','note']}})
    reg_fa.close(); nb_fa.close(); prot_fa.close()
    write_tsv(out/'STAGE04B_R1_CANDIDATE_NEIGHBORHOODS.tsv',cand_rows,list(cand_rows[0].keys()))
    ff=['sequence_id','region_id','region_rank','candidate_start0','candidate_end0','feature_index','feature_type','feature_start0','feature_end0','strand','distance_to_candidate_bp','overlaps_candidate','rt_strict_keyword','gene','locus_tag','product','protein_id','note','db_xref']
    write_tsv(out/'STAGE04B_R1_NEIGHBORHOOD_FEATURES.tsv',feat_rows,ff)
    bf=['sequence_id','region_id','region_rank','start0','end0','best_q_acnv_family','best_q_recurrence','recurrence_best_k','recurrence_best_lag_bp','n_features_within_10kb','n_cds_within_10kb','n_rna_features_within_10kb','n_repeat_region_features_within_10kb','n_rt_keyword_cds_within_10kb','nearest_rt_keyword_distance_bp','rt_keyword_triage_class','annotation_only_no_candidate_filtering']
    write_tsv(out/'STAGE04B_R1_CANDIDATE_BIOLOGY_SUMMARY.tsv',bio_rows,bf)
    rtf=bf+['feature_start0','feature_end0','strand','distance_to_candidate_bp','gene','locus_tag','product','protein_id','note']
    write_tsv(out/'STAGE04B_R1_RT_KEYWORD_TRIAGE.tsv',rt_rows,rtf)
    from collections import Counter
    cc=Counter(r['rt_keyword_triage_class'] for r in bio_rows)
    summary={
        'stage':'Stage04B-R1','n_candidate_regions':len(candidates),'n_candidate_sequences':len(seq_ids),
        'n_exact_sequence_provenance':len(provenance),'fixed_neighborhood_radius_bp':NB,'rt_near_radius_bp':RT_NEAR,
        'rt_keyword_triage_counts':dict(sorted(cc.items())),
        'n_candidates_with_any_rt_keyword_within_10kb':sum(r['n_rt_keyword_cds_within_10kb']>0 for r in bio_rows),
        'n_candidates_with_rt_keyword_overlap_or_within_2kb':sum(r['rt_keyword_triage_class']=='RT_KEYWORD_OVERLAP_OR_WITHIN_2KB' for r in bio_rows),
        'candidate_filtering_performed':False,'score_recomputed':False,'region_ranking_recomputed':False,
        'threshold_changed':False,'retuning_authorized':False,
        'interpretation':'Post-freeze biological annotation and descriptive RT-keyword triage. No-keyword candidates remain active because sequence-based RT/domain screening is required next.'
    }
    gate={
        'stage':'Stage04B-R1','gate':'PASS_STAGE04B_R1_ALL_383_FROZEN_CANDIDATES_ANNOTATED_WITH_EXACT_SEQUENCE_PROVENANCE','gate_pass':True,
        'n_candidate_regions':len(candidates),'n_candidate_sequences':len(seq_ids),'n_exact_sequence_provenance':len(provenance),
        'candidate_registry_sha256':sha256_file(candp),'stage04a_manifest_sha256':sha256_file(manp),'stage04a_board_sha256':sha256_file(boardp),
        'candidate_filtering_performed':False,'score_recomputed':False,'region_ranking_recomputed':False,
        'threshold_changed':False,'truth_used_for_mathematical_scoring':False,'retuning_authorized':False,
        'biological_annotation_authorized':True,'new_ART_identity_calls_authorized':False
    }
    (out/'STAGE04B_R1_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    (out/'STAGE04B_R1_GATE.json').write_text(json.dumps(gate,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2)); print(json.dumps(gate,indent=2)); return 0

if __name__=='__main__':
    sys.exit(main())
