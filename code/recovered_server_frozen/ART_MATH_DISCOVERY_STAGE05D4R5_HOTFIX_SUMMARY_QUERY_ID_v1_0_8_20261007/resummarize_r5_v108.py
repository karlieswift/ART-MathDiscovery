#!/usr/bin/env python3
import argparse,csv,gzip,json,re
from pathlib import Path

RT_PAT=re.compile(r'\b(reverse transcriptase|retrotranscriptase|rna[- ]dependent dna polymerase|rvt|group ii intron.*reverse|retr(on|on-associated).*reverse|diversity[- ]generating retroelement)\b',re.I)

def sp_desc(p):
    d={}; op=gzip.open if str(p).endswith('.gz') else open
    with op(p,'rt',errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                h=line[1:].strip()
                parts=h.split('|')
                acc=parts[1] if '|' in h and len(parts)>2 else h.split()[0]
                d[acc]=h
    return d

def target_acc(t):
    z=t.split()[0]
    if '|' in z:
        p=z.split('|')
        if len(p)>1: return p[1]
    m=re.search(r'(?:AF-)?([A-Z0-9]+)-F\d+',z)
    return m.group(1) if m else z

def map_query(raw, frozen_queries):
    """Map Foldseek structure query header back to frozen AA_<sha> query id.
    Foldseek createdb/easy-search can emit a structure-file-derived id such as
    AA_<sha>_model_0 while Stage05D4-R5 freeze uses AA_<sha>.
    We require a unique mapping to one of the 43 frozen query ids.
    """
    s=raw.strip()
    # Try exact and basename/stem variants first.
    variants=[]
    def add(v):
        if v and v not in variants: variants.append(v)
    add(s)
    add(Path(s).name)
    for v in list(variants):
        if v.endswith('.cif'): add(v[:-4])
        if v.endswith('.mmcif'): add(v[:-6])
    for v in list(variants):
        if v.endswith('_model_0'): add(v[:-8])
    for v in variants:
        if v in frozen_queries: return v, 'exact_or_suffix_strip'
    # Final robust rule: the emitted id may append model/chain text after the
    # frozen AA_<sha> id. AA_<sha> ids are unique and fixed length/unambiguous.
    matches=[q for q in frozen_queries if s.startswith(q) or Path(s).name.startswith(q)]
    if len(matches)==1: return matches[0], 'unique_frozen_prefix'
    raise RuntimeError(f'Cannot uniquely map Foldseek query id {raw!r}; matches={matches}')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--r5',required=True)
    ap.add_argument('--r4',required=True)
    a=ap.parse_args()
    r5=Path(a.r5); r4=Path(a.r4)
    freeze=r5/'frozen'/'STAGE05D4R5_P1_FREEZE.tsv'
    boltz=r5/'boltz_out'
    m8=r5/'raw'/'P1_FULL3D_AFSP.m8'
    sprot=r4/'references'/'uniprot_sprot.fasta.gz'
    out=r5/'outputs'; out.mkdir(exist_ok=True)
    if not freeze.exists(): raise SystemExit(f'[FAIL] missing {freeze}')
    if not m8.exists(): raise SystemExit(f'[FAIL] missing {m8}')
    if not sprot.exists(): raise SystemExit(f'[FAIL] missing {sprot}')

    fr=list(csv.DictReader(open(freeze),delimiter='\t'))
    assert len(fr)==43, len(fr)
    frozen_queries=[x['query'] for x in fr]; fq=set(frozen_queries)

    conf={}
    for p in boltz.glob('predictions/*/confidence_*_model_0.json'):
        q=p.parent.name
        conf[q]=json.load(open(p))
    if len(conf)!=43:
        raise SystemExit(f'[FAIL] expected 43 confidence JSONs, got {len(conf)}')

    desc=sp_desc(sprot)
    hits={}
    qmap_rows=[]
    raw_rows=0
    with open(m8) as f:
        for line_no,line in enumerate(f,1):
            if not line.strip(): continue
            raw_rows += 1
            z=line.rstrip('\n').split('\t')
            if len(z)<11:
                raise RuntimeError(f'm8 line {line_no} has {len(z)} columns, expected >=11: {line.rstrip()}')
            rawq,t=z[0],z[1]
            q,method=map_query(rawq,fq)
            qmap_rows.append({'line_no':line_no,'raw_foldseek_query':rawq,'frozen_query':q,'mapping_method':method})
            a0=target_acc(t); ds=desc.get(a0,'')
            rec={
                'target':t,'target_accession':a0,'evalue':float(z[2]),'bits':float(z[3]),
                'qcov':float(z[4]),'tcov':float(z[5]),'alnlen':int(float(z[6])),
                'qlen':int(float(z[7])),'tlen':int(float(z[8])),'fident':float(z[9]),
                'alntmscore':float(z[10]),'description':ds,'rt_text':bool(RT_PAT.search(ds)),
                'raw_foldseek_query':rawq
            }
            hits.setdefault(q,[]).append(rec)
    # Stable ranking for best hit: Foldseek normally returns ranked hits, but make
    # the postprocessor deterministic independently of row order.
    for q in hits:
        hits[q].sort(key=lambda h:(h['evalue'],-h['bits'],-h['alntmscore']))

    # Audit mapping table.
    map_path=out/'STAGE05D4R5_M8_QUERY_ID_MAP.tsv'
    with open(map_path,'w',newline='') as f:
        fields=['line_no','raw_foldseek_query','frozen_query','mapping_method']
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(qmap_rows)

    rows=[]
    for x in fr:
        q=x['query']; hs=hits.get(q,[]); rhs=[h for h in hs if h['rt_text']]; c=conf.get(q,{})
        best=hs[0] if hs else {}; br=rhs[0] if rhs else {}
        rows.append({
            'query':q,'aa_length':x['aa_length'],'candidate_key_n':x['candidate_key_n'],
            'candidate_keys':x['candidate_keys'],'strata':x['strata'],
            'min_orf_to_candidate_distance_bp':x['min_orf_to_candidate_distance_bp'],
            'confidence_score':c.get('confidence_score',''),'ptm':c.get('ptm',''),
            'complex_plddt':c.get('complex_plddt',''),'n_full3d_hits':len(hs),
            'best_target':best.get('target',''),'best_target_accession':best.get('target_accession',''),
            'best_evalue':best.get('evalue',''),'best_alntmscore':best.get('alntmscore',''),
            'best_qcov':best.get('qcov',''),'best_tcov':best.get('tcov',''),
            'best_description':best.get('description',''),'n_rt_text_full3d_hits':len(rhs),
            'best_rt_target':br.get('target',''),'best_rt_accession':br.get('target_accession',''),
            'best_rt_evalue':br.get('evalue',''),'best_rt_alntmscore':br.get('alntmscore',''),
            'best_rt_qcov':br.get('qcov',''),'best_rt_tcov':br.get('tcov',''),
            'best_rt_description':br.get('description',''),
            'triage_label':'FULL3D_RT_TEXT_HIT_REQUIRES_VALIDATION' if rhs else ('FULL3D_CURATED_NON_RT_TRIAGE' if hs else 'NO_FULL3D_AFSP_HIT_NOT_NEGATIVE')
        })
    triage=out/'STAGE05D4R5_P1_FULL3D_TRIAGE.tsv'
    fields=list(rows[0])
    with open(triage,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(rows)
    cnt={}
    for r in rows: cnt[r['triage_label']]=cnt.get(r['triage_label'],0)+1
    n_hit=sum(1 for q in frozen_queries if hits.get(q))
    n_rt=sum(1 for q in frozen_queries if any(h['rt_text'] for h in hits.get(q,[])))
    unique_raw=len({r['raw_foldseek_query'] for r in qmap_rows})
    unique_mapped=len({r['frozen_query'] for r in qmap_rows})
    summary={
        'stage':'Stage05D4-R5','n_P1':43,'n_boltz_predictions':len(conf),
        'n_foldseek_m8_rows':raw_rows,'n_foldseek_raw_query_ids':unique_raw,
        'n_foldseek_mapped_frozen_queries':unique_mapped,
        'n_queries_with_full3d_hit':n_hit,'n_queries_with_rt_text_full3d_hit':n_rt,
        'triage_counts':cnt,'candidate_filtering_performed':False,
        'foldseek_tm_field':'alntmscore',
        'query_id_normalization':'Map Foldseek structure-file-derived query IDs (e.g. AA_<sha>_model_0) uniquely back to the frozen AA_<sha> IDs before summarization.',
        'claim_boundary':'Boltz-2 full-coordinate predictions plus Foldseek are post-discovery protein interpretation. RT-text structural hits require catalytic-core/topology, DNA architecture, known-family dereplication and ultimately functional validation before ART claims.'
    }
    json.dump(summary,open(out/'STAGE05D4R5_SUMMARY.json','w'),indent=2)
    gate={**summary,
        'gate':'PASS_P1_BOLTZ2_FULL3D_AND_AFSP_FOLDSEEK_QUERY_IDS_NORMALIZED',
        'gate_pass':len(conf)==43,
        'dna_score_recomputed':False,'dna_rank_recomputed':False,
        'p1_compute_priority_only':True,'single_sequence_fallback_used':False,
        'implementation_hotfix':'v1.0.8 postprocessing-only: normalize Foldseek query IDs derived from *_model_0.cif back to the frozen AA_<sha> IDs; reuse existing v1.0.7 m8 and 43 Boltz predictions; no Foldseek rerun, no prediction rerun, no threshold change.',
        'scientific_settings_changed':False,'prediction_workers_rerun':False,'foldseek_rerun':False
    }
    json.dump(gate,open(out/'GATE_STAGE05D4R5.json','w'),indent=2)
    print('[AUDIT] raw m8 rows=',raw_rows)
    print('[AUDIT] unique raw Foldseek query IDs=',unique_raw)
    print('[AUDIT] unique mapped frozen queries=',unique_mapped)
    print(json.dumps(summary,indent=2))
    print(json.dumps(gate,indent=2))

if __name__=='__main__': main()
