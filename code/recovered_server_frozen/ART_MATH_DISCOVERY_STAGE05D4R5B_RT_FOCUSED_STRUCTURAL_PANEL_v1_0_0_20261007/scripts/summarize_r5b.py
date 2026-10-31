#!/usr/bin/env python3
import argparse,csv,json,os,re,collections

def normq(q):
    return re.sub(r'_model_0$','',q)
def normacc(t):
    m=re.search(r'AF-([A-Za-z0-9]+)-F\d+-model',t)
    if m:return m.group(1)
    if '|' in t:
        p=t.split('|')
        if len(p)>=2:return p[1]
    return t

def ff(x):
    try:return float(x)
    except:return float('nan')

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--m8',required=True); ap.add_argument('--freeze',required=True); ap.add_argument('--panel',required=True); ap.add_argument('--outdir',required=True); a=ap.parse_args(); os.makedirs(a.outdir,exist_ok=True)
    qmeta={}; qorder=[]
    with open(a.freeze,newline='') as f:
        r=csv.DictReader(f,delimiter='\t')
        for x in r:qmeta[x['query']]=x;qorder.append(x['query'])
    pacc={}; ptarget={}; pkey={}
    with open(a.panel,newline='') as f:
        r=csv.DictReader(f,delimiter='\t')
        for x in r:
            pacc[x['accession']]=x
            if x.get('afsp_target'): ptarget[x['afsp_target']]=x
            if x.get('afsp_key'): pkey[x['afsp_key']]=x
    fields=['query','target','evalue','bits','qcov','tcov','alnlen','qlen','tlen','fident','alntmscore']
    hits=collections.defaultdict(list); rawrows=[]
    if os.path.exists(a.m8) and os.path.getsize(a.m8)>0:
        with open(a.m8,newline='') as f:
            r=csv.DictReader(f,delimiter='\t',fieldnames=fields)
            for row in r:
                fq=normq(row['query'])
                if fq not in qmeta: raise SystemExit(f'Unmapped query id: {row["query"]}')
                pm=ptarget.get(row['target']) or pkey.get(row['target']) or pacc.get(normacc(row['target']))
                if pm is None: raise SystemExit(f'Unmapped RT-panel target: {row["target"]}')
                row['frozen_query']=fq; row['target_accession']=pm['accession']; row['target_description']=pm['description']; row['target_annotation_class']=pm['annotation_class']; row['target_hmm_models']=pm['hmm_models']; row['target_source_arms']=pm['source_arms']
                tm=ff(row['alntmscore']); qc=ff(row['qcov']); tc=ff(row['tcov'])
                tier='RT_PANEL_WEAK_DESCRIPTIVE'
                if tm>=0.50 and qc>=0.50 and tc>=0.50:tier='RT_FOLD_STRONG_REQUIRES_VALIDATION'
                elif tm>=0.40 and qc>=0.30 and tc>=0.30:tier='RT_FOLD_SUGGESTIVE_REQUIRES_VALIDATION'
                row['structural_tier']=tier
                hits[fq].append(row);rawrows.append(row)
    def ek(h):
        try:return float(h['evalue'])
        except:return 1e300
    def tk(h):
        try:return float(h['alntmscore'])
        except:return -1
    summaries=[]
    for q in qorder:
        hs=hits.get(q,[])
        best_e=min(hs,key=ek) if hs else None
        best_t=max(hs,key=tk) if hs else None
        strong=[h for h in hs if h['structural_tier'].startswith('RT_FOLD_STRONG')]
        sugg=[h for h in hs if h['structural_tier'].startswith('RT_FOLD_SUGGESTIVE')]
        label='NO_RT_PANEL_HIT_NOT_NEGATIVE'
        if hs: label='RT_PANEL_WEAK_DESCRIPTIVE_ONLY'
        if sugg: label='RT_FOLD_SUGGESTIVE_REQUIRES_VALIDATION'
        if strong: label='RT_FOLD_STRONG_REQUIRES_VALIDATION'
        m=qmeta[q]
        summaries.append({
          'query':q,'aa_length':m.get('aa_length',''),'candidate_key_n':m.get('candidate_key_n',''),'candidate_keys':m.get('candidate_keys',''),'strata':m.get('strata',''),'min_orf_to_candidate_distance_bp':m.get('min_orf_to_candidate_distance_bp',''),
          'n_rt_panel_hits':len(hs),'n_rt_fold_strong':len(strong),'n_rt_fold_suggestive':len(sugg),
          'best_evalue_target':best_e['target'] if best_e else '','best_evalue_accession':best_e['target_accession'] if best_e else '','best_evalue':best_e['evalue'] if best_e else '','best_evalue_alntmscore':best_e['alntmscore'] if best_e else '','best_evalue_qcov':best_e['qcov'] if best_e else '','best_evalue_tcov':best_e['tcov'] if best_e else '','best_evalue_description':best_e['target_description'] if best_e else '',
          'best_tm_target':best_t['target'] if best_t else '','best_tm_accession':best_t['target_accession'] if best_t else '','best_alntmscore':best_t['alntmscore'] if best_t else '','best_tm_qcov':best_t['qcov'] if best_t else '','best_tm_tcov':best_t['tcov'] if best_t else '','best_tm_description':best_t['target_description'] if best_t else '','best_tm_hmm_models':best_t['target_hmm_models'] if best_t else '',
          'triage_label':label})
    sp=os.path.join(a.outdir,'STAGE05D4R5B_P1_RT_STRUCTURAL_TRIAGE.tsv')
    with open(sp,'w',newline='') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=list(summaries[0].keys()));w.writeheader();w.writerows(summaries)
    hp=os.path.join(a.outdir,'STAGE05D4R5B_ALL_RT_PANEL_HITS.tsv')
    if rawrows:
        with open(hp,'w',newline='') as f:
            w=csv.DictWriter(f,delimiter='\t',fieldnames=list(rawrows[0].keys()));w.writeheader();w.writerows(rawrows)
    else: open(hp,'w').write('query\ttarget\n')
    counts=collections.Counter(x['triage_label'] for x in summaries)
    S={'stage':'Stage05D4-R5B','n_P1':len(qorder),'n_foldseek_rows':len(rawrows),'n_queries_with_any_rt_panel_hit':sum(x['n_rt_panel_hits']>0 for x in summaries),'n_queries_with_strong_rt_fold':sum(x['n_rt_fold_strong']>0 for x in summaries),'n_queries_with_suggestive_rt_fold':sum(x['n_rt_fold_suggestive']>0 for x in summaries),'triage_counts':dict(counts),'candidate_filtering_performed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,'claim_boundary':'RT-focused structural matches are fold-compatibility evidence only; they do not establish RT activity or ART identity. No-hit is not a negative classification.'}
    with open(os.path.join(a.outdir,'STAGE05D4R5B_SUMMARY.json'),'w') as f:json.dump(S,f,indent=2)
    G=dict(S);G.update({'gate':'PASS_P1_RT_FOCUSED_STRUCTURAL_PANEL_SEARCH_COMPLETE','gate_pass':True,'reference_panel_selected_from_candidate_results':False})
    with open(os.path.join(a.outdir,'GATE_STAGE05D4R5B.json'),'w') as f:json.dump(G,f,indent=2)
    print(json.dumps(G,indent=2))
if __name__=='__main__':main()
