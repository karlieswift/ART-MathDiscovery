#!/usr/bin/env python3
import argparse,csv,json,math,os,re,statistics,collections
from pathlib import Path

def normacc(t):
    m=re.search(r'AF-([A-Za-z0-9]+)-F\d+-model',t)
    return m.group(1) if m else t

def normq(q): return re.sub(r'_model_0$','',q)
def ff(x):
    try:return float(x)
    except:return float('nan')
def tier(tm,qc,tc):
    if tm>=0.50 and qc>=0.50 and tc>=0.50:return 'STRONG'
    if tm>=0.40 and qc>=0.30 and tc>=0.30:return 'SUGGESTIVE'
    return 'WEAK'
def read_m8(path):
    fields=['query','target','evalue','bits','qcov','tcov','alnlen','qlen','tlen','fident','alntmscore']
    if not os.path.exists(path) or os.path.getsize(path)==0:return []
    return list(csv.DictReader(open(path),delimiter='\t',fieldnames=fields))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--freeze',required=True);ap.add_argument('--boltz',required=True);ap.add_argument('--m8',required=True);ap.add_argument('--r5csummary',required=True);ap.add_argument('--outdir',required=True);a=ap.parse_args();os.makedirs(a.outdir,exist_ok=True)
    fr=list(csv.DictReader(open(a.freeze),delimiter='\t')); assert len(fr)==40
    cal=json.load(open(a.r5csummary)); assert cal['calibration_informative'] is True and cal['known_rt_suggestive_or_strong_sensitivity']>=0.75 and cal['matched_nonrt_suggestive_or_strong_rate']<=0.05
    conf={}
    for p in Path(a.boltz).glob('predictions/*/confidence_*_model_0.json'):
        q=p.parent.name;conf[q]=json.load(open(p))
    if len(conf)!=40: raise SystemExit(f'expected 40 confidence jsons, got {len(conf)}')
    frozen={r['query'] for r in fr}
    if set(conf)!=frozen:
        missing=sorted(frozen-set(conf)); extra=sorted(set(conf)-frozen)
        raise SystemExit(f'confidence/freeze query mismatch: missing={len(missing)} extra={len(extra)} missing_examples={missing[:3]} extra_examples={extra[:3]}')
    hits=collections.defaultdict(list)
    for r in read_m8(a.m8):
        q=normq(r['query'])
        if q not in frozen:
            raise SystemExit(f'unmapped Foldseek query not in frozen P2: {r["query"]} -> {q}')
        tm,qc,tc=ff(r['alntmscore']),ff(r['qcov']),ff(r['tcov']);rr=dict(r);rr['target_accession']=normacc(r['target']);rr['tier']=tier(tm,qc,tc);hits[q].append(rr)
    rows=[]
    for r in fr:
        q=r['query'];hs=hits.get(q,[]);best=max(hs,key=lambda h:ff(h['alntmscore'])) if hs else None;c=conf[q]
        rows.append({'query':q,'aa_length':r['aa_length'],'candidate_key_n':r['candidate_key_n'],'candidate_keys':r['candidate_keys'],'strata':r['strata'],'min_orf_to_candidate_distance_bp':r['min_orf_to_candidate_distance_bp'],'confidence_score':c.get('confidence_score',''),'ptm':c.get('ptm',''),'complex_plddt':c.get('complex_plddt',''),'n_clean_rt_hits':len(hs),'best_clean_rt_target':best['target_accession'] if best else '','best_evalue':best['evalue'] if best else '','best_alntmscore':best['alntmscore'] if best else '','best_qcov':best['qcov'] if best else '','best_tcov':best['tcov'] if best else '','best_tier':best['tier'] if best else 'NO_HIT'})
    with open(os.path.join(a.outdir,'STAGE05D4R5D_P2_CLEAN_RT_TRIAGE.tsv'),'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    good=sum(r['best_tier'] in ('STRONG','SUGGESTIVE') for r in rows);scores=[ff(r['best_alntmscore']) for r in rows if r['best_alntmscore']!='']
    tiers=collections.Counter(r['best_tier'] for r in rows)
    S={'stage':'Stage05D4-R5D','n_P2':40,'n_boltz_predictions':40,'n_foldseek_rows':sum(len(v) for v in hits.values()),'n_queries_with_any_clean_rt_hit':sum(r['n_clean_rt_hits']>0 for r in rows),'n_P2_suggestive_or_strong_clean_rt_hits':good,'P2_best_alntmscore_max':max(scores) if scores else None,'triage_counts':dict(tiers),'r5c_known_rt_sensitivity':cal['known_rt_suggestive_or_strong_sensitivity'],'r5c_matched_nonrt_rate':cal['matched_nonrt_suggestive_or_strong_rate'],'r5c_known_rt_p10':cal['known_rt_best_nonself_alntmscore_p10'],'r5c_matched_nonrt_p95':cal['matched_nonrt_best_alntmscore_p95'],'calibration_informative':True,'interpretation':'P2_RT_FOLD_SIGNAL_PRESENT_REQUIRES_VALIDATION' if good else 'P2_INFORMATIVE_FOR_DEPRIORITIZATION','candidate_filtering_performed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,'claim_boundary':'P2 is a compute-priority extension. A calibrated RT-fold match is prioritization evidence only and does not establish RT activity or ART identity; no-hit remains unresolved rather than negative.'}
    json.dump(S,open(os.path.join(a.outdir,'STAGE05D4R5D_SUMMARY.json'),'w'),indent=2)
    G=dict(S);G.update({'gate':'PASS_P2_CALIBRATED_FULL3D_COMPLETE','gate_pass':True,'source_calibration_gate':'PASS_CLEAN_RT_STRUCTURAL_CALIBRATION_COMPLETE','scientific_settings_changed':False,'implementation_hotfix':'v1.0.2 postprocessing-only: import pathlib.Path required by summarize_p2; reuse existing 40 Boltz predictions and existing P2_CLEAN_RT.m8; no Foldseek rerun, no prediction rerun, no threshold change.','prediction_workers_rerun':False,'foldseek_rerun':False})
    json.dump(G,open(os.path.join(a.outdir,'GATE_STAGE05D4R5D.json'),'w'),indent=2)
    print(json.dumps(G,indent=2))
if __name__=='__main__':main()
