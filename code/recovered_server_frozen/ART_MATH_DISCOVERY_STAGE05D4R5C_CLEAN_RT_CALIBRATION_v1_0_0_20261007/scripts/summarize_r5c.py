#!/usr/bin/env python3
import argparse,csv,json,math,os,re,statistics,collections

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

def quant(vals,p):
    vals=sorted(v for v in vals if math.isfinite(v))
    if not vals:return None
    if len(vals)==1:return vals[0]
    x=(len(vals)-1)*p; a=int(math.floor(x)); b=int(math.ceil(x))
    return vals[a] if a==b else vals[a]*(b-x)+vals[b]*(x-a)

def read_m8(path):
    fields=['query','target','evalue','bits','qcov','tcov','alnlen','qlen','tlen','fident','alntmscore']
    if not os.path.exists(path) or os.path.getsize(path)==0:return []
    return list(csv.DictReader(open(path),delimiter='\t',fieldnames=fields))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--p1m8',required=True); ap.add_argument('--posm8',required=True); ap.add_argument('--negm8',required=True)
    ap.add_argument('--p1freeze',required=True); ap.add_argument('--pos',required=True); ap.add_argument('--neg',required=True); ap.add_argument('--outdir',required=True)
    a=ap.parse_args();os.makedirs(a.outdir,exist_ok=True)
    p1=list(csv.DictReader(open(a.p1freeze),delimiter='\t')); pos=list(csv.DictReader(open(a.pos),delimiter='\t')); neg=list(csv.DictReader(open(a.neg),delimiter='\t'))
    pos_by_acc={r['accession']:r for r in pos}; neg_by_acc={r['accession']:r for r in neg}
    p1_by_q={r['query']:r for r in p1}
    def score_rows(rows,mode):
        out=collections.defaultdict(list)
        for r in rows:
            q=normq(r['query']); t=normacc(r['target'])
            tm,qc,tc=ff(r['alntmscore']),ff(r['qcov']),ff(r['tcov'])
            rr=dict(r);rr['query_norm']=q;rr['target_accession']=t;rr['tier']=tier(tm,qc,tc)
            if mode=='pos' and q.startswith('AF-'): q=normacc(q)
            out[q].append(rr)
        return out
    # DB-vs-DB query names are AF-accession model names; normalize accession.
    def dbgroup(rows):
        out=collections.defaultdict(list)
        for r in rows:
            q=normacc(r['query']); t=normacc(r['target']); rr=dict(r); rr['query_accession']=q;rr['target_accession']=t;rr['tier']=tier(ff(r['alntmscore']),ff(r['qcov']),ff(r['tcov']));out[q].append(rr)
        return out
    pg=collections.defaultdict(list)
    for r in read_m8(a.p1m8):
        q=normq(r['query']); t=normacc(r['target']); rr=dict(r);rr['query_norm']=q;rr['target_accession']=t;rr['tier']=tier(ff(r['alntmscore']),ff(r['qcov']),ff(r['tcov']));pg[q].append(rr)
    posg=dbgroup(read_m8(a.posm8)); negg=dbgroup(read_m8(a.negm8))
    # Known RT calibration: exclude self target accession.
    pos_summary=[]
    for r in pos:
        acc=r['accession']; hs=[h for h in posg.get(acc,[]) if h['target_accession']!=acc]
        best=max(hs,key=lambda h:ff(h['alntmscore'])) if hs else None
        pos_summary.append({'accession':acc,'annotation_class':r['annotation_class'],'aa_length':r['aa_length'],'n_nonself_hits':len(hs),'best_nonself_target':best['target_accession'] if best else '','best_nonself_alntmscore':best['alntmscore'] if best else '','best_nonself_qcov':best['qcov'] if best else '','best_nonself_tcov':best['tcov'] if best else '','best_nonself_tier':best['tier'] if best else 'NO_HIT'})
    neg_summary=[]
    for r in neg:
        acc=r['accession'];hs=negg.get(acc,[]);best=max(hs,key=lambda h:ff(h['alntmscore'])) if hs else None
        neg_summary.append({'accession':acc,'matched_rt_accession':r['matched_rt_accession'],'aa_length':r['aa_length'],'length_delta':r['length_delta'],'n_rt_target_hits':len(hs),'best_rt_target':best['target_accession'] if best else '','best_alntmscore':best['alntmscore'] if best else '','best_qcov':best['qcov'] if best else '','best_tcov':best['tcov'] if best else '','best_tier':best['tier'] if best else 'NO_HIT'})
    p1_summary=[]
    for r in p1:
        q=r['query'];hs=pg.get(q,[]);best=max(hs,key=lambda h:ff(h['alntmscore'])) if hs else None
        p1_summary.append({'query':q,'aa_length':r.get('aa_length',''),'candidate_key_n':r.get('candidate_key_n',''),'strata':r.get('strata',''),'min_orf_to_candidate_distance_bp':r.get('min_orf_to_candidate_distance_bp',''),'n_clean_rt_hits':len(hs),'best_clean_rt_target':best['target_accession'] if best else '','best_alntmscore':best['alntmscore'] if best else '','best_qcov':best['qcov'] if best else '','best_tcov':best['tcov'] if best else '','best_tier':best['tier'] if best else 'NO_HIT'})
    def write(path,rows):
        with open(path,'w',newline='') as f:
            w=csv.DictWriter(f,delimiter='\t',fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)
    write(os.path.join(a.outdir,'STAGE05D4R5C_KNOWN_RT_CALIBRATION.tsv'),pos_summary)
    write(os.path.join(a.outdir,'STAGE05D4R5C_MATCHED_NONRT_CALIBRATION.tsv'),neg_summary)
    write(os.path.join(a.outdir,'STAGE05D4R5C_P1_CLEAN_RT_TRIAGE.tsv'),p1_summary)
    pos_good=sum(x['best_nonself_tier'] in ('SUGGESTIVE','STRONG') for x in pos_summary)
    neg_good=sum(x['best_tier'] in ('SUGGESTIVE','STRONG') for x in neg_summary)
    p1_good=sum(x['best_tier'] in ('SUGGESTIVE','STRONG') for x in p1_summary)
    pos_scores=[ff(x['best_nonself_alntmscore']) for x in pos_summary if x['best_nonself_alntmscore']!='']
    neg_scores=[ff(x['best_alntmscore']) for x in neg_summary if x['best_alntmscore']!='']
    p1_scores=[ff(x['best_alntmscore']) for x in p1_summary if x['best_alntmscore']!='']
    sens=pos_good/len(pos_summary); fpr=neg_good/len(neg_summary); informative=(sens>=0.75 and fpr<=0.05)
    S={'stage':'Stage05D4-R5C','n_clean_rt_controls':len(pos_summary),'n_matched_nonrt_controls':len(neg_summary),'n_P1':len(p1_summary),'known_rt_suggestive_or_strong_sensitivity':sens,'matched_nonrt_suggestive_or_strong_rate':fpr,'n_P1_suggestive_or_strong_clean_rt_hits':p1_good,'known_rt_best_nonself_alntmscore_median':statistics.median(pos_scores) if pos_scores else None,'known_rt_best_nonself_alntmscore_p10':quant(pos_scores,0.10),'matched_nonrt_best_alntmscore_p95':quant(neg_scores,0.95),'P1_best_alntmscore_max':max(p1_scores) if p1_scores else None,'calibration_informative':informative,'interpretation':'INFORMATIVE_FOR_DEPRIORITIZATION' if informative and p1_good==0 else ('CALIBRATION_LIMITED' if not informative else 'P1_RT_FOLD_SIGNAL_PRESENT_REQUIRES_VALIDATION'),'candidate_filtering_performed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False,'claim_boundary':'Calibration evaluates the structural-screen operating characteristics against known RT and matched non-RT controls. It does not assign RT activity or ART identity.'}
    json.dump(S,open(os.path.join(a.outdir,'STAGE05D4R5C_SUMMARY.json'),'w'),indent=2)
    G=dict(S);G.update({'gate':'PASS_CLEAN_RT_STRUCTURAL_CALIBRATION_COMPLETE','gate_pass':True,'reference_selection_uses_candidate_results':False})
    json.dump(G,open(os.path.join(a.outdir,'GATE_STAGE05D4R5C.json'),'w'),indent=2)
    print(json.dumps(G,indent=2))
if __name__=='__main__': main()
