#!/usr/bin/env python3
import argparse,csv,gzip,json,re,shutil
from pathlib import Path
RT_PAT=re.compile(r'\b(reverse transcriptase|retrotranscriptase|rna[- ]dependent dna polymerase|rvt|group ii intron.*reverse|retr(on|on-associated).*reverse|diversity[- ]generating retroelement)\b',re.I)
def sp_desc(p):
    d={}; op=gzip.open if str(p).endswith('.gz') else open
    with op(p,'rt',errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                h=line[1:].strip(); acc=h.split('|')[1] if '|' in h and len(h.split('|'))>2 else h.split()[0]; d[acc]=h
    return d
def acc(t):
    z=t.split()[0]
    if '|' in z:
        p=z.split('|');
        if len(p)>1:return p[1]
    m=re.search(r'(?:AF-)?([A-Z0-9]+)-F\d+',z)
    return m.group(1) if m else z

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--freeze',required=True); ap.add_argument('--boltz',required=True); ap.add_argument('--m8',required=True); ap.add_argument('--sprot',required=True); ap.add_argument('--out',required=True); a=ap.parse_args()
    out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    fr=list(csv.DictReader(open(Path(a.freeze)/'STAGE05D4R5_P1_FREEZE.tsv'),delimiter='\t')); assert len(fr)==43
    conf={}
    for p in Path(a.boltz).glob('predictions/*/confidence_*_model_0.json'):
        q=p.parent.name; x=json.load(open(p)); conf[q]=x
    desc=sp_desc(a.sprot); hits={}
    with open(a.m8) as f:
        for line in f:
            z=line.rstrip('\n').split('\t');
            if len(z)<11: continue
            q,t=z[0],z[1]; a0=acc(t); ds=desc.get(a0,'')
            r={'target':t,'target_accession':a0,'evalue':float(z[2]),'bits':float(z[3]),'qcov':float(z[4]),'tcov':float(z[5]),'alnlen':int(float(z[6])),'qlen':int(float(z[7])),'tlen':int(float(z[8])),'fident':float(z[9]),'tmscore':float(z[10]),'description':ds,'rt_text':bool(RT_PAT.search(ds))}
            hits.setdefault(q,[]).append(r)
    rows=[]
    for x in fr:
        q=x['query']; hs=hits.get(q,[]); rhs=[h for h in hs if h['rt_text']]; c=conf.get(q,{})
        best=hs[0] if hs else {}; br=rhs[0] if rhs else {}
        rows.append({'query':q,'aa_length':x['aa_length'],'candidate_key_n':x['candidate_key_n'],'candidate_keys':x['candidate_keys'],'strata':x['strata'],'min_orf_to_candidate_distance_bp':x['min_orf_to_candidate_distance_bp'],'confidence_score':c.get('confidence_score',''),'ptm':c.get('ptm',''),'complex_plddt':c.get('complex_plddt',''),'n_full3d_hits':len(hs),'best_target':best.get('target',''),'best_evalue':best.get('evalue',''),'best_tmscore':best.get('tmscore',''),'best_qcov':best.get('qcov',''),'best_tcov':best.get('tcov',''),'best_description':best.get('description',''),'n_rt_text_full3d_hits':len(rhs),'best_rt_target':br.get('target',''),'best_rt_evalue':br.get('evalue',''),'best_rt_tmscore':br.get('tmscore',''),'best_rt_qcov':br.get('qcov',''),'best_rt_description':br.get('description',''),'triage_label':'FULL3D_RT_TEXT_HIT_REQUIRES_VALIDATION' if rhs else ('FULL3D_CURATED_NON_RT_TRIAGE' if hs else 'NO_FULL3D_AFSP_HIT_NOT_NEGATIVE')})
    fields=list(rows[0]);
    with open(out/'STAGE05D4R5_P1_FULL3D_TRIAGE.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(rows)
    cnt={}
    for r in rows: cnt[r['triage_label']]=cnt.get(r['triage_label'],0)+1
    s={'stage':'Stage05D4-R5','n_P1':43,'n_boltz_predictions':len(conf),'n_queries_with_full3d_hit':sum(bool(hits.get(r['query'])) for r in fr),'n_queries_with_rt_text_full3d_hit':sum(any(h['rt_text'] for h in hits.get(r['query'],[])) for r in fr),'triage_counts':cnt,'candidate_filtering_performed':False,'claim_boundary':'Boltz-2 full-coordinate predictions plus Foldseek are post-discovery protein interpretation. RT-text structural hits require catalytic-core/topology, DNA architecture, known-family dereplication and ultimately functional validation before ART claims.'}
    json.dump(s,open(out/'STAGE05D4R5_SUMMARY.json','w'),indent=2)
    gate={**s,'gate':'PASS_P1_BOLTZ2_FULL3D_AND_AFSP_FOLDSEEK_COMPLETE' if len(conf)==43 else 'FAIL_INCOMPLETE_BOLTZ_P1','gate_pass':len(conf)==43,'dna_score_recomputed':False,'dna_rank_recomputed':False,'p1_compute_priority_only':True,'single_sequence_fallback_used':False}
    json.dump(gate,open(out/'GATE_STAGE05D4R5.json','w'),indent=2)
    print(json.dumps(s,indent=2))
if __name__=='__main__':main()
