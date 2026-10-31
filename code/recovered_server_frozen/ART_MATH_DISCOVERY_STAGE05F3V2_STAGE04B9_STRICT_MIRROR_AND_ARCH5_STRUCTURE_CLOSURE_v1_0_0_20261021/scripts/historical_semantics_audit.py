#!/usr/bin/env python3
import argparse,csv,gzip,json,math
from collections import defaultdict,Counter
from pathlib import Path

QUAL_DEFS=[
 'BOTH_1E2_INCLUSIVE',
 'IE_1E2_INCLUSIVE',
 'FULL_1E2_INCLUSIVE',
 'BOTH_1E2_NONFORMAL_ONLY',
 'IE_1E2_NONFORMAL_ONLY',
 'FULL_1E2_NONFORMAL_ONLY',
]

def qpass(name,full,ie):
    formal=(full<=1e-4 and ie<=1e-4)
    if name=='BOTH_1E2_INCLUSIVE': return full<=1e-2 and ie<=1e-2
    if name=='IE_1E2_INCLUSIVE': return ie<=1e-2
    if name=='FULL_1E2_INCLUSIVE': return full<=1e-2
    if name=='BOTH_1E2_NONFORMAL_ONLY': return full<=1e-2 and ie<=1e-2 and not formal
    if name=='IE_1E2_NONFORMAL_ONLY': return ie<=1e-2 and not formal
    if name=='FULL_1E2_NONFORMAL_ONLY': return full<=1e-2 and not formal
    raise KeyError(name)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--r4-map',required=True)
    ap.add_argument('--r3-domtblout',required=True)
    ap.add_argument('--r3-pair-hits',required=True)
    ap.add_argument('--out',required=True)
    a=ap.parse_args()

    expected={}
    with open(a.r4_map,newline='') as f:
        r=csv.DictReader(f,delimiter='\t')
        need={'candidate_id','formal_model_n','qualifying_model_n','best_domain_i_evalue'}
        if not need.issubset(set(r.fieldnames or [])):
            raise RuntimeError('Historical R4 map missing '+str(sorted(need-set(r.fieldnames or []))))
        for row in r: expected[row['candidate_id']]=row
    if len(expected)!=10598:
        raise RuntimeError(f'Historical map expected 10598 candidates, got {len(expected)}')

    formal_profiles=defaultdict(set)
    qual_profiles={name:defaultdict(set) for name in QUAL_DEFS}
    bestdom={}
    with open(a.r3_domtblout,errors='ignore') as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            x=line.split()
            if len(x)<22: continue
            profile=x[0]; q=x[3]
            if q not in expected: continue
            try:
                full=float(x[6]); ie=float(x[12]); score=float(x[13])
                ali_from=int(x[17]); ali_to=int(x[18]); env_from=int(x[19]); env_to=int(x[20])
            except Exception:
                continue
            if full<=1e-4 and ie<=1e-4:
                formal_profiles[q].add(profile)
                rec={'profile':profile,'full_evalue':full,'domain_i_evalue':ie,'score':score,
                     'ali_from':ali_from,'ali_to':ali_to,'env_from':env_from,'env_to':env_to}
                key=(ie,full,profile,ali_from,ali_to)
                old=bestdom.get(q)
                if old is None or key < (old['domain_i_evalue'],old['full_evalue'],old['profile'],old['ali_from'],old['ali_to']):
                    bestdom[q]=rec
            for name in QUAL_DEFS:
                if qpass(name,full,ie): qual_profiles[name][q].add(profile)

    formal_mismatch=[]; bestie_mismatch=[]
    for q,row in expected.items():
        got=len(formal_profiles.get(q,set())); want=int(row['formal_model_n'])
        if got!=want:
            formal_mismatch.append((q,want,got))
        if q not in bestdom:
            bestie_mismatch.append((q,row['best_domain_i_evalue'],'NO_FORMAL_DOMAIN'))
        else:
            want=float(row['best_domain_i_evalue']); gotie=float(bestdom[q]['domain_i_evalue'])
            if want!=gotie:
                bestie_mismatch.append((q,want,gotie))

    qual_mismatch_n={}
    qual_examples={}
    matching=[]
    for name in QUAL_DEFS:
        bad=[]
        for q,row in expected.items():
            got=len(qual_profiles[name].get(q,set())); want=int(row['qualifying_model_n'])
            if got!=want:
                if len(bad)<20: bad.append((q,want,got))
        nbad=sum(len(qual_profiles[name].get(q,set()))!=int(row['qualifying_model_n']) for q,row in expected.items())
        qual_mismatch_n[name]=nbad; qual_examples[name]=bad
        if nbad==0: matching.append(name)

    # Inspect historical pair-hit semantics as an independent audit descriptor.
    pair_classes=Counter(); pair_fields=[]; pair_formal_counts=defaultdict(set); pair_qual_counts=defaultdict(set)
    with gzip.open(a.r3_pair_hits,'rt') as f:
        r=csv.DictReader(f,delimiter='\t'); pair_fields=list(r.fieldnames or [])
        qcol=next((c for c in ['query_id','candidate_id','query'] if c in pair_fields),None)
        ccol='class' if 'class' in pair_fields else None
        pcol=next((c for c in ['target_id','target','profile','model_id','best_rvt_target'] if c in pair_fields),None)
        for i,row in enumerate(r):
            if ccol: pair_classes[row.get(ccol,'')]+=1
            if qcol and ccol and pcol and row.get(qcol) in expected:
                q=row[qcol]; p=row[pcol]; cls=row.get(ccol,'')
                if cls=='FORMAL': pair_formal_counts[q].add(p)
                if cls in ('FORMAL','QUALIFYING'): pair_qual_counts[q].add(p)
    pair_formal_mismatch=None; pair_qual_mismatch=None
    if pair_formal_counts:
        pair_formal_mismatch=sum(len(pair_formal_counts.get(q,set()))!=int(row['formal_model_n']) for q,row in expected.items())
        pair_qual_mismatch=sum(len(pair_qual_counts.get(q,set()))!=int(row['qualifying_model_n']) for q,row in expected.items())

    # The historical scan protocol explicitly used BOTH thresholds for exploratory rows.
    canonical='BOTH_1E2_INCLUSIVE'
    passed=(len(formal_mismatch)==0 and len(bestie_mismatch)==0 and canonical in matching)
    out={
      'stage':'Stage05F3V2-Stage04B2-historical-model-semantics-audit',
      'gate':'PASS_F3B3A_MODEL_COUNT_SEMANTICS_RECONSTRUCTED_FROM_HISTORICAL_DOMTBLOUT' if passed else 'HOLD_F3B3A_MODEL_COUNT_SEMANTICS_MISMATCH',
      'gate_pass':passed,
      'historical_candidate_n':len(expected),
      'formal_definition':'distinct RVT profile targets having >=1 domain with full E<=1e-4 AND i-E<=1e-4',
      'formal_model_count_mismatch_n':len(formal_mismatch),
      'formal_model_count_mismatch_examples':formal_mismatch[:20],
      'best_formal_domain_ie_mismatch_n':len(bestie_mismatch),
      'best_formal_domain_ie_mismatch_examples':bestie_mismatch[:20],
      'qualifying_definition_candidates_mismatch_n':qual_mismatch_n,
      'qualifying_definition_matching_candidates':matching,
      'qualifying_definition_used_for_delta':canonical if canonical in matching else None,
      'historical_pair_hit_fields':pair_fields,
      'historical_pair_hit_class_counts':dict(pair_classes),
      'pair_hit_formal_model_count_mismatch_n':pair_formal_mismatch,
      'pair_hit_formal_plus_qualifying_count_mismatch_n':pair_qual_mismatch,
      'retuning_authorized':False,
      'science_change':False,
      'claim_boundary':'This audit reconstructs the pre-existing F3B3A model-count semantics from historical frozen R3/R4 evidence before applying the historical rule to the native 54G delta operational-family universe.'
    }
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    with open(a.out,'w') as f: json.dump(out,f,indent=2); f.write('\n')
    print(json.dumps(out,indent=2))
    if not passed: raise SystemExit(2)

if __name__=='__main__': main()
