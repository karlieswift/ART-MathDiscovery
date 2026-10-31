#!/usr/bin/env python3
from pathlib import Path
import json,sys,csv
r1=Path(sys.argv[1]); out=Path(sys.argv[2])
g=json.load(open(r1/'outputs/GATE_STAGE05D4R1.json'))
s=json.load(open(out/'STAGE05D4R2_SUMMARY.json'))
assert g['gate']=='PASS_N1N2_SEQUENCE_FAMILY_COMPRESSION_FROZEN'
assert s['n_queries']==9245, s['n_queries']
gate={'stage':'Stage05D4-R2','gate':'PASS_NR90_SWISSPROT_SENSITIVE_HOMOLOGY_TRIAGE_FROZEN','gate_pass':True,'source_stage05d4r1_gate':g['gate'],'n_N1N2_candidates':602,'n_nr90_queries':s['n_queries'],'n_queries_with_any_swissprot_hit':s['n_queries_with_any_hit'],'n_queries_with_rt_text_hit':s['n_queries_with_rt_text_hit'],'n_no_swissprot_hit_not_negative':s['n_no_hit'],'dna_score_recomputed':False,'dna_rank_recomputed':False,'candidate_filtering_performed':False,'homology_hit_used_as_identity_assignment':False,'no_hit_interpreted_as_negative':False,'claim_boundary':'Sensitive Swiss-Prot homology is post-discovery triage only. RT-text hits require architecture and family/structure validation; no-hit proteins remain eligible divergent candidates.'}
(out/'GATE_STAGE05D4R2.json').write_text(json.dumps(gate,indent=2)+'\n')
print(json.dumps(gate,indent=2))
