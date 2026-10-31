#!/usr/bin/env python3
import json,sys,os
out=sys.argv[1]
S=json.load(open(os.path.join(out,'STAGE05D4R4_SUMMARY.json')))
G={
 'stage':'Stage05D4-R4',
 'gate':'PASS_PROSTT5_FOLDSEEK_REMOTE_STRUCTURE_TRIAGE_FROZEN',
 'gate_pass': S['n_tier_A_queries']==350,
 'n_tier_A_queries':S['n_tier_A_queries'],
 'n_queries_with_any_structural_hit':S['n_queries_with_any_structural_hit'],
 'n_queries_with_rt_text_structural_hit':S['n_queries_with_rt_text_structural_hit'],
 'dna_score_recomputed':False,
 'dna_rank_recomputed':False,
 'candidate_filtering_performed':False,
 'prostt5_3di_is_full_3d_structure':False,
 'structural_hit_used_as_identity_assignment':False,
 'no_hit_interpreted_as_negative':False,
 'claim_boundary':'ProstT5/Foldseek hits are fast post-discovery structural-alphabet triage only. RT-text hits require full-structure, architecture, known-family dereplication, and functional validation before ART claims.'
}
json.dump(G,open(os.path.join(out,'GATE_STAGE05D4R4.json'),'w'),indent=2)
print(json.dumps(G,indent=2))
