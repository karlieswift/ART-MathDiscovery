#!/usr/bin/env python3
import csv,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'outputs'
a=list(csv.DictReader(open(OUT/'REMAINING_BROAD_RT_HIT_AUDIT.tsv'),delimiter='\t'))
o=list(csv.DictReader(open(OUT/'REMAINING_OFFICIAL_ANNOTATION_AUDIT.tsv'),delimiter='\t'))
passn=sum(r.get('status')=='PASS' for r in o)
known=sum(r.get('status')=='PASS' and r.get('official_product','').lower() not in ('','hypothetical protein','uncharacterized protein') for r in o)
g={
 'stage':'Stage05D3',
 'gate':'PASS_REMAINING_STAGE05D1_BROAD_RT_HITS_AUDITED_WITH_NOVELTY_HOLD',
 'gate_pass': len(a)==3 and len(o)==3,
 'n_frozen_remaining_broad_rt_targets':len(a),
 'official_annotation_pass_n':passn,
 'official_known_product_n':known,
 'dna_score_recomputed':False,
 'dna_rank_recomputed':False,
 'candidate_filtering_performed':False,
 'myrt_threshold_changed':False,
 'new_ART_confirmed_n':0,
 'wet_lab_ready_n':0,
 'claim_boundary':'Stage05D3 closes the remaining formal broad-RVT hits. HMM, motif, distance, and official annotation are post-discovery interpretation only; divergent RTs may exist among HMM-negative candidates.'
}
(OUT/'GATE_STAGE05D3.json').write_text(json.dumps(g,indent=2)+'\n')
print(json.dumps(g,indent=2))
