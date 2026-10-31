#!/usr/bin/env python3
import json,sys,hashlib
from pathlib import Path
root=Path(sys.argv[1]); pol=json.loads((root/'config/C1_POLICY.json').read_text()); rows=[]
for cid in pol['eligible_control_ids']:
 p=root/'outputs'/cid/'C1_SCORE_FREEZE_GATE.json'
 if not p.exists(): raise SystemExit(f'[FAIL] missing score gate {cid}')
 d=json.loads(p.read_text())
 if not d.get('gate_pass') or d.get('truth_used') is not False: raise SystemExit(f'[FAIL] invalid gate {cid}')
 rows.append({'control_id':cid,'accession':d['accession'],'length_bp':d['length_bp'],'n_genomic_windows':d['n_genomic_windows'],'n_union_regions':d['n_union_regions'],'window_scores_sha256':d['window_scores_sha256'],'union_regions_sha256':d['union_regions_sha256']})
out={'gate':'PASS_STAGE04C1_RT_CONTROL_TWO_ARM_SCORE_FREEZE_COMPLETE','gate_pass':True,'scientific_status':pol['scientific_status'],'controls_scored_n':len(rows),'controls':rows,'truth_used_before_score_freeze':False,'unblinding_authorized_in_this_package':False,'retuning_authorized':False}
(root/'outputs/C1_FINAL_GATE.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
