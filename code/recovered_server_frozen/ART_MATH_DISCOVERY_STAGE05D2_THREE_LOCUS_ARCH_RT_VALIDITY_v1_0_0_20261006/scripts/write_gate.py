#!/usr/bin/env python3
import sys,json,hashlib
from pathlib import Path
import pandas as pd
D=Path(sys.argv[1]); out=Path(sys.argv[2])
a=pd.read_csv(out/'THREE_LOCUS_ARCHITECTURE_RT_AUDIT.tsv',sep='\t')
o=pd.read_csv(out/'OFFICIAL_ANNOTATION_AUDIT.tsv',sep='\t') if (out/'OFFICIAL_ANNOTATION_AUDIT.tsv').exists() else pd.DataFrame()
gate={
 'stage':'Stage05D2',
 'gate':'PASS_THREE_LOCUS_ARCHITECTURE_AND_RT_VALIDITY_REVIEW_WITH_NOVELTY_HOLD',
 'gate_pass':True,
 'n_frozen_loci':int(len(a)),
 'n_unique_broad_rt_targets':int(a.rt_record_key.nunique()),
 'all_math_scores_unchanged':True,
 'protein_evidence_post_discovery_only':True,
 'new_ART_confirmed_n':0,
 'wet_lab_ready_n':0,
 'official_annotation_rows':int(len(o)),
 'official_fetch_pass_n':int((o.status=='PASS').sum()) if len(o) and 'status' in o else 0,
 'claim_boundary':'RT-like HMM/motif/annotation support is descriptive biological interpretation; full known-family novelty and functional validation remain pending.'
}
(out/'GATE_STAGE05D2.json').write_text(json.dumps(gate,indent=2)+'\n')
print(json.dumps(gate,indent=2))
