#!/usr/bin/env python3
import json,sys
from pathlib import Path
p=Path(sys.argv[1])/'outputs/GATE_STAGE05D4.json'
g=json.loads(p.read_text())
assert g['gate']=='PASS_N1N2_REMOTE_HOMOLOGY_UNIVERSE_FROZEN'
assert g['n_N1N2_candidates']==602 and g['n_nonformal_broad_rt_unique_aa']==13062
print('[PASS] Stage05D4-R1 source gate/count selftest')
