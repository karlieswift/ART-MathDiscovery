#!/usr/bin/env python3
import json,csv
from pathlib import Path
R=Path(__file__).resolve().parents[1]
p=json.loads((R/'config/C0_POLICY.json').read_text())
assert p['primary_method_for_future_C1']==['ACNV_FAMILY','GENERIC_RECURRENCE']
assert p['proposal_q_threshold']==0.005
assert p['top_k_survival_cutoff'] is False
assert p['retuning_authorized'] is False
assert p['score_execution_in_C0'] is False
rows=list(csv.DictReader((R/'registry/RT_SYSTEM_CONTROL_REGISTRY.tsv').open(),delimiter='\t'))
assert len(rows)==10
assert sum(r['genomic_repeat_array_expected']=='YES' for r in rows)==5
assert sum(r['genomic_repeat_array_expected']=='NO' for r in rows)==5
print('[PASS] Stage04C0 policy selftest')
