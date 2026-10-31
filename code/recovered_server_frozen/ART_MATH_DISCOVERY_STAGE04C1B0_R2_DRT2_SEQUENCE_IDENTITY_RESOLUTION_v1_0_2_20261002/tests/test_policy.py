#!/usr/bin/env python3
import json,sys
from pathlib import Path
r=Path(__file__).resolve().parents[1]
p=json.loads((r/'config/HOTFIX_POLICY.json').read_text())
assert p['score_inputs_read'] is False
assert p['score_mutation_authorized'] is False
assert p['retuning_authorized'] is False
assert p['network_fetch_performed'] is False
assert 'YADD' in p['ec78_primary_rule']
assert p['drt2_source_product_annotation_required'] is False
assert p['drt2_source_protein_id_required'] is False
assert 'YCDD' in p['drt2_primary_rule']
assert '425' in p['drt2_primary_rule']
print('[PASS] Stage04C1B0-R2 no-score-access / EC78+DRT2 sequence-identity policy selftest')
