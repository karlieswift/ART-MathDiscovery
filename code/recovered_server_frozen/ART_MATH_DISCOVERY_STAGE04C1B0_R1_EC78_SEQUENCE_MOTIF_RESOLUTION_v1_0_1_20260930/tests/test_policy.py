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
print('[PASS] Stage04C1B0-R1 no-score-access / EC78 sequence-motif policy selftest')
