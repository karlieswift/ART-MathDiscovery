import json, pathlib
r=pathlib.Path(__file__).resolve().parents[1]
p=json.loads((r/'config/C1B0_POLICY.json').read_text())
assert p['score_inputs_read'] is False
assert p['score_mutation_authorized'] is False
assert p['retuning_authorized'] is False
assert p['truth_resolution_occurs_after_C1_score_freeze'] is True
s=(r/'scripts/resolve_truth.py').read_text()
assert 'C1_WINDOW_FAMILY_SCORES' not in s and 'C1_UNION_REGIONS' not in s
print('[PASS] Stage04C1B0 no-score-access policy selftest')
