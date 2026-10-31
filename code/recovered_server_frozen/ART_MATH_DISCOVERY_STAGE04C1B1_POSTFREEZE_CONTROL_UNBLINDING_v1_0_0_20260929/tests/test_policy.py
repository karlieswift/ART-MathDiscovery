import json,pathlib
r=pathlib.Path(__file__).resolve().parents[1]
p=json.loads((r/'config/C1B1_POLICY.json').read_text())
assert p['retuning_authorized'] is False
assert p['score_recomputation_authorized'] is False
assert p['ranking_recomputation_authorized'] is False
s=(r/'scripts/unblind.py').read_text()
assert 'q<=0.005' not in s and 'q <=' not in s
print('[PASS] Stage04C1B1 frozen-score/no-retune policy selftest')
