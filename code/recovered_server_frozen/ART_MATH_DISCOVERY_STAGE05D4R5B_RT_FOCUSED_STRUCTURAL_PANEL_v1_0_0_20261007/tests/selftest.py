#!/usr/bin/env python3
import json,re
p=json.load(open('config/POLICY.json'))
assert p['query_universe_expected']==43
assert p['hmm_top_per_profile']==10
assert p['reference_panel_selected_from_candidate_results'] is False
assert p['strong_rt_fold_rule'].startswith('alntmscore >= 0.50')
assert p['suggestive_rt_fold_rule'].startswith('alntmscore >= 0.40')
print('[PASS] Stage05D4-R5B policy selftest')
