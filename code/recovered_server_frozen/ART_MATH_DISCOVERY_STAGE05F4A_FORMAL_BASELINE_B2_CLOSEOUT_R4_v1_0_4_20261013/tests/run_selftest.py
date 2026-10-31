#!/usr/bin/env python3
from pathlib import Path
import ast, json
here=Path(__file__).resolve().parents[1]
for p in list((here/'scripts').glob('*.py')):
    ast.parse(p.read_text(),filename=str(p))
for p in list(here.glob('*.sh')):
    t=p.read_text(); assert '01_run_score_freeze' not in t and '00_preflight_and_freeze' not in t
c=json.loads((here/'config/BENCHMARK_CONTRACT.json').read_text())
assert c['primary_rt_neighborhood_radius_bp']==10000
assert len(c['primary_positive_set'])==8
pt=(here/'scripts/postfreeze_unblind.py').read_text()
assert 'truth_dual_arm' in pt and 'BASELINE_EXACT_ARRAY_SUBSET_SUMMARY.tsv' in pt and 'BASELINE_STRICT_EXTERNAL_ARRAY_SUBSET_SUMMARY.tsv' in pt
rt=(here/'scripts/resolve_truth.py').read_text()
assert 'EXACT_CORE_PRIMARY_DESCRIPTION' in rt and 'RT_rrtT_STM14_4641' in rt and 'RT_Ec78' in rt and 'RT_DRT2' in rt
print('[PASS] R4 B2-only package syntax/contract/truth-resolution selftest')
print('[PASS] no B0/B1 scoring runners are present; frozen R2 B1 must be imported')
print('[PASS] pre-registered 8-case primary panel preserved; 7-case exact-array subset is descriptive only')
