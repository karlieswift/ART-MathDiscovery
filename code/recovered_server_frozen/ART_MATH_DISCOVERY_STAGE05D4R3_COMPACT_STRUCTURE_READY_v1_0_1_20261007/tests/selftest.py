import json, sys
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"config/POLICY.json"
x=json.loads(p.read_text())
assert x["source_n1_n2_candidates"]==602
assert x["source_nr90_queries"]==9245
assert x["source_no_swissprot_hit_queries"]==8062
assert x["primary_tier"]["aa_length_min"]==180
assert x["primary_tier"]["aa_length_max"]==700
assert x["primary_tier"]["max_orf_to_math_candidate_distance_bp"]==1000
assert x["extended_tier"]["max_orf_to_math_candidate_distance_bp"]==2500
assert x["rules"]["dna_score_recomputed"] is False
assert x["rules"]["candidate_filtering_performed"] is False
print("[PASS] Stage05D4-R3 policy freeze selftest")
