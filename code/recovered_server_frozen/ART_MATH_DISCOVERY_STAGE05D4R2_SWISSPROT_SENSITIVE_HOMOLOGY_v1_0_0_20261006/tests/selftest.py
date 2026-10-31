from pathlib import Path
import json,re,tempfile
p=Path(__file__).resolve().parents[1]
j=json.load(open(p/'config/POLICY.json'))
assert j['dna_score_recomputed'] is False and j['candidate_filtering_authorized'] is False
assert j['mmseqs_sensitivity']==7.5 and j['evalue']==0.001
pat=re.compile(r'reverse transcript|RNA[- ]directed DNA polymerase|retron',re.I)
assert pat.search('Reverse transcriptase') and pat.search('RNA-directed DNA polymerase')
print('[PASS] Stage05D4-R2 policy/text-triage selftest')
