#!/usr/bin/env python3
import csv,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
pol=json.loads((root/'config/HOTFIX_POLICY.json').read_text())
rows={r['control_id']:r for r in csv.DictReader((root/'config/ACCESSION_MAP.tsv').open(),delimiter='\t')}
assert rows['DRT2_KPN']['fetch_accession']=='AP006725.1'
assert int(rows['DRT2_KPN']['expected_length_bp'])==5248520
assert int(rows['DRT2_KPN']['expected_length_bp'])+224152==5472672
assert pol['expected_exact_lengths_bp']['DRT2_KPN']==5248520
assert not pol['run_acnv_or_recurrence'] and not pol['retuning_authorized']
print('[PASS] Stage04C0-R3 corrected DRT2 chromosome-length policy selftest')
