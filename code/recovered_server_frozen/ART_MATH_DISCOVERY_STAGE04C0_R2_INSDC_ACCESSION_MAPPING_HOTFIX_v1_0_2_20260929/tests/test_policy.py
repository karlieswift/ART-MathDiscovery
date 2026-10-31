#!/usr/bin/env python3
import csv,json,sys,tempfile,subprocess
from pathlib import Path
R=Path(__file__).resolve().parents[1]
p=json.loads((R/'config/HOTFIX_POLICY.json').read_text())
assert p['provider']=='DDBJ_GETENTRY_INSDC' and not p['run_acnv_or_recurrence'] and not p['retuning_authorized']
rows=list(csv.DictReader((R/'config/ACCESSION_MAP.tsv').open(),delimiter='\t'))
assert len(rows)==4
M={r['control_id']:r for r in rows}
assert M['SEN2_14028S']['registry_accession']=='NC_016856.1' and M['SEN2_14028S']['fetch_accession']=='CP001363.1'
assert M['DRT2_KPN']['registry_accession']=='NC_012731.1' and M['DRT2_KPN']['fetch_accession']=='AP006725.1'
for r in rows: assert not r['fetch_accession'].startswith(('NC_','NZ_'))
print('[PASS] Stage04C0-R2 accession mapping/policy selftest')
