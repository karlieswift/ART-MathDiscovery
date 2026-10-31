#!/usr/bin/env python3
from pathlib import Path
import csv, hashlib, json, sys

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

root=Path(sys.argv[1] if len(sys.argv)>1 else 'boltz_out').resolve()
out=Path('outputs'); out.mkdir(exist_ok=True)
cifs=sorted(root.rglob('*.cif'))+sorted(root.rglob('*.mmcif'))
jsons=sorted(root.rglob('*.json'))
rows=[]
for p in cifs:
    rows.append({'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)})
with (out/'STAGE05F2C_BOLTZ_STRUCTURE_MANIFEST.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['path','bytes','sha256'],delimiter='\t'); w.writeheader(); w.writerows(rows)
status={
 'stage':'Stage05F2C','version':'v1.0.0','phase':'BOLTZ2_FULL3D',
 'gate':'PASS_STAGE05F2C_BOLTZ2_FULL3D_PREDICTION' if len(cifs)>=1 else 'FAIL_STAGE05F2C_NO_STRUCTURE_OUTPUT',
 'gate_pass':len(cifs)>=1,'structure_file_n':len(cifs),'json_file_n':len(jsons),'boltz_out':str(root),
 'target_query_id':'SF_00172800','target_aa_sha256':'73f4e1517b9e8df706aefcf7fa60df9ecda1e0ed210b45eaf84834a0e2e42239',
 'structural_rt_compatibility_evaluated':False,
 'claim_boundary':'Boltz-2 coordinate prediction completed. RT identity remains unresolved until frozen calibrated RT-structure comparison.'
}
(out/'STAGE05F2C_BOLTZ_GATE.json').write_text(json.dumps(status,indent=2)+'\n')
print(json.dumps(status,indent=2))
if not status['gate_pass']: sys.exit(2)
