#!/usr/bin/env python3
import argparse,csv,json
from pathlib import Path
ap=argparse.ArgumentParser(); ap.add_argument('--registry',required=True); ap.add_argument('--map',required=True); ap.add_argument('--out',required=True); a=ap.parse_args()
reg={r['control_id']:r for r in csv.DictReader(open(a.registry),delimiter='\t')}
rows=list(csv.DictReader(open(a.map),delimiter='\t'))
seen=set(); outrows=[]
for r in rows:
    cid=r['control_id'];
    if cid in seen: raise SystemExit(f'[FAIL] duplicate mapping {cid}')
    seen.add(cid)
    if cid not in reg: raise SystemExit(f'[FAIL] mapping control absent C0 registry: {cid}')
    ra=reg[cid].get('whole_sequence_accession','')
    if ra != r['registry_accession']:
        raise SystemExit(f'[FAIL] registry accession mismatch {cid}: registry={ra!r} map={r["registry_accession"]!r}')
    fa=r['fetch_accession']
    if not fa or fa.startswith(('NC_','NZ_','NM_','NR_','NP_','WP_','YP_')):
        raise SystemExit(f'[FAIL] fetch accession is not an INSDC source accession: {cid} {fa}')
    outrows.append({'control_id':cid,'registry_accession':ra,'fetch_accession':fa,'relation':r['relation']})
out={'gate':'PASS_STAGE04C0_R2_ACCESSION_MAP_VALIDATED','gate_pass':True,'mapped_records_n':len(rows),'records':outrows}
Path(a.out).write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
