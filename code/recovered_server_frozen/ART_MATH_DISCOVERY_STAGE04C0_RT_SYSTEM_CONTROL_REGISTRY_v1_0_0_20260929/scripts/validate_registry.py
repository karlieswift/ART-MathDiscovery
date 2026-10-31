#!/usr/bin/env python3
import csv, json, hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
reg=ROOT/'registry/RT_SYSTEM_CONTROL_REGISTRY.tsv'
evd=ROOT/'registry/EVIDENCE_MATRIX.tsv'
truth=ROOT/'truth/CONTROL_LOCUS_TRUTH_SEALED.tsv'
required=['control_id','panel_group','system_class','system_name','control_role','source_resolution_status','genomic_repeat_array_expected','experimental_evidence_tier','c1_scan_eligible_now','unblinding_truth_ready']
with reg.open() as f:
    rows=list(csv.DictReader(f, delimiter='\t'))
assert len(rows)==10, len(rows)
assert len({r['control_id'] for r in rows})==10
for r in rows:
    for k in required: assert r.get(k,'')!='', (r['control_id'],k)
    if r['c1_scan_eligible_now']=='1': assert r['whole_sequence_accession'], r['control_id']
with evd.open() as f: erows=list(csv.DictReader(f,delimiter='\t'))
assert {r['control_id'] for r in erows}=={r['control_id'] for r in rows}
with truth.open() as f: trows=list(csv.DictReader(f,delimiter='\t'))
assert {r['control_id'] for r in trows}=={r['control_id'] for r in rows}
summary={
 'gate':'PASS_STAGE04C0_REGISTRY_SCHEMA_AND_HOLD_POLICY',
 'panel_rows':len(rows),
 'groups':{g:sum(r['panel_group']==g for r in rows) for g in sorted({r['panel_group'] for r in rows})},
 'whole_sequence_accession_present_n':sum(bool(r['whole_sequence_accession']) for r in rows),
 'c1_scan_eligible_now_n':sum(r['c1_scan_eligible_now']=='1' for r in rows),
 'unblinding_truth_ready_n':sum(r['unblinding_truth_ready']=='1' for r in rows),
 'pending_parent_or_supplement_resolution_n':sum(r['source_resolution_status'].startswith('PENDING') for r in rows),
 'retuning_authorized':False
}
print(json.dumps(summary,indent=2))
(ROOT/'outputs').mkdir(exist_ok=True)
(ROOT/'outputs/STAGE04C0_REGISTRY_GATE.json').write_text(json.dumps(summary,indent=2)+'\n')
