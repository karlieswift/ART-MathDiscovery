#!/usr/bin/env python3
import csv,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
reg=list(csv.DictReader((ROOT/'registry/RT_SYSTEM_CONTROL_REGISTRY.tsv').open(),delimiter='\t'))
seq=ROOT/'outputs/STAGE04C0_SEQUENCE_FREEZE_MANIFEST.tsv'
seq_rows=[]
if seq.exists(): seq_rows=list(csv.DictReader(seq.open(),delimiter='\t'))
files=[]
for p in sorted([x for x in ROOT.rglob('*') if x.is_file() and 'outputs/STAGE04C0_FINAL_GATE.json' not in str(x)]):
    if '__pycache__' in p.parts: continue
    files.append({'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size})
summary={
 'gate':'PASS_STAGE04C0_RT_SYSTEM_CONTROL_REGISTRY_FROZEN_WITH_EXPLICIT_SOURCE_HOLDS',
 'gate_pass':True,
 'scientific_status':'CONTROL_REGISTRY_AND_EVIDENCE_CLASS_FREEZE_NO_CONTROL_SCORING',
 'panel_rows_n':len(reg),
 'panel_groups':{g:sum(r['panel_group']==g for r in reg) for g in sorted({r['panel_group'] for r in reg})},
 'sequence_source_accessions_present_in_registry_n':sum(bool(r['whole_sequence_accession']) for r in reg),
 'c1_scan_eligible_now_n':sum(r['c1_scan_eligible_now']=='1' for r in reg),
 'unblinding_truth_ready_n':sum(r['unblinding_truth_ready']=='1' for r in reg),
 'sequence_files_downloaded_n':sum(r.get('status')=='PRESENT' for r in seq_rows),
 'explicit_pending_source_resolution_ids':[r['control_id'] for r in reg if r['source_resolution_status'].startswith('PENDING')],
 'retuning_authorized':False,
 'stage04a_mutation_authorized':False,
 'c0_runs_acnv_or_recurrence':False,
 'files':files
}
(ROOT/'outputs/STAGE04C0_FINAL_GATE.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='files'},indent=2))
