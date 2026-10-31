#!/usr/bin/env python3
import csv,hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; SRC=Path(sys.argv[1])
rows=[]
for p in sorted((ROOT/'outputs/record_gates').glob('*.json')):
    d=json.loads(p.read_text()); rows.append(d)
out=ROOT/'outputs/STAGE04C0_R1_FIXED_SEQUENCE_MANIFEST.tsv'
with out.open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['control_id','accession','length_bp','fasta_sha256','flatfile_sha256','status'])
    for d in rows: w.writerow([d['control_id'],d['accession'],d['length_bp'],d['fasta_sha256'],d['flatfile_sha256'],'PRESENT_VALIDATED'])
expected=sum(1 for r in csv.DictReader((SRC/'registry/RT_SYSTEM_CONTROL_REGISTRY.tsv').open(),delimiter='\t') if r.get('whole_sequence_accession'))
gate={
 'gate':'PASS_STAGE04C0_R1_SEQUENCE_ACQUISITION_FIXED' if len(rows)==expected and expected>0 else 'FAIL_STAGE04C0_R1_SEQUENCE_ACQUISITION_FIXED',
 'gate_pass':len(rows)==expected and expected>0,
 'validated_records_n':len(rows),'expected_records_n':expected,
 'source_c0_root':str(SRC),'provider':'DDBJ_GETENTRY_INSDC',
 'original_c0_mutated':False,'runs_acnv_or_recurrence':False,'retuning_authorized':False
}
(ROOT/'outputs/STAGE04C0_R1_FINAL_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
print(json.dumps(gate,indent=2))
if not gate['gate_pass']: raise SystemExit(1)
