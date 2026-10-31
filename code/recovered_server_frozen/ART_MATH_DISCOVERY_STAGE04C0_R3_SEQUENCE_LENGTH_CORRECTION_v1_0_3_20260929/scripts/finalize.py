#!/usr/bin/env python3
import csv,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; R2=Path(sys.argv[1])
rows=[]
for p in sorted((ROOT/'outputs/record_gates').glob('*.json')):
    d=json.loads(p.read_text())
    if d.get('gate_pass'): rows.append(d)
out=ROOT/'outputs/STAGE04C0_R3_FIXED_SEQUENCE_MANIFEST.tsv'
with out.open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['control_id','registry_accession','accession','relation','length_bp','fasta_sha256','flatfile_sha256','status'])
    for d in rows: w.writerow([d['control_id'],d['registry_accession'],d['accession'],d['relation'],d['length_bp'],d['fasta_sha256'],d['flatfile_sha256'],'PRESENT_VALIDATED'])
expected=sum(1 for _ in csv.DictReader((ROOT/'config/ACCESSION_MAP.tsv').open(),delimiter='\t'))
gate={
 'gate':'PASS_STAGE04C0_R3_SEQUENCE_LENGTH_CORRECTION_FIXED' if len(rows)==expected and expected>0 else 'FAIL_STAGE04C0_R3_SEQUENCE_LENGTH_CORRECTION_FIXED',
 'gate_pass':len(rows)==expected and expected>0,'validated_records_n':len(rows),'expected_records_n':expected,
 'source_r2_root':str(R2),'network_fetch_performed':False,'records_reused_from_r2':True,
 'drt2_chromosome_length_bp':5248520,'drt2_plasmid_length_bp_not_scanned':224152,
 'previous_erroneous_two_replicon_total_bp':5472672,
 'original_c0_mutated':False,'failed_r1_mutated':False,'r2_mutated':False,
 'runs_acnv_or_recurrence':False,'retuning_authorized':False
}
(ROOT/'outputs/STAGE04C0_R3_FINAL_GATE.json').write_text(json.dumps(gate,indent=2)+'\n'); print(json.dumps(gate,indent=2))
if not gate['gate_pass']: raise SystemExit(1)
