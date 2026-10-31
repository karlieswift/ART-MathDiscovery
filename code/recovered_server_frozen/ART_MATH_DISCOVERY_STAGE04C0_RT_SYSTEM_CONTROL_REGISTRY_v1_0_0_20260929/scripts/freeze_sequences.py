#!/usr/bin/env python3
import csv, hashlib, json, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
in_dir=ROOT/'inputs/sequences'
out=ROOT/'outputs/STAGE04C0_SEQUENCE_FREEZE_MANIFEST.tsv'
reg=list(csv.DictReader((ROOT/'registry/RT_SYSTEM_CONTROL_REGISTRY.tsv').open(),delimiter='\t'))
def seq_len(p):
    n=0
    for line in p.open():
        if not line.startswith('>'): n+=len(line.strip())
    return n
rows=[]
for r in reg:
    acc=r['whole_sequence_accession']
    if not acc: continue
    p=in_dir/f'{r["control_id"]}__{acc}.fasta'
    if not p.exists():
        rows.append([r['control_id'],acc,'MISSING','',''])
        continue
    b=p.read_bytes(); h=hashlib.sha256(b).hexdigest(); n=seq_len(p)
    rows.append([r['control_id'],acc,'PRESENT',str(n),h])
out.parent.mkdir(exist_ok=True)
with out.open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['control_id','accession','status','length_bp','sha256']); w.writerows(rows)
print(out)
