#!/usr/bin/env python3
from __future__ import annotations
import csv, hashlib, json, sys
from pathlib import Path

TARGET_Q='SF_00172800'
TARGET_SHA='73f4e1517b9e8df706aefcf7fa60df9ecda1e0ed210b45eaf84834a0e2e42239'
EXPECTED_ACCS={'PX404140.1','PX418225.1'}

def read_tsv(p):
    with open(p,newline='') as f: return list(csv.DictReader(f,delimiter='\t'))

def read_fasta(p):
    d={}; name=None; seq=[]
    with open(p) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if name is not None: d[name]=''.join(seq)
                name=line[1:].split()[0]; seq=[]
            else: seq.append(line)
    if name is not None: d[name]=''.join(seq)
    return d

def find(base,name):
    p=Path(base)/'outputs'/name
    if p.exists(): return p
    hits=list(Path(base).rglob(name))
    if len(hits)==1: return hits[0]
    raise SystemExit(f'cannot resolve {name} under {base}; hits={len(hits)}')

def main():
    if len(sys.argv)!=3:
        raise SystemExit('usage: prepare_stage05f2c.py F1_DIR F2B_DIR')
    f1=Path(sys.argv[1]).resolve(); f2b=Path(sys.argv[2]).resolve()
    out=Path('outputs'); out.mkdir(exist_ok=True)
    inp=Path('inputs/boltz2'); inp.mkdir(parents=True,exist_ok=True)
    g1=json.loads(find(f1,'STAGE05F1_GATE.json').read_text())
    g2=json.loads(find(f2b,'STAGE05F2B_GATE.json').read_text())
    if not g1.get('gate_pass') or not g2.get('gate_pass'):
        raise SystemExit('upstream F1/F2B gate not PASS')
    fa=read_fasta(find(f1,'STAGE05F1_UNIQUE_SIXFRAME_PEPTIDES.faa'))
    if TARGET_Q not in fa: raise SystemExit(f'{TARGET_Q} absent from F1 FASTA')
    seq=fa[TARGET_Q]
    sha=hashlib.sha256(seq.encode()).hexdigest()
    if sha!=TARGET_SHA or len(seq)!=205:
        raise SystemExit(f'frozen target mismatch: len={len(seq)} sha={sha}')
    q=read_tsv(find(f2b,'STAGE05F2B_F2C_FOLLOWUP_QUEUE.tsv'))
    rows=[r for r in q if r.get('target_query_id')==TARGET_Q]
    accs={r.get('sequence_id') for r in rows}
    if accs!=EXPECTED_ACCS:
        raise SystemExit(f'expected accessions {sorted(EXPECTED_ACCS)}, got {sorted(accs)}')
    for r in rows:
        if r.get('f2b_locus_class')!='LONG_SINGLE_FRAGMENT_RT_LIKE':
            raise SystemExit('frozen F2B class mismatch')
    (out/'STAGE05F2C_TARGET.faa').write_text(f'>{TARGET_Q} sha256={sha}\n{seq}\n')
    yml=(
        'version: 1\n'
        'sequences:\n'
        '  - protein:\n'
        '      id: A\n'
        f'      sequence: {seq}\n'
    )
    (inp/f'{TARGET_Q}.yaml').write_text(yml)
    cols=['sequence_id','candidate_keys','target_query_id','target_aa_length','target_rvt_profile','target_full_evalue','target_query_coverage','target_model_coverage','target_yxdd_like','f2b_locus_class']
    with (out/'STAGE05F2C_LOCUS_INSTANCES.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(rows)
    gate={
        'stage':'Stage05F2C','version':'v1.0.0','phase':'PREPARE','gate':'PASS_STAGE05F2C_FROZEN_TARGET_PREPARED','gate_pass':True,
        'target_query_id':TARGET_Q,'target_aa_sha256':sha,'target_aa_length':len(seq),'target_instance_n':len(rows),'target_accessions':sorted(accs),
        'upstream_f1_gate':g1.get('gate'),'upstream_f2b_gate':g2.get('gate'),
        'structure_prediction_performed':False,'structural_rt_compatibility_evaluated':False,
        'claim_boundary':'Preparation only. No structural or biochemical RT claim.'
    }
    (out/'STAGE05F2C_PREPARE_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate,indent=2))
if __name__=='__main__': main()
