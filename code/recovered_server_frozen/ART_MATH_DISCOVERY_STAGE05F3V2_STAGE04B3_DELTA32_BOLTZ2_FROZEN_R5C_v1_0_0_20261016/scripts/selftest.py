#!/usr/bin/env python3
import csv, hashlib, json, subprocess, sys, tempfile
from pathlib import Path
root=Path(__file__).resolve().parent.parent
with tempfile.TemporaryDirectory() as td:
    d=Path(td); tsv=d/'s.tsv'; faa=d/'s.faa'; gate=d/'g.json'; ofa=d/'o.faa'; aud=d/'a.tsv'; inp=d/'inp'; man=d/'m.tsv'
    cols=['family_rep_id','priority_candidate_id','priority_candidate_accession','priority_candidate_aa_length','shortlisted','shortlist_rank']
    rows=[]
    with tsv.open('w',newline='') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=cols); w.writeheader()
        for i in range(1,33):
            cid=f'ACC{i}|q{i}|+|f0|0-180|aa60|o{i}'
            r={'family_rep_id':f'G{i:08d}','priority_candidate_id':cid,'priority_candidate_accession':f'ACC{i}','priority_candidate_aa_length':'60','shortlisted':'True','shortlist_rank':str(i)}; rows.append(r); w.writerow(r)
    with faa.open('w') as f:
        for r in rows: f.write('>'+r['priority_candidate_id']+'\n'+'M'*60+'\n')
    def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
    g={'gate_pass':True,'gate':'PASS_F3V2_STAGE04B2_HISTORICAL_F3B3A_EXACT_RULE_PROJECTED_TO_DELTA_ONLY_FAMILIES','structure_priority_shortlist_family_n':32,'known_rt_annotation_used_for_ranking':False,'novelty_used_for_ranking':False,'structure_used_for_ranking':False,'truth_used_for_ranking':False,'retuning_authorized':False,'scientific_threshold_changed':False,'shortlist_tsv_sha256':sha(tsv),'shortlist_fasta_sha256':sha(faa)}
    gate.write_text(json.dumps(g))
    subprocess.check_call([sys.executable,str(root/'scripts/prepare_delta32_handoff.py'),str(gate),str(tsv),str(faa),str(ofa),str(aud)])
    subprocess.check_call([sys.executable,str(root/'scripts/prepare_inputs_historical_exact.py'),str(ofa),str(aud),str(inp),str(man)])
    with man.open() as f: m=list(csv.DictReader(f,delimiter='\t'))
    assert len(m)==32 and int(m[0]['shortlist_rank'])==1 and int(m[-1]['shortlist_rank'])==32
print('PASS_STAGE04B3_HANDOFF_AND_HISTORICAL_INPUT_SELFTEST')
