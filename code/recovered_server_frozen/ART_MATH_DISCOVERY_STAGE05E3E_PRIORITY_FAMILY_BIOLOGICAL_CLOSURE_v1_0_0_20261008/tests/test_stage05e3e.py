#!/usr/bin/env python3
import importlib.util,tempfile
from pathlib import Path
spec=importlib.util.spec_from_file_location('e3e',Path(__file__).resolve().parents[1]/'scripts/stage05e3e.py')
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
assert m.revcomp('ACGTN')=='NACGT'
assert m.topology_and_gap(100,200,'+',50,80)==('UPSTREAM',20)
assert m.topology_and_gap(100,200,'-',220,230)==('UPSTREAM',20)
assert m.topology_and_gap(100,200,'+',150,160)==('OVERLAP',0)
seq='AAAACGTTTTACGTGGG'; hits=m.scan_motif(seq,'ACGT'); assert len(hits)>=2
ident,qc,tc,aln=m.global_identity('AAAA','AAAA'); assert ident==1 and qc==1 and tc==1 and aln==4
with tempfile.TemporaryDirectory() as td:
    p=Path(td)/'x.tblout'
    p.write_text('target - query - 1e-10 50.0 0.0 1e-10 50.0 0.0 1.0 1 0 0 1 1 1 1 desc here\n')
    r=m.parse_tblout(p); assert len(r)==1 and r[0]['target_name']=='target' and r[0]['full_evalue']==1e-10
print('[PASS] Stage05E3E primitive/schema selftest')
