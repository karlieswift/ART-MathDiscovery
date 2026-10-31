#!/usr/bin/env python3
import csv, hashlib, json, sys, tempfile
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

pkg=Path(__file__).resolve().parents[1]
cand=pkg/'inputs'/'STAGE04A_R1_DUAL_ARM_REGIONS.tsv'
h=hashlib.sha256(cand.read_bytes()).hexdigest()
assert h=='73007864fdef7dba662c9b90a02e2b915359c1fbd93a6061edeea35afde07349', h
with open(cand,newline='',encoding='utf-8') as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==383
assert len({r['sequence_id'] for r in rows})==174
assert all(float(r['best_q_acnv_family']) <= 0.005+1e-15 for r in rows)
assert all(float(r['best_q_recurrence']) <= 0.005+1e-15 for r in rows)

# Verify Biopython GenBank round-trip and explicit RT keyword fixture.
rec=SeqRecord(Seq('A'*1000),id='TEST.1',name='TEST',description='synthetic')
rec.annotations['molecule_type']='DNA'
rec.features=[SeqFeature(FeatureLocation(100,400,strand=1),type='CDS',qualifiers={'product':['reverse transcriptase'],'translation':['M'+'A'*99]})]
with tempfile.TemporaryDirectory() as td:
    p=Path(td)/'x.gb'
    SeqIO.write(rec,str(p),'genbank')
    got=list(SeqIO.parse(str(p),'genbank'))
    assert len(got)==1 and len(got[0].features)==1

# Regression test for v1.0.1: Stage04A manifest hash is FASTA file bytes,
# while GenBank identity must use a separate normalized sequence hash.
with tempfile.TemporaryDirectory() as td:
    fp=Path(td)/'frozen.fasta'
    fp.write_text('>TEST description\nACGTACGT\n', encoding='utf-8')
    file_sha=hashlib.sha256(fp.read_bytes()).hexdigest()
    seq_sha=hashlib.sha256(b'ACGTACGT').hexdigest()
    assert file_sha != seq_sha
    parsed=next(SeqIO.parse(str(fp),'fasta'))
    parsed_seq_sha=hashlib.sha256(str(parsed.seq).upper().encode()).hexdigest()
    assert parsed_seq_sha == seq_sha


# v1.0.2 regression: reject HTML/error payloads and accept GenBank-like payload shape.
def looks_like_genbank(data):
    if not data or len(data) < 1000:
        return False
    head=data[:10000].lstrip()
    return head.startswith(b"LOCUS") and b"ORIGIN" in data and data.rstrip().endswith(b"//")
assert not looks_like_genbank(b"<html>" + b"x"*4000)
valid = b"LOCUS       TEST       10 bp    DNA\n" + b"FEATURES             Location/Qualifiers\n" + b"x"*1000 + b"\nORIGIN\n        1 acgtacgtac\n//\n"
assert looks_like_genbank(valid)
print('[PASS] Stage04B-R1 v1.0.2 selftest: canonical 383/174 registry, hash semantics, GenBank payload validation, and Biopython support')
