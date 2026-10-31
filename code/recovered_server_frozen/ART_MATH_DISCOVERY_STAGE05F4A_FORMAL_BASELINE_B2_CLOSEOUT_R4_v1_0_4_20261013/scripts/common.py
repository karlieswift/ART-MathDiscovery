#!/usr/bin/env python3
from __future__ import annotations
import csv, hashlib, json, os, re, sys
from pathlib import Path

def sha256_file(p: Path) -> str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()

def read_tsv(p: Path):
    with p.open(newline='') as f: return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(p: Path, rows, fields):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('w', newline='') as f:
        w=csv.DictWriter(f, delimiter='\t', fieldnames=fields, lineterminator='\n'); w.writeheader(); w.writerows(rows)

def norm_accession(s:str)->str:
    return s.strip().upper()

def fasta_records(p: Path):
    name=None; seq=[]
    with p.open(errors='replace') as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if name is not None: yield name, ''.join(seq).upper()
                name=line[1:]; seq=[]
            else: seq.append(re.sub('[^A-Za-z]','',line))
    if name is not None: yield name, ''.join(seq).upper()

def seq_sha(seq:str)->str:
    return hashlib.sha256(seq.encode()).hexdigest()

def safe_json(p:Path, obj):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, indent=2, sort_keys=True)+"\n")
