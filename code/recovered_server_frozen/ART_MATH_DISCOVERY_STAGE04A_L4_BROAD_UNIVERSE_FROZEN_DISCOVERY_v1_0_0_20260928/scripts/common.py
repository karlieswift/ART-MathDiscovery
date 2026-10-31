#!/usr/bin/env python3
from pathlib import Path
import hashlib, json, re

def sha256_file(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def sha256_text(s): return hashlib.sha256(s.encode()).hexdigest()

def write_json(path,obj):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True)+"\n")

def parse_fasta(path):
    header=None; chunks=[]; nrec=0
    with open(path,errors='strict') as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                nrec += 1
                if nrec>1: raise ValueError(f'multi-record FASTA not allowed: {path}')
                header=line[1:].strip()
            else:
                if header is None: raise ValueError(f'sequence before FASTA header: {path}')
                chunks.append(line.upper())
    if nrec!=1: raise ValueError(f'expected exactly one FASTA record: {path}')
    seq=''.join(chunks)
    return header,seq

def safe_id(s):
    x=re.sub(r'[^A-Za-z0-9._-]+','_',str(s)).strip('._-')
    if not x: raise ValueError('empty safe sequence id')
    return x[:160]

def accession_root(s):
    x=str(s).split()[0]
    return x.split('.')[0]
