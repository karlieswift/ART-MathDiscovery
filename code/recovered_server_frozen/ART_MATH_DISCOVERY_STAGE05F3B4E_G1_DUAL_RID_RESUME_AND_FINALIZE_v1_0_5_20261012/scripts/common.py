#!/usr/bin/env python3
import csv, hashlib, json, re
from pathlib import Path

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''): h.update(b)
    return h.hexdigest()

def read_fasta(path):
    out={}; h=None; seq=[]; desc={}
    with open(path,encoding='utf-8',errors='ignore') as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if h is not None: out[h]=''.join(seq).upper()
                raw=line[1:]; h=raw.split()[0]; desc[h]=raw; seq=[]
            else: seq.append(line)
        if h is not None: out[h]=''.join(seq).upper()
    return out,desc

def wrap(s,w=80): return '\n'.join(s[i:i+w] for i in range(0,len(s),w))+'\n'

def read_tsv(path):
    with open(path,encoding='utf-8',errors='ignore') as f: return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(path, rows, fields=None):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    if fields is None:
        fields=list(rows[0]) if rows else []
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)

def b(x): return str(x).strip().lower()=='true'

def fnum(x, default=0.0):
    try: return float(x)
    except: return default

def inum(x, default=0):
    try: return int(float(x))
    except: return default

def family_from_desc(raw):
    m=re.search(r'(?:^|\s)family=([^\s]+)',raw)
    return m.group(1) if m else ''
