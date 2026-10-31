from pathlib import Path
import csv, hashlib, json, re

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):
            h.update(b)
    return h.hexdigest()

def md5(path):
    h=hashlib.md5()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):
            h.update(b)
    return h.hexdigest()

def read_tsv(path):
    with open(path,newline='',encoding='utf-8') as f:
        return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(path, rows, fields):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

def dump_json(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n',encoding='utf-8')

def strict_overlap(a0,a1,b0,b1):
    return max(a0,b0) < min(a1,b1)

def overlap_bp(a0,a1,b0,b1):
    return max(0,min(a1,b1)-max(a0,b0))

def gap_bp(a0,a1,b0,b1):
    if strict_overlap(a0,a1,b0,b1): return 0
    if a1 <= b0: return b0-a1
    return a0-b1

def parse_fasta(path):
    name=None; seq=[]
    with open(path,encoding='utf-8',errors='ignore') as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if name is not None: yield name,''.join(seq).upper()
                name=line[1:].split()[0]; seq=[]
            else: seq.append(line)
    if name is not None: yield name,''.join(seq).upper()

def revcomp(s):
    return s.translate(str.maketrans('ACGTNacgtn','TGCANtgcan'))[::-1]

def normalize_system_label(s):
    s=s.strip()
    aliases={
      'Alistipes-putredinis':'Alistipes-putredinis',
      'Bacteroides-coprophilus':'Bacteroides-coprophilus',
      'Bacteroides-sp-9-1-42FAA':'Bacteroides-sp-9-1-42FAA',
      'Bacteroides-sp-2-1-16':'Bacteroides-sp-2-1-16',
      'Blautia-hansenii':'Blautia-hansenii'
    }
    if s in aliases: return aliases[s]
    for a in aliases:
        if a in s: return aliases[a]
    return s
