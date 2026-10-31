"""Small, dependency-free IO and coordinate helpers for a frozen evidence run."""
from __future__ import annotations
import csv, hashlib, json, os, re, subprocess, sys
from pathlib import Path
from datetime import datetime, timezone

def sha_file(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576), b''): h.update(b)
    return h.hexdigest()

def sha_text(s): return hashlib.sha256(s.encode('utf-8')).hexdigest()
def canonical_json(x): return json.dumps(x, sort_keys=True, separators=(',',':'))
def read_json(p): return json.loads(Path(p).read_text(encoding='utf8'))
def now(): return datetime.now(timezone.utc).isoformat()
def atomic_text(p,s):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_name(p.name+'.partial');tmp.write_text(s,encoding='utf8');os.replace(tmp,p)
def write_json(p,x): atomic_text(p,json.dumps(x,sort_keys=True,indent=2)+'\n')
def read_tsv(p):
    with Path(p).open(encoding='utf8',newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(p,rows,fields):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+'.partial')
    with tmp.open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='raise');w.writeheader();w.writerows(rows)
    os.replace(tmp,p)

def fasta(p):
    desc=None; chunks=[]
    with Path(p).open(encoding='ascii') as f:
        for line in f:
            s=line.strip()
            if not s: continue
            if s.startswith('>'):
                if desc is not None: yield desc.split()[0],desc,''.join(chunks).upper()
                desc=s[1:];chunks=[]
                if not desc: raise ValueError('Empty FASTA header')
            else:
                if desc is None: raise ValueError('Sequence before header')
                chunks.append(''.join(s.split()))
    if desc is not None: yield desc.split()[0],desc,''.join(chunks).upper()

def fasta_dict(p):
    d={}
    for key,desc,seq in fasta(p):
        if key in d or not seq: raise ValueError('Duplicate/empty FASTA: '+key)
        d[key]=(desc,seq)
    return d

def write_fasta(p,rows):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+'.partial')
    with tmp.open('w',encoding='ascii',newline='\n') as f:
        for name,seq in rows:
            if not seq or any(c.isspace() for c in name): raise ValueError('Bad FASTA row '+name)
            f.write('>'+name+'\n')
            for i in range(0,len(seq),80): f.write(seq[i:i+80]+'\n')
    os.replace(tmp,p)

def rc(s): return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def overlaps(a,b): return int(a[0])<int(b[1]) and int(b[0])<int(a[1])
def distance(a,b): return max(0,int(a[0])-int(b[1]),int(b[0])-int(a[1]))
def merge_intervals(xs):
    out=[]
    for a,b in sorted((int(a),int(b)) for a,b in xs):
        if not a<b: raise ValueError('Bad interval')
        if out and a<=out[-1][1]: out[-1][1]=max(b,out[-1][1])
        else: out.append([a,b])
    return out

def covered_bp(interval,mask):
    a,b=interval
    return sum(max(0,min(b,y)-max(a,x)) for x,y in mask)

def parts_from_registry(row):
    return [(int(p['start0']),int(p['end0'])) for p in json.loads(row['cds_parts_json'])]

def run_cmd(cmd,log,timeout=None):
    log=Path(log);log.parent.mkdir(parents=True,exist_ok=True)
    with log.open('w',encoding='utf8') as f:
        f.write('COMMAND '+json.dumps([str(x) for x in cmd])+'\n');f.flush()
        cp=subprocess.run([str(x) for x in cmd],stdout=f,stderr=subprocess.STDOUT,timeout=timeout)
    if cp.returncode: raise RuntimeError('Command failed rc=%s; see %s'%(cp.returncode,log))

def assert_manifest(root, manifest):
    errors=[]
    for r in read_tsv(manifest):
        p=Path(root)/r['path']
        if not p.is_file():errors.append('missing:'+r['path']);continue
        if p.stat().st_size!=int(r['size_bytes']) or sha_file(p)!=r['sha256']:errors.append('changed:'+r['path'])
    if errors: raise RuntimeError('Input lock failed: '+'; '.join(errors[:10]))

def output_manifest(root,paths,dest):
    rr=[]
    for p in sorted(paths):
        if p.is_file():rr.append({'path':str(p.relative_to(root)),'size_bytes':p.stat().st_size,'sha256':sha_file(p)})
    write_tsv(dest,rr,['path','size_bytes','sha256'])
    return sha_file(dest)
