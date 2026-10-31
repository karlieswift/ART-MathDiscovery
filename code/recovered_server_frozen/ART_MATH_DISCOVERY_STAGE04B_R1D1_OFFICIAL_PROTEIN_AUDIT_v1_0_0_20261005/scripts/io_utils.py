from pathlib import Path
import contextlib,csv,datetime,gzip,hashlib,json,os,platform,sys,tempfile

def digest(p,algorithm='sha256'):
    h=hashlib.new(algorithm)
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def shatext(t):return hashlib.sha256(t.encode('ascii')).hexdigest()

def read_tsv(p):
    with Path(p).open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(p,rows,fields):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)

def js(p,obj):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile('w',dir=p.parent,delete=False,encoding='utf-8') as f:
        json.dump(obj,f,indent=2,sort_keys=True);f.write('\n');tmp=f.name
    os.replace(tmp,p)

def read_json(p):return json.loads(Path(p).read_text())

def check_hash(p,expected,label='file'):
    got=digest(p)
    if got!=expected:raise RuntimeError(f'{label} hash mismatch: {p}; expected={expected}; observed={got}')
    return got

def single_fasta(p):
    records=fasta_text(Path(p).read_text())
    if len(records)!=1:raise ValueError(f'Expected 1 record in {p}, got {len(records)}')
    return next(iter(records.items()))

def fasta_text(text):
    result={};key=None;buf=[]
    def save():
        if key is not None:
            if key in result:raise ValueError('Duplicate FASTA identifier: '+key)
            value=''.join(buf).upper()
            if not value:raise ValueError('Empty sequence: '+key)
            result[key]=value
    for line in text.splitlines():
        line=line.strip()
        if not line:continue
        if line.startswith('>'):
            save();key=line[1:].split()[0];buf=[]
        else:
            if key is None:raise ValueError('Sequence before FASTA header')
            if any(x.isspace() for x in line):raise ValueError('Whitespace inside sequence')
            buf.append(line)
    save()
    return result

def manifest_tree(root,out,exclude=()):
    root=Path(root);out=Path(out)
    rows=[]
    for p in sorted(root.rglob('*')):
        if not p.is_file() or p==out:continue
        rel=p.relative_to(root).as_posix()
        if any(rel==x or rel.startswith(x+'/') for x in exclude):continue
        rows.append(dict(path=rel,bytes=p.stat().st_size,sha256=digest(p)))
    write_tsv(out,rows,['path','bytes','sha256'])
    return digest(out)

def runtime_info():
    d=dict(python=sys.version,executable=sys.executable,platform=platform.platform())
    for name in ('numpy','pandas','Bio'):
        try:
            m=__import__(name);d[name]=m.__version__
        except ImportError:pass
    return d

@contextlib.contextmanager
def lock_dir(p):
    p=Path(p)
    try:p.mkdir()
    except FileExistsError:raise RuntimeError(f'Run lock exists: {p}; verify no process is active before removing this lock')
    try:yield
    finally:p.rmdir()
