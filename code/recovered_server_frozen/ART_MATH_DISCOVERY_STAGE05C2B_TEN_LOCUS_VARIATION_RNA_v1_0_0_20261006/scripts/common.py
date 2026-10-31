from pathlib import Path
import csv,hashlib,json,os

def sha_bytes(b):return hashlib.sha256(b).hexdigest()
def sha_text(s):return sha_bytes(s.encode('ascii'))
def sha_file(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def read_tsv(p):
 with Path(p).open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(p,rows,fields=None):
 rows=list(rows);p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 if fields is None:
  if not rows:raise ValueError('empty rows require explicit fields')
  fields=list(rows[0])
 with p.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def read_json(p):return json.loads(Path(p).read_text())
def write_json(p,obj):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.tmp')
 t.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n');os.replace(t,p)
def fasta(p):
 out={};name=None;buf=[]
 for line in Path(p).read_text().splitlines():
  line=line.strip()
  if not line:continue
  if line.startswith('>'):
   if name is not None:out[name]=''.join(buf)
   name=line[1:].split()[0];buf=[]
   if name in out:raise ValueError('duplicate FASTA id '+name)
  else:
   if name is None:raise ValueError('sequence before header')
   buf.append(line.upper())
 if name is not None:
  if name in out:raise ValueError('duplicate FASTA id '+name)
  out[name]=''.join(buf)
 if not out:raise ValueError('empty FASTA '+str(p))
 return out
def write_fasta(p,seqs):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w') as f:
  for k,s in seqs.items():
   f.write('>'+k+'\n')
   for i in range(0,len(s),80):f.write(s[i:i+80]+'\n')
def rc(s):return s.translate(str.maketrans('ACGTU','TGCAA'))[::-1]
def intervals_union(xs):
 ans=[]
 for a,b in sorted(xs):
  if b<=a:continue
  if ans and a<=ans[-1][1]:ans[-1][1]=max(ans[-1][1],b)
  else:ans.append([a,b])
 return ans
def complement(a,b,xs):
 out=[];s=a
 for c,d in intervals_union((max(a,c),min(b,d)) for c,d in xs):
  if c>s:out.append((s,c))
  s=max(s,d)
 if s<b:out.append((s,b))
 return out
def windows(a,b,lengths=(120,240,480)):
 if 60<=b-a<120:return [(a,b,'SHORT_INTERGENIC')]
 out=[]
 for n in lengths:
  if b-a<n:continue
  ss=list(range(a,b-n+1,n//2))
  if ss[-1]!=b-n:ss.append(b-n)
  out.extend((s,s+n,'INTERGENIC_WINDOW_'+str(n)) for s in ss)
 return out
def validate_manifest(root,name='PACKAGE_MANIFEST.tsv'):
 root=Path(root);rows=read_tsv(root/name)
 for r in rows:
  p=root/r['path']
  if not p.is_file() or p.stat().st_size!=int(r['size_bytes']) or sha_file(p)!=r['sha256']:raise ValueError('manifest mismatch: '+r['path'])
 return len(rows)
