#!/usr/bin/env python3
import sys,csv,hashlib
inp,outfa,maptsv=sys.argv[1:4]
seq2id={}; rows=[]; name=None; seq=[]
def flush():
  if name is None:return
  s=''.join(seq).rstrip('*').upper(); h=hashlib.sha256(s.encode()).hexdigest()
  uid=seq2id.get(s)
  if uid is None:
    uid=f'F3V2U{len(seq2id)+1:08d}'; seq2id[s]=uid
  rows.append((name,uid,h,len(s)))
for raw in open(inp):
  line=raw.strip()
  if not line:continue
  if line.startswith('>'):
    flush(); name=line[1:].split()[0];seq=[]
  else:seq.append(line)
flush()
with open(outfa,'w') as f:
  for s,uid in sorted(((s,u) for s,u in seq2id.items()),key=lambda x:x[1]): f.write(f'>{uid}\n{s}\n')
with open(maptsv,'w',newline='') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['formal_query_id','exact_unique_id','aa_sha256','aa_length']);w.writerows(rows)
print('formal_instance_n =',len(rows));print('exact_aa_unique_n =',len(seq2id))
