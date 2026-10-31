#!/usr/bin/env python3
import sys,csv
ids=set()
with open(sys.argv[1]) as f:
  r=csv.DictReader(f,delimiter='\t')
  for x in r:
    if x['formal']=='True': ids.add(x['query_id'])
out=open(sys.argv[3],'w')
name=None; seq=[]
def flush():
  if name in ids:
    out.write('>'+name+'\n'+''.join(seq)+'\n')
for raw in open(sys.argv[2]):
  line=raw.strip()
  if not line: continue
  if line.startswith('>'):
    flush(); name=line[1:].split()[0]; seq=[]
  else: seq.append(line)
flush(); out.close()
