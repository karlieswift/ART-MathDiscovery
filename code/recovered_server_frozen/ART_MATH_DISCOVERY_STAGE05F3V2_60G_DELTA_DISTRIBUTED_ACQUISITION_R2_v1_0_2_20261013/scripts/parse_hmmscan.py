#!/usr/bin/env python3
import sys,csv
p,out=sys.argv[1:3]
# HMMER domtblout: target/accession/tlen query/accession/qlen full-E/full-score/full-bias # of c-E i-E dom-score dom-bias ...
best={}
with open(p,errors='ignore') as f:
  for line in f:
    if not line.strip() or line.startswith('#'): continue
    z=line.split()
    if len(z)<23: continue
    target=z[0]; q=z[3]
    try: full=float(z[6]); ie=float(z[12]); score=float(z[13])
    except: continue
    formal=(full<=1e-4 and ie<=1e-4); exploratory=(full<=1e-2 and ie<=1e-2)
    if not exploratory: continue
    key=q
    cand=(0 if formal else 1, ie, full, -score, target)
    if key not in best or cand<best[key][0]: best[key]=(cand,(target,full,ie,score,formal,exploratory))
with open(out,'w',newline='') as f:
  w=csv.writer(f,delimiter='\t'); w.writerow(['query_id','best_rvt_target','full_evalue','best_domain_i_evalue','domain_score','formal','exploratory'])
  for q in sorted(best):
    target,full,ie,score,formal,exploratory=best[q][1]
    w.writerow([q,target,full,ie,score,str(formal),str(exploratory)])
