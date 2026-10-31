#!/usr/bin/env python3
import argparse,csv,json
from collections import Counter
from pathlib import Path
ORDER={'NO_HIT':0,'WEAK':1,'SUGGESTIVE':2,'STRONG':3}
def tier(tm,q,t):
 if tm>=.50 and q>=.50 and t>=.50:return 'STRONG'
 if tm>=.40 and q>=.30 and t>=.30:return 'SUGGESTIVE'
 return 'WEAK'
def read_tsv(p):
 with open(p,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(p,rows,cols):
 with open(p,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def classify(p):
 if not p.exists() or p.stat().st_size==0:return 'NO_HIT',None,[]
 hs=[]
 with open(p) as f:
  for line in f:
   z=line.rstrip('\n').split('\t')
   if len(z)<11:continue
   h={'target':z[1],'evalue':float(z[2]),'qcov':float(z[4]),'tcov':float(z[5]),'tm':float(z[10])};h['tier']=tier(h['tm'],h['qcov'],h['tcov']);hs.append(h)
 if not hs:return 'NO_HIT',None,[]
 top=max((h['tier'] for h in hs),key=lambda x:ORDER[x]);same=[h for h in hs if h['tier']==top];rep=max(same,key=lambda h:(h['tm'],h['qcov'],h['tcov'],-h['evalue'],h['target']))
 return top,rep,hs
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--manifest',required=True);ap.add_argument('--foldseek-dir',required=True);ap.add_argument('--out-tsv',required=True);ap.add_argument('--out-gate',required=True);ap.add_argument('--stage',required=True);ap.add_argument('--expected-n',type=int,required=True);a=ap.parse_args()
 m=read_tsv(Path(a.manifest));
 if len(m)!=a.expected_n:raise SystemExit(f'HOLD manifest n {len(m)} expected {a.expected_n}')
 out=[]
 for r in m:
  q=r['query_id'];t,rep,hs=classify(Path(a.foldseek_dir)/f'{q}.m8');z=dict(r);z.update({'tier':t,'representative_target':'' if rep is None else rep['target'],'alntmscore':'' if rep is None else f'{rep["tm"]:.6g}','qcov':'' if rep is None else f'{rep["qcov"]:.6g}','tcov':'' if rep is None else f'{rep["tcov"]:.6g}','evalue':'' if rep is None else f'{rep["evalue"]:.6g}','returned_hit_n':len(hs)});out.append(z)
 write_tsv(Path(a.out_tsv),out,list(out[0]));c=Counter(x['tier'] for x in out)
 gate={'stage':a.stage,'gate':'PASS_'+a.stage.replace('-','_')+'_STRUCTURAL_COMPATIBILITY_COMPLETE','gate_pass':True,'n':len(out),'strong_n':c['STRONG'],'suggestive_n':c['SUGGESTIVE'],'weak_n':c['WEAK'],'no_hit_n':c['NO_HIT'],'positive_n':c['STRONG']+c['SUGGESTIVE'],'classification_rule':'per-hit frozen thresholds then query=max compatibility tier over all returned clean-RT hits','threshold_changed':False,'reference_set_changed':False,'foldseek_search_changed':False,'retuning_performed':False}
 with open(a.out_gate,'w') as f:json.dump(gate,f,indent=2);f.write('\n')
 print(json.dumps(gate,indent=2))
if __name__=='__main__':main()
