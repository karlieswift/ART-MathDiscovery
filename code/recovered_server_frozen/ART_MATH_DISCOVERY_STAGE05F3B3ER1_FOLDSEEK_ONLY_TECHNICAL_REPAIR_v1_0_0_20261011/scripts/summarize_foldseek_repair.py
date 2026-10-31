#!/usr/bin/env python3
import csv,hashlib,json,re,sys
from pathlib import Path
prior=Path(sys.argv[1]); fstatp=Path(sys.argv[2]); folddir=Path(sys.argv[3]); lockp=Path(sys.argv[4]); outdir=Path(sys.argv[5]); outdir.mkdir(parents=True,exist_ok=True)
def h(p):
 x=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): x.update(b)
 return x.hexdigest()
def ff(x):
 try:return float(x)
 except:return float('nan')
def normacc(t):
 m=re.search(r'AF-([A-Za-z0-9]+)-F\d+-model',t); return m.group(1) if m else t
def tier(tm,qc,tc):
 if tm>=.50 and qc>=.50 and tc>=.50:return 'STRONG'
 if tm>=.40 and qc>=.30 and tc>=.30:return 'SUGGESTIVE'
 return 'WEAK'
rows=list(csv.DictReader(open(prior),delimiter='\t')); assert len(rows)==32
fs={r['candidate_id']:r for r in csv.DictReader(open(fstatp),delimiter='\t')}
targets=[r for r in rows if r['f3b3e_final_best_tier']=='FOLDSEEK_TECHNICAL_FAILURE_AFTER_SERIAL_RESOURCE_CLOSURE']
assert len(targets)==5 and set(fs)=={r['candidate_id'] for r in targets}
fields=['query','target','evalue','bits','qcov','tcov','alnlen','qlen','tlen','fident','alntmscore']
out=[]; hitsall=[]
for p in rows:
 r=dict(p); cid=r['candidate_id']; r['b3er1_resolution_source']='PRE_B3ER1'; r['b3er1_final_best_tier']=r['f3b3e_final_best_tier']
 if r['f3b3e_final_best_tier']!='FOLDSEEK_TECHNICAL_FAILURE_AFTER_SERIAL_RESOURCE_CLOSURE': out.append(r); continue
 r['b3er1_resolution_source']='FOLDSEEK_ONLY_TECHNICAL_REPAIR'; s=fs[cid]
 if s['foldseek_status']!='PASS':
  r['foldseek_status']='TECHNICAL_FAILURE'; r['b3er1_final_best_tier']='FOLDSEEK_TECHNICAL_FAILURE_AFTER_DIRECT_REPAIR'; out.append(r); continue
 r['foldseek_status']='PASS'; m8=folddir/f'{cid}.m8'; hs=[]
 if m8.exists() and m8.stat().st_size:
  for x in csv.DictReader(open(m8),delimiter='\t',fieldnames=fields):
   x['candidate_id']=cid; x['target_accession']=normacc(x['target']); x['tier']=tier(ff(x['alntmscore']),ff(x['qcov']),ff(x['tcov'])); hs.append(x); hitsall.append(x)
 r['r5c_hit_n']=str(len(hs))
 if not hs:
  r['best_clean_rt_target']='';r['best_alntmscore']='';r['best_qcov']='';r['best_tcov']='';r['b3er1_final_best_tier']='NO_HIT'
 else:
  b=max(hs,key=lambda x:ff(x['alntmscore'])); r['best_clean_rt_target']=b['target_accession'];r['best_alntmscore']=b['alntmscore'];r['best_qcov']=b['qcov'];r['best_tcov']=b['tcov'];r['b3er1_final_best_tier']=b['tier']
 out.append(r)
base=['shortlist_rank','family_rep_id','candidate_id','accession','aa_length','known_sequence_reference_hit_n','known_top_reference_id','known_top_reference_original_headers','f3b3e_resolution_source','pre_f3b3e_best_tier','serial_exclusive_boltz_status','serial_failure_category','serial_cif_path','foldseek_status','r5c_hit_n','best_clean_rt_target','best_alntmscore','best_qcov','best_tcov','f3b3e_final_best_tier']
cols=base+['b3er1_resolution_source','b3er1_final_best_tier']
final=outdir/'F3B3ER1_FINAL_STRUCTURE_R5C_SUMMARY_32.tsv'
with final.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n',extrasaction='ignore'); w.writeheader();w.writerows(out)
hitcols=['candidate_id']+fields+['target_accession','tier']
with (outdir/'F3B3ER1_REPAIRED_CLEAN_RT_HITS.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=hitcols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(hitsall)
keys=['STRONG','SUGGESTIVE','WEAK','NO_HIT','RESOURCE_LIMITED_BOLTZ_UNRESOLVED','FOLDSEEK_TECHNICAL_FAILURE_AFTER_DIRECT_REPAIR']
c={k:sum(r['b3er1_final_best_tier']==k for r in out) for k in keys}
unres=c['RESOURCE_LIMITED_BOLTZ_UNRESOLVED']+c['FOLDSEEK_TECHNICAL_FAILURE_AFTER_DIRECT_REPAIR']
gate={
 'stage':'Stage05F3B3E-R1','gate':'PASS_F3B3ER1_FOLDSEEK_ONLY_TECHNICAL_REPAIR_COMPLETE' if c['FOLDSEEK_TECHNICAL_FAILURE_AFTER_DIRECT_REPAIR']==0 else 'PASS_F3B3ER1_FOLDSEEK_ONLY_REPAIR_ATTEMPT_COMPLETE_WITH_TECHNICAL_FAILURES','gate_pass':True,
 'structure_handoff_family_n':32,'foldseek_repair_target_n':5,'foldseek_repair_pass_n':5-c['FOLDSEEK_TECHNICAL_FAILURE_AFTER_DIRECT_REPAIR'],'foldseek_repair_failure_n':c['FOLDSEEK_TECHNICAL_FAILURE_AFTER_DIRECT_REPAIR'],
 'final_resolved_structure_n':32-unres,'final_r5c_strong_n':c['STRONG'],'final_r5c_suggestive_n':c['SUGGESTIVE'],'final_r5c_weak_n':c['WEAK'],'final_r5c_no_hit_n':c['NO_HIT'],'final_strong_or_suggestive_n':c['STRONG']+c['SUGGESTIVE'],
 'resource_limited_boltz_unresolved_n':c['RESOURCE_LIMITED_BOLTZ_UNRESOLVED'],'foldseek_technical_failure_n':c['FOLDSEEK_TECHNICAL_FAILURE_AFTER_DIRECT_REPAIR'],'technical_structure_closure_complete':unres==0,
 'boltz_recomputed':False,'scientific_threshold_changed':False,'shortlist_changed':False,'sequence_changed':False,
 'implementation_lock_sha256':h(lockp),'final_summary_sha256':h(final),
 'claim_boundary':'Foldseek runtime/handoff repair only; frozen R5C thresholds unchanged. Residual failures remain technical/resource-censored and are not structural negatives.'}
with (outdir/'STAGE05F3B3ER1_GATE.json').open('w') as f: json.dump(gate,f,indent=2);f.write('\n')
print(json.dumps(gate,indent=2))
