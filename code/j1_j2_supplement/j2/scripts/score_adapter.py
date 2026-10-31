#!/usr/bin/env python3
"""Truth-free bridges to byte-identical project engines; invoked in subprocesses.

Only paths below frozen/ and scores/ are opened. The sealed directory is never
read. Import isolation avoids collision between the historical common.py files.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,json,sys
from pathlib import Path

def read(p):
 with Path(p).open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec);sys.modules[name]=mod;spec.loader.exec_module(mod);return mod
def board(run,cid):
 import pandas as pd
 frozen=run/'frozen';cd=run/'scores'/cid
 mod=load('clean_two_arm_common',frozen/'vendor/clean_two_arm_replay/scripts/common.py')
 q=cd/'arm/outputs'/cid/'STAGE03A_WINDOW_FAMILY_SCORES.tsv.gz'
 df=pd.read_csv(q,sep='\t')
 if set(df.sequence_id.astype(str))!={cid}:raise ValueError('wrong sequence identity in window matrix')
 sel,regs=mod.build_two_arm_board(df,thr=.005,prefix='AR')
 sel.to_csv(cd/'CLEAN_TWO_ARM_SELECTED_PROPOSALS.tsv',sep='\t',index=False)
 regs.to_csv(cd/'CLEAN_TWO_ARM_UNION_REGIONS.tsv',sep='\t',index=False)
 rows=[]
 for r in regs.to_dict('records'):
  rows.append({'rank':int(r['region_rank']),'candidate_id':r['region_id'],'start0':int(r['start0']),
   'end0':int(r['end0']),'strand':'.','support_family_n':int(r['support_family_n'])})
 write(cd/'ARM_M_FULL_FROZEN.tsv',rows,['rank','candidate_id','start0','end0','strand','support_family_n'])
 # Frozen baseline loader requires the full recurrence proposal set, not only
 # those recurrence proposals inside retained dual-arm regions.
 proposals=[]
 for r in sel.to_dict('records'):
  proposals.append({'family':r['family'],'start0':int(r['start0']),'end0':int(r['end0']),
   'qvalue':r['family_q'],'source_path':'frozen-primary-window-matrix'})
 write(cd/'baseline_context/frozen/canonical'/cid/'CANONICAL_PROPOSALS.tsv',proposals,['family','start0','end0','qvalue','source_path'])
 propfile=cd/'baseline_context/frozen/canonical'/cid/'CANONICAL_PROPOSALS.tsv'
 (cd/'RECURRENCE_CONTEXT_READY.json').write_text(json.dumps({'status':'COMPLETE_ORIGINAL_RECURRENCE_PROPOSALS',
   'proposal_sha256':hashlib.sha256(propfile.read_bytes()).hexdigest(),'proposal_n':len(proposals)},indent=2)+'\n')
 print(json.dumps({'status':'FROZEN_TWO_ARM_BOARD_BUILT','opaque_id':cid,'candidate_n':len(rows),'truth_read':False}))

def validate_domtbl(dom,known_queries):
 """Fail on malformed HMM output rather than silently calling it zero hits."""
 lines=dom.read_text().splitlines();footer=False
 for line in lines:
  if line.startswith('# [ok]'):footer=True
  if not line or line.startswith('#'):continue
  fields=line.split()
  if len(fields)<22:raise ValueError('malformed HMMER domain record in '+str(dom))
  if fields[3] not in known_queries:raise ValueError('unmapped HMMER query '+fields[3])
  for i in [6,7,12,13]:float(fields[i])
 if not footer:raise ValueError('HMMER completion marker absent in '+str(dom))

def baseline(run,cid):
 frozen=run/'frozen';cd=run/'scores'/cid;src=frozen/'vendor/formal_baseline/scripts'
 sys.path.insert(0,str(src));sys.modules.pop('common',None)
 mod=load('frozen_formal_score',src/'run_score_freeze.py')
 evaluate=load('frozen_formal_orders',src/'postfreeze_unblind.py')
 rows=read(frozen/'SCORING_MANIFEST.tsv');s=next(r for r in rows if r['opaque_id']==cid)
 pol=json.loads((frozen/'POLICY.json').read_text());tools=json.loads((frozen/'TOOLCHAIN.json').read_text())
 summary=mod.score_case({'case_id':cid},{'chosen_fasta':str(run/s['fasta_path']),'length_bp':s['length_bp']},
  frozen/'db/RVT-All.hmm',tools['prodigal']['path'],tools['hmmscan']['path'],cd/'baseline',pol['cpu'],cd/'baseline_context')
 bdir=cd/'baseline'/cid
 for route,protein in [('prodigal','proteins.faa'),('sixframe','sixframe.faa')]:
  queries={h.split()[0] for h,_ in mod.fasta_records(bdir/route/protein)}
  validate_domtbl(bdir/route/'rvt.domtblout',queries)
 loci=read(bdir/'BASELINE_RT_LOCI.tsv')
 for method in ['PRODIGAL_RVT','SIXFRAME_RVT']:
  rt,rr=evaluate.rt_orders(loci,method)
  for suffix,ordered in [('_ONLY',rt),('_PLUS_FROZEN_REC',rr)]:
   normalized=[]
   for i,r in enumerate(ordered,1):
    normalized.append(dict(r,rank=i,candidate_id=r['locus_id']))
   write(cd/(method+suffix+'.tsv'),normalized,
    ['rank','candidate_id','start0','end0','strand','local_recurrence_hit_n','best_local_recurrence_q','full_evalue','domain_i_evalue','hmm_score'])
 (cd/'BASELINE_RESOURCE_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
 print(json.dumps({'status':'FROZEN_BASELINE_ROUTES_COMPLETE','opaque_id':cid,'truth_read':False,**summary}))

def independent_route(run,cid,route):
 """Reorder orchestration only: one frozen RT route can survive another failure.

Every scientific helper and parameter is imported from the original baseline.
Unlike score_case's sequential all-routes wrapper, this adapter does not require
successful Arm M/other-RT execution to compute an RT-only result. Recurrence-
dependent ranking is emitted only when its genuine frozen proposal file exists.
"""
 frozen=run/'frozen';cd=run/'scores'/cid;src=frozen/'vendor/formal_baseline/scripts'
 sys.path.insert(0,str(src));sys.modules.pop('common',None)
 mod=load('frozen_independent_rt',src/'run_score_freeze.py')
 ev=load('frozen_independent_rt_orders',src/'postfreeze_unblind.py')
 s=next(r for r in read(frozen/'SCORING_MANIFEST.tsv') if r['opaque_id']==cid)
 pol=json.loads((frozen/'POLICY.json').read_text());tool=json.loads((frozen/'TOOLCHAIN.json').read_text())
 recfile=cd/'baseline_context/frozen/canonical'/cid/'CANONICAL_PROPOSALS.tsv'
 recseal=cd/'RECURRENCE_CONTEXT_READY.json';rec_available=False
 if recfile.is_file() and recseal.is_file():
  seal=json.loads(recseal.read_text())
  if hashlib.sha256(recfile.read_bytes()).hexdigest()!=seal['proposal_sha256']:
   raise ValueError('completed recurrence proposal context changed')
  rec_available=True
 rec=mod.load_rec(cd/'baseline_context',cid) if rec_available else []
 fa=run/s['fasta_path'];n=int(s['length_bp']);method='PRODIGAL_RVT' if route=='prodigal' else 'SIXFRAME_RVT'
 dst=cd/'independent_rt'/route;dst.mkdir(parents=True,exist_ok=True)
 faa=dst/('proteins.faa' if route=='prodigal' else 'sixframe.faa');coord={}
 if route=='prodigal':
  mod.run([tool['prodigal']['path'],'-i',str(fa),'-o',str(dst/'genes.gff'),'-a',str(faa),'-f','gff','-p','meta','-q'])
  coord=mod.parse_prodigal_faa_coords(faa)
 else:
  records=list(mod.fasta_records(fa))
  if len(records)!=1:raise ValueError('expected one frozen native contig')
  with faa.open('w') as f:
   for i,(strand,frame,a,b,start,end,pep) in enumerate(mod.sixframe(records[0][1],mod.MINAA),1):
    q=f'SF{i:08d}';coord[q]=(start,end,strand)
    f.write(f'>{q} start0={start} end0={end} strand={strand} frame={frame}\n{pep}\n')
 dom=dst/'rvt.domtblout';tbl=dst/'rvt.tblout'
 mod.run([tool['hmmscan']['path'],'--noali','--cpu',str(pol['cpu']),'-E',str(mod.EMAX),'--domE',str(mod.DOMMAX),
   '--tblout',str(tbl),'--domtblout',str(dom),str(frozen/'db/RVT-All.hmm'),str(faa)])
 validate_domtbl(dom,set(coord));rows=[]
 for q,model,fe,ie,fs,ds in mod.parse_domtblout(dom):
  if q not in coord:raise ValueError('HMM-positive peptide lacks frozen genomic coordinates: '+q)
  start,end,strand=coord[q];nh,rq,rp=mod.local_rec(start,end,rec)
  rows.append({'case_id':cid,'method':method,'candidate_id':q,'start0':start,'end0':end,'strand':strand,
    'hmm_model':model,'full_evalue':fe,'domain_i_evalue':ie,'hmm_score':fs,'domain_score':ds,
    'local_recurrence_hit_n':nh,'best_local_recurrence_q':rq,'recurrence_source':rp})
 fields=['case_id','method','candidate_id','start0','end0','strand','hmm_model','full_evalue','domain_i_evalue','hmm_score','domain_score','local_recurrence_hit_n','best_local_recurrence_q','recurrence_source']
 write(dst/'RAW_RT_CANDIDATES.tsv',rows,fields)
 loci=mod.collapse_rt_loci(rows,method,rec)
 write(dst/'RT_LOCI.tsv',loci,['case_id','method','locus_id','start0','end0','strand','member_n','member_ids','hmm_model','full_evalue','domain_i_evalue','hmm_score','domain_score','local_recurrence_hit_n','best_local_recurrence_q','recurrence_source'])
 rt,rr=ev.rt_orders(loci,method)
 choices=[('_ONLY',rt)]+([('_PLUS_FROZEN_REC',rr)] if rec_available else [])
 for suffix,ordered in choices:
  normalized=[dict(r,rank=i,candidate_id=r['locus_id']) for i,r in enumerate(ordered,1)]
  write(cd/(method+suffix+'.tsv'),normalized,['rank','candidate_id','start0','end0','strand','local_recurrence_hit_n','best_local_recurrence_q','full_evalue','domain_i_evalue','hmm_score'])
 (dst/'ROUTE_STATUS.json').write_text(json.dumps({'route':route,'status':'RT_ONLY_AND_REC_COMPLETE' if rec_available else 'RT_ONLY_COMPLETE_REC_DEPENDENCY_MISSING',
   'rt_locus_n':len(loci),'raw_candidate_n':len(rows),'recurrence_available':rec_available,'truth_read':False},indent=2)+'\n')
 print(json.dumps({'method':method,'rt_locus_n':len(loci),'recurrence_available':rec_available,'truth_read':False}))
def main():
 p=argparse.ArgumentParser();p.add_argument('--run',required=True);p.add_argument('--id',required=True);p.add_argument('--phase',choices=['board','baseline','prodigal','sixframe'],required=True);a=p.parse_args()
 if a.phase in ['prodigal','sixframe']:independent_route(Path(a.run).resolve(),a.id,a.phase)
 else:{'board':board,'baseline':baseline}[a.phase](Path(a.run).resolve(),a.id)
if __name__=='__main__':main()
