#!/usr/bin/env python3
from __future__ import annotations
import argparse,hashlib,importlib.util,json,sys,time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from common_stage03a_exact import revcomp, robust_score, collapse_oriented_q, empirical_upper_q, recurrence_scores, merge_overlapping, sha256_file
_ACNV=None
def load_module(name,path):
 spec=importlib.util.spec_from_file_location(name,str(path)); m=importlib.util.module_from_spec(spec); sys.modules[name]=m; spec.loader.exec_module(m); return m
def init_acnv(src):
 global _ACNV; _ACNV=load_module('frozen_acnv_c1',src)
def acnv_one(seq):
 n,mu,_=_ACNV.calculate_values(seq); vals=list(n.values())+list(mu.values())
 for k in range(1,7):
  d=_ACNV.calculate_all_NV_Kmer_sliding_window(seq,mu,n,k); vals.extend(float(d.get(t,0.0)) for t in _ACNV.generate_kmers(_ACNV.standard_aa,k))
 x=np.asarray(vals,dtype=np.float32)
 if x.shape!=(5468,) or not np.isfinite(x).all(): raise RuntimeError(f'bad ACNV vector {x.shape}')
 return x
def parse_fasta(p):
 header=None; parts=[]
 for line in open(p):
  line=line.strip()
  if not line: continue
  if line.startswith('>'):
   if header is not None: raise ValueError('expected exactly one FASTA record')
   header=line[1:].split()[0]
  else: parts.append(line.upper())
 if header is None: raise ValueError('missing FASTA header')
 return header,''.join(parts)
def make_windows(seq,sid,lens,stride_fraction):
 rows=[]
 for L in lens:
  if len(seq)<L: continue
  stride=max(1,round(L*stride_fraction)); starts=list(range(0,len(seq)-L+1,stride))
  if starts[-1]!=len(seq)-L: starts.append(len(seq)-L)
  for st in starts:
   s=seq[st:st+L]
   for strand,z in [('+',s),('-',revcomp(s))]: rows.append((sid,L,st,st+L,strand,hashlib.sha256(z.encode()).hexdigest(),z))
 return pd.DataFrame(rows,columns=['sequence_id','window_len','start0','end0','strand','sha256','sequence'])
def score_block(df,mm,cols,prefix):
 scores=np.zeros(len(df),dtype=float); qc=[]
 for (L,strand),ii0 in df.groupby(['window_len','strand'],sort=True).groups.items():
  ii=np.asarray(list(ii0),dtype=int); s,k=robust_score(mm[np.ix_(ii,np.asarray(cols,dtype=int))]); scores[ii]=s; qc.append({'window_len':int(L),'strand':strand,'n':len(ii),'informative_features':k,'total_features':len(cols)})
 return collapse_oriented_q(df,scores,prefix),qc
def merge_two(selected,d):
 fams=[('ACNV_FAMILY','q_acnv_family'),('GENERIC_RECURRENCE','q_recurrence')]; regs=merge_overlapping(selected); rows=[]
 for a0,a1 in regs:
  selov=selected[(selected.start0<a1)&(selected.end0>a0)]; allov=d[(d.start0<a1)&(d.end0>a0)]; support=[]; qs=[]; row={'start0':a0,'end0':a1,'span_bp':a1-a0,'selected_window_rows':int(len(selov))}
  for fam,qc in fams:
   if len(selov[selov.family==fam]): support.append(fam)
   qbest=float(allov[qc].min()) if len(allov) else 1.0; safe='acnv_family' if fam=='ACNV_FAMILY' else 'recurrence'; row['best_q_'+safe]=qbest; qs.append(qbest)
  row['support_family_n']=len(support); row['support_families']=';'.join(support); row['sum_two_family_surprisal']=float(np.sum(-np.log10(np.clip(qs,1e-300,1)))); row['best_single_family_q']=float(min(qs))
  if len(allov):
   ia=allov.q_acnv_family.idxmin(); ir=allov.q_recurrence.idxmin(); row['acnv_best_member']=str(allov.loc[ia,'acnv_family_best_member']); row['acnv_best_window_start0']=int(allov.loc[ia,'start0']); row['acnv_best_window_end0']=int(allov.loc[ia,'end0']); row['acnv_best_window_len']=int(allov.loc[ia,'window_len']); row['recurrence_best_k']=int(allov.loc[ir,'recurrence_best_k']); row['recurrence_best_lag_bp']=int(allov.loc[ir,'recurrence_best_lag_bp']); row['recurrence_best_window_start0']=int(allov.loc[ir,'start0']); row['recurrence_best_window_end0']=int(allov.loc[ir,'end0']); row['recurrence_best_window_len']=int(allov.loc[ir,'window_len'])
  rows.append(row)
 out=pd.DataFrame(rows)
 if len(out):
  out=out.sort_values(['support_family_n','sum_two_family_surprisal','best_single_family_q','start0'],ascending=[False,False,True,True]).reset_index(drop=True); out.insert(0,'region_rank',np.arange(1,len(out)+1)); out.insert(1,'region_id',[f'AR{i:04d}' for i in range(1,len(out)+1)])
 return out
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); ap.add_argument('--control-id',required=True); ap.add_argument('--accession',required=True); ap.add_argument('--fasta',required=True); ap.add_argument('--threads',type=int,default=32); ap.add_argument('--keep-raw',type=int,default=0); a=ap.parse_args()
 root=Path(a.root); pol=json.loads((root/'config/C1_POLICY.json').read_text()); out=root/'outputs'/a.control_id; work=root/'work'/a.control_id; out.mkdir(parents=True,exist_ok=True); work.mkdir(parents=True,exist_ok=True)
 header,seq=parse_fasta(a.fasta)
 if set(seq)-set('ACGT'): raise SystemExit('non-ACGT genome; frozen engine has no N policy')
 df=make_windows(seq,header,pol['window_lengths_bp'],pol['stride_fraction']); print(f'[{a.control_id}] length={len(seq)} oriented_windows={len(df)}',flush=True)
 mm=np.lib.format.open_memmap(work/'acnv_5468.npy',mode='w+',dtype='float32',shape=(len(df),5468)); src=str(root/'frozen_sources/ACNV_SOURCE__legacy_acnv_extractor_frozen.py'); t0=time.time()
 with ProcessPoolExecutor(max_workers=a.threads,initializer=init_acnv,initargs=(src,)) as ex:
  for i,x in enumerate(ex.map(acnv_one,df.sequence.astype(str).tolist(),chunksize=4)):
   mm[i]=x
   if (i+1)%1000==0: print(f'[{a.control_id} ACNV] {i+1}/{len(df)} elapsed={time.time()-t0:.1f}s',flush=True)
 mm.flush(); ac_sec=time.time()-t0
 blocks={'a0':np.arange(0,92),'a4':np.arange(92,348),'a5':np.arange(348,1372),'a6':np.arange(1372,5468)}; genomic=None; qc={}
 for name,cols in blocks.items():
  c,q=score_block(df,mm,cols,name); qc[name]=q; genomic=c if genomic is None else genomic.merge(c,on=['sequence_id','window_len','start0','end0'],validate='one_to_one')
 qcols=['q_a0','q_a4','q_a5','q_a6']; arr=genomic[qcols].to_numpy(float); surpr=-np.log10(np.clip(arr,1e-300,1)); genomic['acnv_family_stat']=surpr.max(axis=1); bi=surpr.argmax(axis=1); members=np.asarray(['ACNV92','ACNV4MER_BLOCK','ACNV5MER_BLOCK','ACNV6MER_BLOCK']); genomic['acnv_family_best_member']=members[bi]; genomic['acnv_family_best_member_q']=arr[np.arange(len(arr)),bi]; pct,q=empirical_upper_q(genomic.acnv_family_stat,genomic.window_len); genomic['acnv_family_pct']=pct; genomic['q_acnv_family']=q
 rec=pol['recurrence']; recstat,bk,bl,bm=recurrence_scores(seq,genomic,tuple(rec['k_values']),rec['min_lag_bp'],rec['max_lag_bp'],rec['min_valid_pairs'],rec['min_recurrence_match_count']); genomic['recurrence_stat']=recstat; genomic['recurrence_best_k']=bk; genomic['recurrence_best_lag_bp']=bl; genomic['recurrence_best_match_count']=bm; pct,q=empirical_upper_q(recstat,genomic.window_len); genomic['recurrence_pct']=pct; genomic['q_recurrence']=q
 fp=out/'C1_WINDOW_FAMILY_SCORES.tsv.gz'; genomic.to_csv(fp,sep='\t',index=False,compression={'method':'gzip','mtime':0})
 thr=float(pol['proposal_q_threshold']); chunks=[]
 for fam,qcname in [('ACNV_FAMILY','q_acnv_family'),('GENERIC_RECURRENCE','q_recurrence')]:
  z=genomic.loc[genomic[qcname]<=thr,['sequence_id','window_len','start0','end0',qcname]].copy(); z['family']=fam; z=z.rename(columns={qcname:'family_q'}); chunks.append(z)
 selected=pd.concat(chunks,ignore_index=True).sort_values(['start0','end0','family']); sp=out/'C1_SELECTED_PROPOSALS.tsv'; selected.to_csv(sp,sep='\t',index=False); regions=merge_two(selected,genomic); rp=out/'C1_UNION_REGIONS.tsv'; regions.to_csv(rp,sep='\t',index=False)
 gate={'gate':'PASS_STAGE04C1_TWO_ARM_SCORE_BOARD_FROZEN','gate_pass':True,'scientific_status':pol['scientific_status'],'control_id':a.control_id,'accession':a.accession,'length_bp':len(seq),'n_oriented_windows':len(df),'n_genomic_windows':len(genomic),'proposal_q_threshold':thr,'proposal_families':['ACNV_FAMILY','GENERIC_RECURRENCE'],'n_union_regions':len(regions),'truth_used':False,'annotation_used':False,'rt_coordinates_used':False,'partner_coordinates_used':False,'retuning_authorized':False,'acnv_seconds':ac_sec,'window_scores_sha256':sha256_file(fp),'selected_sha256':sha256_file(sp),'union_regions_sha256':sha256_file(rp),'acnv_qc':qc}
 (out/'C1_SCORE_FREEZE_GATE.json').write_text(json.dumps(gate,indent=2)+'\n'); print(json.dumps(gate,indent=2))
 if not a.keep_raw: (work/'acnv_5468.npy').unlink(missing_ok=True)
if __name__=='__main__': main()
