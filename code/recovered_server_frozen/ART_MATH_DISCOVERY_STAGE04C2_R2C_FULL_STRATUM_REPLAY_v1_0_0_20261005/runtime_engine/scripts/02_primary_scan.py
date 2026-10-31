#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,gzip,hashlib,importlib.util,json,math,os,shutil,sys,time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np,pandas as pd
from common import *

_ACNV=None; _IRR=None

def init_acnv(src):
    global _ACNV
    _ACNV=load_module('frozen_acnv_stage03a',src)

def acnv_one(seq):
    n,mu,_=_ACNV.calculate_values(seq); vals=list(n.values())+list(mu.values())
    for k in range(1,7):
        d=_ACNV.calculate_all_NV_Kmer_sliding_window(seq,mu,n,k)
        vals.extend(float(d.get(t,0.0)) for t in _ACNV.generate_kmers(_ACNV.standard_aa,k))
    x=np.asarray(vals,dtype=np.float32)
    if x.shape!=(5468,) or not np.isfinite(x).all(): raise RuntimeError(f'bad ACNV vector {x.shape}')
    return x

def init_irr(src,vocab):
    global _IRR
    sys.argv=['irr','--input','DUMMY.csv','--output','DUMMY.csv','--max-distance','20','--vocab-file',vocab,'--combo-k','2','--combo-include-smaller','--features','V3','--workers','1']
    _IRR=load_module('frozen_irr_stage03a',src)

def irr_one(item):
    idx,seq=item; row=_IRR.compute_irregular_for_sequence(str(seq),'0',idx); x=np.asarray(row[2:],dtype=np.float32)
    if x.shape!=(3570,) or not np.isfinite(x).all(): raise RuntimeError(f'bad irregular vector {x.shape}')
    return idx,x

def make_windows(seq,sid,lens,stride_fraction):
    rows=[]
    for L in lens:
        if len(seq)<L: continue
        stride=max(1,round(L*stride_fraction)); starts=list(range(0,len(seq)-L+1,stride))
        if starts[-1]!=len(seq)-L: starts.append(len(seq)-L)
        for st in starts:
            s=seq[st:st+L]
            for strand,z in [('+',s),('-',revcomp(s))]:
                rows.append((sid,L,st,st+L,strand,hashlib.sha256(z.encode()).hexdigest(),z))
    return pd.DataFrame(rows,columns=['sequence_id','window_len','start0','end0','strand','sha256','sequence'])

def score_block(df,mm,cols,prefix):
    scores=np.zeros(len(df),dtype=float); qc=[]
    for (L,strand),ii0 in df.groupby(['window_len','strand'],sort=True).groups.items():
        ii=np.asarray(list(ii0),dtype=int); s,k=robust_score(mm[np.ix_(ii,np.asarray(cols,dtype=int))]); scores[ii]=s
        qc.append({'window_len':int(L),'strand':strand,'n':len(ii),'informative_features':k,'total_features':len(cols)})
    return collapse_oriented_q(df,scores,prefix),qc

def merge_regions(selected,d):
    fams=[('ACNV_FAMILY','q_acnv_family'),('IRREGULARITY_V3_D20','q_irr'),('GENERIC_RECURRENCE','q_recurrence')]
    regs=merge_overlapping(selected); rows=[]
    for a0,a1 in regs:
        selov=selected[(selected.start0<a1)&(selected.end0>a0)]; allov=d[(d.start0<a1)&(d.end0>a0)]
        support=[]; row={'start0':a0,'end0':a1,'span_bp':a1-a0,'selected_window_rows':int(len(selov))}; qs=[]
        for fam,qc in fams:
            ss=selov[selov.family==fam]
            if len(ss): support.append(fam)
            qbest=float(allov[qc].min()) if len(allov) else 1.0; safe={'ACNV_FAMILY':'acnv_family','IRREGULARITY_V3_D20':'irr','GENERIC_RECURRENCE':'recurrence'}[fam]
            row['best_q_'+safe]=qbest; qs.append(qbest)
        row['support_family_n']=len(support); row['support_families']=';'.join(support); row['sum_family_surprisal']=float(np.sum(-np.log10(np.clip(qs,1e-300,1)))); row['best_single_family_q']=float(min(qs))
        if len(allov):
            ia=allov.q_acnv_family.idxmin(); ir=allov.q_recurrence.idxmin()
            row['acnv_best_member']=str(allov.loc[ia,'acnv_family_best_member']); row['acnv_best_window_start0']=int(allov.loc[ia,'start0']); row['acnv_best_window_end0']=int(allov.loc[ia,'end0']); row['acnv_best_window_len']=int(allov.loc[ia,'window_len'])
            row['recurrence_best_k']=int(allov.loc[ir,'recurrence_best_k']); row['recurrence_best_lag_bp']=int(allov.loc[ir,'recurrence_best_lag_bp']); row['recurrence_best_window_start0']=int(allov.loc[ir,'start0']); row['recurrence_best_window_end0']=int(allov.loc[ir,'end0']); row['recurrence_best_window_len']=int(allov.loc[ir,'window_len'])
        rows.append(row)
    out=pd.DataFrame(rows)
    if len(out):
        out=out.sort_values(['support_family_n','sum_family_surprisal','best_single_family_q','start0'],ascending=[False,False,True,True]).reset_index(drop=True)
        out.insert(0,'region_rank',np.arange(1,len(out)+1)); out.insert(1,'region_id',[f'R{i:04d}' for i in range(1,len(out)+1)])
    return out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); ap.add_argument('--accession',required=True); ap.add_argument('--threads',type=int,default=32); ap.add_argument('--keep-raw',type=int,default=0); a=ap.parse_args()
    root=Path(a.root).resolve(); acc=a.accession; pol=read_json(root/'config/ENGINE_POLICY.json'); gr=root/'work'/acc; out=root/'outputs'/acc; gr.mkdir(parents=True,exist_ok=True); out.mkdir(parents=True,exist_ok=True)
    header,seq=parse_fasta(root/'inputs'/f'{acc}.frozen.fasta')
    if set(seq)-set('ACGT'): raise SystemExit('non-ACGT validation genome; frozen engine has no N policy')
    df=make_windows(seq,header,pol['window_lengths_bp'],pol['stride_fraction']); wp=gr/'windows.tsv.gz'; df.to_csv(wp,sep='\t',index=False,compression={'method':'gzip','mtime':0})
    mm_ac=np.lib.format.open_memmap(gr/'acnv_5468.npy',mode='w+',dtype='float32',shape=(len(df),5468)); t0=time.time(); src=str(root/'frozen_sources/ACNV_SOURCE__legacy_acnv_extractor_frozen.py')
    with ProcessPoolExecutor(max_workers=a.threads,initializer=init_acnv,initargs=(src,)) as ex:
        for i,x in enumerate(ex.map(acnv_one,df.sequence.astype(str).tolist(),chunksize=4)):
            mm_ac[i]=x
            if (i+1)%1000==0: print(f'[{acc} ACNV] {i+1}/{len(df)} elapsed={time.time()-t0:.1f}s',flush=True)
    mm_ac.flush(); ac_sec=time.time()-t0
    mm_ir=np.lib.format.open_memmap(gr/'irr_v3_d20.npy',mode='w+',dtype='float32',shape=(len(df),3570)); t0=time.time(); isrc=str(root/'frozen_sources/IRR_LEGACY_V3__hu_irregular_legacy_v3_features.py'); vocab=str(root/'frozen_sources/DNA_VOCAB_1TO3__vocab_dna_1to3mer.tsv')
    with ProcessPoolExecutor(max_workers=a.threads,initializer=init_irr,initargs=(isrc,vocab)) as ex:
        for idx,x in ex.map(irr_one,enumerate(df.sequence.astype(str).tolist()),chunksize=4):
            mm_ir[idx]=x
            if (idx+1)%500==0: print(f'[{acc} IRR] {idx+1}/{len(df)} elapsed={time.time()-t0:.1f}s',flush=True)
    mm_ir.flush(); ir_sec=time.time()-t0
    # exact frozen proposal arms
    blocks={'a0':np.arange(0,92),'a4':np.arange(92,348),'a5':np.arange(348,1372),'a6':np.arange(1372,5468)}; genomic=None; qc={}
    for name,cols in blocks.items():
        c,q=score_block(df,mm_ac,cols,name); qc[name]=q; genomic=c if genomic is None else genomic.merge(c,on=['sequence_id','window_len','start0','end0'],validate='one_to_one')
    scores=np.zeros(len(df)); iq=[]
    for (L,strand),ii0 in df.groupby(['window_len','strand'],sort=True).groups.items():
        ii=np.asarray(list(ii0),dtype=int); s,k=robust_score(mm_ir[ii]); scores[ii]=s; iq.append({'window_len':int(L),'strand':strand,'n':len(ii),'informative_features':k,'total_features':3570})
    irr=collapse_oriented_q(df,scores,'irr'); genomic=genomic.merge(irr,on=['sequence_id','window_len','start0','end0'],validate='one_to_one')
    # family calibrated max
    qcols=['q_a0','q_a4','q_a5','q_a6']; arr=genomic[qcols].to_numpy(float); surpr=-np.log10(np.clip(arr,1e-300,1)); genomic['acnv_family_stat']=surpr.max(axis=1); bi=surpr.argmax(axis=1); members=np.asarray(['ACNV92','ACNV4MER_BLOCK','ACNV5MER_BLOCK','ACNV6MER_BLOCK']); genomic['acnv_family_best_member']=members[bi]; genomic['acnv_family_best_member_q']=arr[np.arange(len(arr)),bi]
    pct,q=empirical_upper_q(genomic.acnv_family_stat,genomic.window_len); genomic['acnv_family_pct']=pct; genomic['q_acnv_family']=q
    # generic recurrence
    recstat,bk,bl,bm=recurrence_scores(seq,genomic,tuple(pol['recurrence']['k_values']),pol['recurrence']['min_lag_bp'],pol['recurrence']['max_lag_bp'],pol['recurrence']['min_valid_pairs'],pol['recurrence']['min_recurrence_match_count'])
    genomic['recurrence_stat']=recstat; genomic['recurrence_best_k']=bk; genomic['recurrence_best_lag_bp']=bl; genomic['recurrence_best_match_count']=bm; pct,q=empirical_upper_q(recstat,genomic.window_len); genomic['recurrence_pct']=pct; genomic['q_recurrence']=q
    fp=out/'STAGE03A_WINDOW_FAMILY_SCORES.tsv.gz'; genomic.to_csv(fp,sep='\t',index=False,compression={'method':'gzip','mtime':0})
    thr=float(pol['proposal_q_threshold']); fams=[('ACNV_FAMILY','q_acnv_family'),('IRREGULARITY_V3_D20','q_irr'),('GENERIC_RECURRENCE','q_recurrence')]; chunks=[]
    for fam,qcname in fams:
        z=genomic.loc[genomic[qcname]<=thr,['sequence_id','window_len','start0','end0',qcname]].copy(); z['family']=fam; z=z.rename(columns={qcname:'family_q'}); chunks.append(z)
    selected=pd.concat(chunks,ignore_index=True).sort_values(['start0','end0','family']) if chunks else pd.DataFrame()
    sp=out/'STAGE03A_SELECTED_PROPOSAL_WINDOWS.tsv'; selected.to_csv(sp,sep='\t',index=False)
    regions=merge_regions(selected,genomic) if len(selected) else pd.DataFrame(columns=['region_rank','region_id','start0','end0','span_bp','selected_window_rows','support_family_n','support_families','sum_family_surprisal','best_single_family_q'])
    rp=out/'STAGE03A_UNION_REGIONS.tsv'; regions.to_csv(rp,sep='\t',index=False)
    gate={'gate':'PASS_STAGE03A_BLIND_PRIMARY_REGION_BOARD_FROZEN','gate_pass':True,'scientific_status':pol['scientific_status'],'accession':acc,'fasta_header':header,'length_bp':len(seq),'n_oriented_windows':len(df),'n_genomic_windows':len(genomic),'proposal_q_threshold':thr,'proposal_families':[x[0] for x in fams],'n_union_regions':len(regions),'same_window_AND_used':False,'top_k_cutoff_used':False,'all_union_regions_retained':True,'truth_used':False,'postcutoff_ART_family_used':False,'acnv_seconds':ac_sec,'irregular_seconds':ir_sec,'window_family_scores_sha256':sha256_file(fp),'selected_windows_sha256':sha256_file(sp),'union_regions_sha256':sha256_file(rp),'acnv_qc':qc,'irregular_qc':iq}
    write_json(out/'STAGE03A_PRIMARY_GATE.json',gate); print(json.dumps(gate,indent=2))
    if not a.keep_raw:
        for p in [gr/'acnv_5468.npy',gr/'irr_v3_d20.npy']: p.unlink(missing_ok=True)
if __name__=='__main__': main()
