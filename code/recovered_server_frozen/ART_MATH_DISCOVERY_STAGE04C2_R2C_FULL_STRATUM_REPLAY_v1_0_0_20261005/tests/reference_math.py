"""Numerical fixtures transcribed from the user's displayed frozen common.py.
These are *test references*, not replacement production scientific functions.
"""
import numpy as np
import pandas as pd

def robust_score(X):
    X=np.asarray(X,dtype=np.float32)
    if not np.isfinite(X).all():raise ValueError('nonfinite matrix')
    med=np.median(X,axis=0);dev=np.abs(X-med)
    scale=1.4826*np.median(dev,axis=0);fb=np.mean(dev,axis=0)
    scale=np.where(scale>1e-12,scale,fb);keep=scale>1e-12
    if not np.any(keep):return np.zeros(X.shape[0],dtype=float),0
    z=dev[:,keep]/scale[keep];np.clip(z,0,20,out=z)
    return np.sqrt(np.mean(z.astype(np.float64)**2,axis=1)),int(keep.sum())

def collapse_oriented_q(df,oriented_score,prefix):
    x=df[['sequence_id','window_len','start0','end0','strand']].copy();x[prefix+'_score']=oriented_score
    x[prefix+'_oriented_pct']=x.groupby(['window_len','strand'])[prefix+'_score'].rank(method='average',pct=True)
    key=['sequence_id','window_len','start0','end0']
    best=x.loc[x.groupby(key)[prefix+'_oriented_pct'].idxmax(),key+[prefix+'_score',prefix+'_oriented_pct','strand']].rename(columns={prefix+'_score':prefix+'_strandmax_score',prefix+'_oriented_pct':prefix+'_strandmax_pct','strand':prefix+'_best_strand'})
    best[prefix+'_pct']=best.groupby('window_len')[prefix+'_strandmax_pct'].rank(method='average',pct=True)
    nscale=best.groupby('window_len')['start0'].transform('size').astype(float)
    best['q_'+prefix]=np.maximum(0.5/nscale,1-best[prefix+'_pct']+0.5/nscale)
    return best

def empirical_upper_q(values,groups):
    values=np.asarray(values,dtype=float);groups=np.asarray(groups)
    pct=np.zeros(len(values),dtype=float);q=np.zeros(len(values),dtype=float)
    for g in sorted(set(groups.tolist())):
        ii=np.where(groups==g)[0];s=pd.Series(values[ii]);p=s.rank(method='average',pct=True).to_numpy();n=float(len(ii))
        pct[ii]=p;q[ii]=np.maximum(0.5/n,1-p+0.5/n)
    return pct,q

def merge_overlapping(selected):
    xs=selected[['start0','end0']].drop_duplicates().sort_values(['start0','end0']).to_numpy(int);regs=[]
    for a,b in xs:
        if not regs or int(a)>=int(regs[-1][1]):regs.append([int(a),int(b)])
        else:regs[-1][1]=max(int(regs[-1][1]),int(b))
    return regs
