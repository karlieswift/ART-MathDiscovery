from __future__ import annotations
import hashlib, importlib.util, json, gzip, io, sys
from pathlib import Path
import numpy as np
import pandas as pd

RC_TRANS=str.maketrans('ACGTN','TGCAN')

def sha256_file(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):
            h.update(b)
    return h.hexdigest()

def write_json(path,obj):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    with open(path,'w') as f: json.dump(obj,f,indent=2,sort_keys=True)
    with open(path,'a') as f: f.write('\n')

def read_json(path):
    with open(path) as f: return json.load(f)

def revcomp(s): return s.translate(RC_TRANS)[::-1]

def parse_fasta(path):
    header=None; seq=[]
    with open(path) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if header is not None: raise ValueError('expected exactly one FASTA record')
                header=line[1:].strip().split()[0]
            else: seq.append(line.upper())
    if header is None: raise ValueError('missing FASTA header')
    s=''.join(seq)
    return header,s

def write_fasta(path,header,seq,width=80):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    with open(path,'w') as f:
        f.write('>'+header+'\n')
        for i in range(0,len(seq),width): f.write(seq[i:i+width]+'\n')

def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,str(path))
    if spec is None or spec.loader is None: raise RuntimeError(f'cannot load {path}')
    m=importlib.util.module_from_spec(spec); sys.modules[name]=m; spec.loader.exec_module(m); return m

def robust_score(X):
    X=np.asarray(X,dtype=np.float32)
    if not np.isfinite(X).all(): raise ValueError('nonfinite matrix')
    med=np.median(X,axis=0); dev=np.abs(X-med)
    scale=1.4826*np.median(dev,axis=0); fb=np.mean(dev,axis=0)
    scale=np.where(scale>1e-12,scale,fb); keep=scale>1e-12
    if not np.any(keep): return np.zeros(X.shape[0],dtype=float),0
    z=dev[:,keep]/scale[keep]; np.clip(z,0,20,out=z)
    return np.sqrt(np.mean(z.astype(np.float64)**2,axis=1)),int(keep.sum())

def collapse_oriented_q(df,oriented_score,prefix):
    x=df[['sequence_id','window_len','start0','end0','strand']].copy(); x[prefix+'_score']=oriented_score
    x[prefix+'_oriented_pct']=x.groupby(['window_len','strand'])[prefix+'_score'].rank(method='average',pct=True)
    key=['sequence_id','window_len','start0','end0']
    best=x.loc[x.groupby(key)[prefix+'_oriented_pct'].idxmax(),key+[prefix+'_score',prefix+'_oriented_pct','strand']].rename(columns={prefix+'_score':prefix+'_strandmax_score',prefix+'_oriented_pct':prefix+'_strandmax_pct','strand':prefix+'_best_strand'})
    best[prefix+'_pct']=best.groupby('window_len')[prefix+'_strandmax_pct'].rank(method='average',pct=True)
    nscale=best.groupby('window_len')['start0'].transform('size').astype(float)
    best['q_'+prefix]=np.maximum(0.5/nscale,1-best[prefix+'_pct']+0.5/nscale)
    return best

def empirical_upper_q(values,groups):
    values=np.asarray(values,dtype=float); groups=np.asarray(groups)
    pct=np.zeros(len(values),dtype=float); q=np.zeros(len(values),dtype=float)
    for g in sorted(set(groups.tolist())):
        ii=np.where(groups==g)[0]; s=pd.Series(values[ii])
        p=s.rank(method='average',pct=True).to_numpy(); n=float(len(ii))
        pct[ii]=p; q[ii]=np.maximum(0.5/n,1-p+0.5/n)
    return pct,q

def encode_bases(seq):
    lut=np.full(256,255,dtype=np.uint8)
    for i,b in enumerate(b'ACGT'): lut[b]=i
    x=lut[np.frombuffer(seq.encode('ascii'),dtype=np.uint8)]
    if np.any(x>3): raise ValueError('non-ACGT')
    return x

def rolling_codes(base,k):
    n=len(base)-k+1; out=np.zeros(n,dtype=np.uint32)
    for j in range(k): out=out*4+base[j:j+n].astype(np.uint32)
    return out

def recurrence_scores(seq,windows,k_values=(4,6,8),min_lag=8,max_lag=512,min_pairs=32,min_match_count=2):
    base=encode_bases(seq); N=len(base); nwin=len(windows)
    starts=windows.start0.to_numpy(np.int64); lens=windows.window_len.to_numpy(np.int64); ends=windows.end0.to_numpy(np.int64)
    pref=np.zeros((4,N+1),dtype=np.int32)
    for b in range(4): np.cumsum((base==b).astype(np.int32),out=pref[b,1:])
    cnt=np.column_stack([pref[b,ends]-pref[b,starts] for b in range(4)]).astype(float)
    pbase=np.sum((cnt/lens[:,None])**2,axis=1)
    best_z=np.zeros(nwin); best_k=np.zeros(nwin,dtype=np.int16); best_lag=np.zeros(nwin,dtype=np.int16); best_m=np.zeros(nwin,dtype=np.int32)
    for k in k_values:
        codes=rolling_codes(base,int(k)); nk=lens-int(k)+1; p0=np.power(pbase,int(k))
        d0=max(int(min_lag),int(k)); d1=min(int(max_lag),len(codes)-1)
        for lag in range(d0,d1+1):
            npairs=nk-lag; mask=npairs>=int(min_pairs)
            if not np.any(mask): continue
            m=(codes[:-lag]==codes[lag:]).astype(np.int32); cs=np.empty(len(m)+1,dtype=np.int32); cs[0]=0; np.cumsum(m,dtype=np.int32,out=cs[1:])
            ii=np.where(mask)[0]; st=starts[ii]; npair=npairs[ii]; obs=cs[st+npair]-cs[st]
            pp=np.clip(p0[ii],1e-12,1-1e-12); exp=npair*pp; var=npair*pp*(1-pp)
            z=(obs-exp)/np.sqrt(var+1e-9); z=np.where(obs>=int(min_match_count),np.maximum(z,0.0),0.0)
            upd=z>best_z[ii]
            if np.any(upd):
                jj=ii[upd]; best_z[jj]=z[upd]; best_k[jj]=int(k); best_lag[jj]=int(lag); best_m[jj]=obs[upd]
    return best_z,best_k,best_lag,best_m

def merge_overlapping(selected):
    xs=selected[['start0','end0']].drop_duplicates().sort_values(['start0','end0']).to_numpy(int)
    regs=[]
    for a,b in xs:
        if not regs or int(a)>=int(regs[-1][1]): regs.append([int(a),int(b)])
        else: regs[-1][1]=max(int(regs[-1][1]),int(b))
    return regs

def overlaps(a0,a1,b0,b1): return int(a0)<int(b1) and int(a1)>int(b0)
