from __future__ import annotations
import csv,gzip,hashlib,io,json,tarfile
from pathlib import Path
import numpy as np
import pandas as pd

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def read_json(p): return json.loads(Path(p).read_text())
def write_json(p,obj): Path(p).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
def overlap(a0,a1,b0,b1): return a0 < b1 and a1 > b0

def merge_intervals(selected):
    a=selected[['start0','end0']].drop_duplicates().sort_values(['start0','end0']).to_numpy(int)
    out=[]
    for x,y in a:
        if not out or x >= out[-1][1]: out.append([int(x),int(y)])
        else: out[-1][1]=max(out[-1][1],int(y))
    return out

def validate_scores(d):
    req=['sequence_id','window_len','start0','end0','q_acnv_family','q_recurrence',
         'acnv_family_best_member','recurrence_best_k','recurrence_best_lag_bp']
    miss=[x for x in req if x not in d.columns]
    if miss: raise SystemExit(f'FAIL missing required frozen-score columns: {miss}')
    if len(d)==0: raise SystemExit('FAIL empty frozen score table')

def build_two_arm_board(d,thr=0.005,prefix='TWO_ARM'):
    validate_scores(d)
    chunks=[]
    for fam,qc in [('ACNV_FAMILY','q_acnv_family'),('GENERIC_RECURRENCE','q_recurrence')]:
        z=d.loc[d[qc]<=thr,['sequence_id','window_len','start0','end0',qc]].copy()
        z['family']=fam; z=z.rename(columns={qc:'family_q'}); chunks.append(z)
    sel=pd.concat(chunks,ignore_index=True).sort_values(['start0','end0','family']).reset_index(drop=True)
    rows=[]
    for a0,a1 in merge_intervals(sel):
        so=sel[(sel.start0<a1)&(sel.end0>a0)]
        ao=d[(d.start0<a1)&(d.end0>a0)]
        qa=float(ao.q_acnv_family.min()); qr=float(ao.q_recurrence.min()); fams=[]
        if (so.family=='ACNV_FAMILY').any(): fams.append('ACNV_FAMILY')
        if (so.family=='GENERIC_RECURRENCE').any(): fams.append('GENERIC_RECURRENCE')
        ia=ao.q_acnv_family.idxmin(); ir=ao.q_recurrence.idxmin()
        rows.append({
            'start0':int(a0),'end0':int(a1),'span_bp':int(a1-a0),
            'selected_window_rows':int(len(so)),
            'support_family_n':len(fams),'support_families':';'.join(fams),
            'best_q_acnv_family':qa,'best_q_recurrence':qr,
            'sum_two_family_surprisal':float(-np.log10(max(qa,1e-300))-np.log10(max(qr,1e-300))),
            'best_single_family_q':min(qa,qr),
            'acnv_best_member':str(ao.loc[ia,'acnv_family_best_member']),
            'acnv_best_window_start0':int(ao.loc[ia,'start0']),
            'acnv_best_window_end0':int(ao.loc[ia,'end0']),
            'acnv_best_window_len':int(ao.loc[ia,'window_len']),
            'recurrence_best_k':int(float(ao.loc[ir,'recurrence_best_k'])),
            'recurrence_best_lag_bp':int(float(ao.loc[ir,'recurrence_best_lag_bp'])),
            'recurrence_best_window_start0':int(ao.loc[ir,'start0']),
            'recurrence_best_window_end0':int(ao.loc[ir,'end0']),
            'recurrence_best_window_len':int(ao.loc[ir,'window_len'])
        })
    r=pd.DataFrame(rows)
    if len(r):
        r=r.sort_values(['support_family_n','sum_two_family_surprisal','best_single_family_q','start0'],
                        ascending=[False,False,True,True]).reset_index(drop=True)
        r.insert(0,'region_rank',np.arange(1,len(r)+1))
        r.insert(1,'region_id',[f'{prefix}{i:04d}' for i in range(1,len(r)+1)])
    else:
        r=pd.DataFrame(columns=['region_rank','region_id','start0','end0','span_bp','selected_window_rows',
                                'support_family_n','support_families','best_q_acnv_family','best_q_recurrence',
                                'sum_two_family_surprisal','best_single_family_q'])
    return sel,r

def best_truth_window(d,s0,e0,col,extra=None):
    ws=d[(d.start0<e0)&(d.end0>s0)]
    if len(ws)==0: return None
    b=ws.sort_values([col,'window_len','start0','end0'],ascending=[True,True,True,True]).iloc[0]
    q=float(b[col]); allq=d[col].astype(float).to_numpy(); rank=1+int(np.sum(allq<q))
    out={'q':q,'global_rank':rank,'n_windows':int(len(allq)),'start0':int(b.start0),'end0':int(b.end0),'window_len':int(b.window_len)}
    if extra:
        for dst,src,kind in extra:
            v=b[src]
            out[dst]=int(float(v)) if kind=='int' else str(v)
    return out

def read_tsv_from_tar(tf,suffix):
    hits=[n for n in tf.getnames() if n.endswith(suffix)]
    if len(hits)!=1: raise SystemExit(f'FAIL tar member {suffix}: {hits}')
    raw=tf.extractfile(hits[0]).read()
    if suffix.endswith('.gz'):
        raw=gzip.decompress(raw)
    return pd.read_csv(io.BytesIO(raw),sep='\t'),hits[0]
