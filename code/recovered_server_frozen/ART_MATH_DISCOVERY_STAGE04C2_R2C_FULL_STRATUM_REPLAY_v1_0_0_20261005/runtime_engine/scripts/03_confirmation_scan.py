#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,math,subprocess,sys,time
from functools import lru_cache
from multiprocessing import Pool
from pathlib import Path
import numpy as np,pandas as pd
from common import *

_W={}
SALT='STAGE02B_C_BG_V1'

def bg_key(r):
    s=f"{SALT}|{r.sequence_id}|{int(r.window_len)}|{int(r.start0)}|{int(r.end0)}"
    return hashlib.sha256(s.encode()).hexdigest()

def init_glmy(root):
    root=Path(root); B=load_module('val_glmy_builder',root/'frozen_sources/GLMY_DNA_BUILDER__79A2_build_full_glmy5_fresh_features.py'); C=load_module('val_glmy_core',root/'frozen_sources/GLMY_CORE__glmy5_core.py'); D=C.load_digraph_class(root/'frozen_sources/GLMY_DIGRAPH__Digraph_Hu.py')
    tokens=B.read_vocab(root/'frozen_sources/DNA_VOCAB_1TO2__vocab_dna_1to2mer.tsv'); groups=B.make_token_groups(tokens); lens=sorted({len(x) for x in tokens}); sizes=np.asarray([len(g) for g in groups],dtype=np.int8); dists=(2,3,4,5); maxd=max(dists)
    @lru_cache(maxsize=50000)
    def core_cached(n,edges): return tuple(C.compute_glmy5_from_count(D,n,edges).vector)
    _W.update(B=B,C=C,D=D,tokens=tokens,groups=groups,lens=lens,sizes=sizes,dists=dists,maxd=maxd,core=core_cached)

def glmy_one(task):
    oid,seq=task; B=_W['B']; groups=_W['groups']; dists=_W['dists']; maxd=_W['maxd']; core=_W['core']; pm=B.get_position_maps(seq,_W['lens']); tensor=np.zeros((4,210,5),dtype=np.float32)
    for gi,g in enumerate(groups):
        positions,lens_at=B.group_positions(pm,g); nodes=sorted(positions); n=len(nodes)
        if n==0: continue
        edges=[]
        for ui,u in enumerate(nodes):
            lu=lens_at.get(u,1)
            for vi,v in enumerate(nodes):
                if ui==vi: continue
                lv=lens_at.get(v,1); d=B.adjusted_distance(u,v,lu,lv)
                if d<=maxd: edges.append((ui,vi,int(d)))
        for di,thr in enumerate(dists):
            e=tuple((u,v,w) for u,v,w in edges if w<=thr); tensor[di,gi,:]=core(n,e)
    summ,names=B.summarize_tokens(tensor,_W['sizes']); return oid,summ,names

def fit(X):
    X=np.asarray(X,dtype=float); med=np.nanmedian(X,axis=0); mad=np.nanmedian(np.abs(X-med),axis=0); scale=1.4826*mad; bad=~np.isfinite(scale)|(scale<=1e-12)
    if bad.any():
        meanad=np.nanmean(np.abs(X-med),axis=0); scale[bad]=meanad[bad]
    info=np.isfinite(scale)&(scale>1e-12)&np.isfinite(med); return med,scale,info

def score_rows(X,med,scale,info):
    if info.sum()==0: return np.zeros(len(X)),0
    Z=np.abs((np.asarray(X,dtype=float)[:,info]-med[info])/scale[info]); Z=np.clip(Z,0,20); return np.sqrt(np.mean(Z*Z,axis=1)),int(info.sum())

def window_scores(meta,mat):
    out={}; info_by_len={}
    for L in sorted(meta.window_len.unique()):
        bg_or=meta[(meta.window_len==L)&(meta.panel_role=='background')]; bg_windows=bg_or.panel_id.unique().tolist()
        if len(bg_windows)!=24: raise SystemExit(f'expected 24 backgrounds at length {L}, got {len(bg_windows)}')
        med,sc,info=fit(mat[bg_or.index.to_list()]); info_by_len[int(L)]=int(info.sum())
        for pid in meta.loc[meta.window_len==L,'panel_id'].unique():
            rr=meta[(meta.window_len==L)&(meta.panel_id==pid)]; ridx=rr.index.to_list(); role=rr.iloc[0].panel_role
            if role=='background':
                train=bg_or[bg_or.panel_id!=pid].index.to_list(); m,s,i=fit(mat[train]); scores,ninf=score_rows(mat[ridx],m,s,i)
            else: scores,ninf=score_rows(mat[ridx],med,sc,info)
            best=int(np.argmax(scores)); out[pid]={'score':float(scores[best]),'best_strand':str(rr.iloc[best].strand),'informative_features':ninf}
    return out,info_by_len

def add_q(records,key):
    for L in sorted({int(r['window_len']) for r in records}):
        rr=[r for r in records if int(r['window_len'])==L]; bg=[x[key]['score'] for x in rr if x['panel_role']=='background']
        if len(bg)!=24: raise SystemExit(f'bad background n {L}')
        for x in rr: x[key]['empirical_q']=(1+sum(v>=x[key]['score'] for v in bg))/25.0

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); ap.add_argument('--accession',required=True); ap.add_argument('--workers',type=int,default=16); ap.add_argument('--keep-raw',type=int,default=0); a=ap.parse_args()
    root=Path(a.root).resolve(); acc=a.accession; pol=read_json(root/'config/ENGINE_POLICY.json'); out=root/'outputs'/acc; work=root/'work'/acc; out.mkdir(parents=True,exist_ok=True); work.mkdir(parents=True,exist_ok=True)
    header,seq=parse_fasta(root/'inputs'/f'{acc}.frozen.fasta'); reg=pd.read_csv(out/'STAGE03A_UNION_REGIONS.tsv',sep='\t'); sel=pd.read_csv(out/'STAGE03A_SELECTED_PROPOSAL_WINDOWS.tsv',sep='\t'); uni=pd.read_csv(out/'STAGE03A_WINDOW_FAMILY_SCORES.tsv.gz',sep='\t')
    if len(reg)==0:
        gate={'gate':'PASS_STAGE03A_BLIND_CONFIRMATION_BOARD_FROZEN','gate_pass':True,'accession':acc,'region_n':0,'truth_used':False,'note':'No primary union regions; confirmation board empty.'}; write_json(out/'STAGE03A_CONFIRMATION_GATE.json',gate); pd.DataFrame().to_csv(out/'STAGE03A_CONFIRMATION_BOARD.tsv',sep='\t',index=False); print(json.dumps(gate,indent=2)); return
    anchors=[]
    for _,r in reg.sort_values('region_rank').iterrows():
        x=sel[(sel.start0<r.end0)&(sel.end0>r.start0)].copy()
        if x.empty: raise SystemExit(f'no proposal window overlaps {r.region_id}')
        x=x.sort_values(['family_q','window_len','start0','end0','family']); q=x.iloc[0]
        anchors.append({'panel_id':f"A_{r.region_id}",'panel_role':'anchor','region_id':r.region_id,'region_rank':int(r.region_rank),'region_start0':int(r.start0),'region_end0':int(r.end0),'window_len':int(q.window_len),'start0':int(q.start0),'end0':int(q.end0),'proposal_family':str(q.family),'proposal_family_q':float(q.family_q)})
    anchor_lengths=sorted({x['window_len'] for x in anchors}); regions=[(int(r.start0),int(r.end0)) for _,r in reg.iterrows()]; bgs=[]
    for L in anchor_lengths:
        e=uni[uni.window_len==L].copy(); keep=[]
        for _,r in e.iterrows():
            if any(overlaps(r.start0,r.end0,a0,a1) for a0,a1 in regions): continue
            keep.append(r)
        if len(keep)<24: raise SystemExit(f'not enough nonunion backgrounds length {L}: {len(keep)}')
        keep=sorted(keep,key=bg_key)[:24]
        for j,r in enumerate(keep,1): bgs.append({'panel_id':f'B_L{L}_{j:02d}','panel_role':'background','region_id':'','region_rank':'','region_start0':'','region_end0':'','window_len':int(L),'start0':int(r.start0),'end0':int(r.end0),'proposal_family':'','proposal_family_q':''})
    panel=anchors+bgs; fields=['panel_id','panel_role','region_id','region_rank','region_start0','region_end0','window_len','start0','end0','proposal_family','proposal_family_q']; pp=out/'STAGE03A_ANCHOR_BACKGROUND_PANEL.tsv'
    with open(pp,'w',newline='') as f: w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(panel)
    oriented=[]
    for r in panel:
        frag=seq[int(r['start0']):int(r['end0'])]
        for strand,ss in [('+',frag),('-',revcomp(frag))]:
            z=dict(r); z.update(oriented_id=f"{r['panel_id']}|{strand}",strand=strand,sequence=ss); oriented.append(z)
    op=work/'confirmation_oriented.tsv.gz'
    with gzip.open(op,'wt',newline='') as f: w=csv.DictWriter(f,fieldnames=fields+['oriented_id','strand','sequence'],delimiter='\t'); w.writeheader(); w.writerows(oriented)
    # GLMY
    t0=time.time(); tasks=[(x['oriented_id'],x['sequence']) for x in oriented]
    with Pool(processes=a.workers,initializer=init_glmy,initargs=(str(root),)) as pool: gres=list(pool.imap_unordered(glmy_one,tasks,chunksize=1))
    gres=sorted(gres,key=lambda x:x[0]); ids=[x[0] for x in gres]; local_names=gres[0][2]; gl_names=[f'd{d}_{n}' for d in (2,3,4,5) for n in local_names]; glmat=np.stack([x[1].reshape(-1) for x in gres])
    glp=out/'STAGE03A_GLMY_SUMMARY.tsv.gz'
    with gzip.open(glp,'wt',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['oriented_id']+gl_names)
        for oid,v in zip(ids,glmat): w.writerow([oid]+[f'{float(x):.9g}' for x in v])
    gl_sec=time.time()-t0
    # V3+LGBETA d18
    inp=work/'irr_confirm_input.csv'; raw=work/'irr_confirm_raw.csv'
    with open(inp,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['sequence','label']); w.writeheader();
        for x in oriented: w.writerow({'sequence':x['sequence'],'label':x['oriented_id']})
    cmd=[sys.executable,str(root/'frozen_sources/IRR_LEGACY_V3__hu_irregular_legacy_v3_features.py'),'--input',str(inp),'--output',str(raw),'--vocab-file',str(root/'frozen_sources/DNA_VOCAB_1TO3__vocab_dna_1to3mer.tsv'),'--max-distance','18','--features','V3,LGBETA','--combo-k','2','--combo-include-smaller','--workers',str(a.workers)]
    t0=time.time(); log=root/'logs'/f'{acc}_confirmation_irregular.log'
    with open(log,'w') as lf: p=subprocess.run(cmd,stdout=lf,stderr=subprocess.STDOUT,text=True)
    if p.returncode!=0: raise SystemExit(f'confirmation irregular failed rc={p.returncode}; see {log}')
    ir=pd.read_csv(raw); feat=[c for c in ir.columns if c not in ('sequence','label')]; v3=[c for c in feat if c.endswith('_V3')]; lg=[c for c in feat if c.endswith('_LGBETA')]
    if len(feat)!=7140 or len(v3)!=3570 or len(lg)!=3570: raise SystemExit(f'bad confirmation irregular features {len(feat)} {len(v3)} {len(lg)}')
    ir=ir.drop(columns=['sequence']).rename(columns={'label':'oriented_id'}); irp=out/'STAGE03A_IRR_V3_LGBETA_D18.tsv.gz'; ir.to_csv(irp,sep='\t',index=False,compression={'method':'gzip','mtime':0}); ir_sec=time.time()-t0
    # score
    meta=pd.DataFrame(oriented).drop(columns=['sequence']); base=meta.merge(pd.read_csv(glp,sep='\t'),on='oriented_id',validate='one_to_one').merge(pd.read_csv(irp,sep='\t'),on='oriented_id',validate='one_to_one'); meta=base[['oriented_id','panel_id','panel_role','region_id','region_rank','region_start0','region_end0','window_len','start0','end0','proposal_family','proposal_family_q','strand']].copy(); meta.index=np.arange(len(meta))
    glcols=gl_names; ircols=feat; methods={'glmy_all':glcols,'irr_joint':ircols,'irr_v3_d18':v3,'irr_lgbeta_d18':lg}
    for d in (2,3,4,5): methods[f'glmy_d{d}']=[x for x in glcols if x.startswith(f'd{d}_')]
    first=meta.groupby('panel_id',sort=False).first().reset_index(); records=[]
    for _,r in first.iterrows(): records.append({'panel_id':r.panel_id,'panel_role':r.panel_role,'region_id':r.region_id if pd.notna(r.region_id) else '','region_rank':int(r.region_rank) if str(r.region_rank) not in ('','nan') and pd.notna(r.region_rank) else '','region_start0':int(r.region_start0) if str(r.region_start0) not in ('','nan') and pd.notna(r.region_start0) else '','region_end0':int(r.region_end0) if str(r.region_end0) not in ('','nan') and pd.notna(r.region_end0) else '','window_len':int(r.window_len),'start0':int(r.start0),'end0':int(r.end0),'proposal_family':r.proposal_family if pd.notna(r.proposal_family) else '','proposal_family_q':float(r.proposal_family_q) if str(r.proposal_family_q) not in ('','nan') and pd.notna(r.proposal_family_q) else ''})
    bypid={r['panel_id']:r for r in records}; info={}
    for name,cols in methods.items():
        ws,inf=window_scores(meta,base[cols].to_numpy(dtype=float)); info[name]=inf
        for pid,v in ws.items(): bypid[pid][name]=v
        add_q(records,name)
    flat=[]
    for r in records:
        z={k:r[k] for k in ['panel_id','panel_role','region_id','region_rank','region_start0','region_end0','window_len','start0','end0','proposal_family','proposal_family_q']}
        for name in methods:
            z[f'{name}_score']=r[name]['score']; z[f'{name}_best_strand']=r[name]['best_strand']; z[f'{name}_q']=r[name]['empirical_q']; z[f'{name}_informative_features']=r[name]['informative_features']
        z['orthogonal_consensus']=-math.log10(max(r['glmy_all']['empirical_q']*r['irr_joint']['empirical_q'],1e-300)) if r['panel_role']=='anchor' else np.nan; flat.append(z)
    adf=pd.DataFrame(flat); allp=out/'STAGE03A_ALL_CONFIRMATION_SCORES.tsv'; adf.to_csv(allp,sep='\t',index=False); board=adf[adf.panel_role=='anchor'].copy().sort_values(['orthogonal_consensus','glmy_all_q','irr_joint_q','glmy_all_score','irr_joint_score','region_rank'],ascending=[False,True,True,False,False,True]).reset_index(drop=True); board.insert(0,'confirmation_rank',np.arange(1,len(board)+1)); bp=out/'STAGE03A_CONFIRMATION_BOARD.tsv'; board.to_csv(bp,sep='\t',index=False)
    ok=bool(len(board)==len(reg) and all(v>0 for m in info.values() for v in m.values()))
    gate={'gate':'PASS_STAGE03A_BLIND_CONFIRMATION_BOARD_FROZEN' if ok else 'FAIL_STAGE03A_BLIND_CONFIRMATION_BOARD_FROZEN','gate_pass':ok,'scientific_status':pol['scientific_status'],'accession':acc,'region_n':len(board),'all_union_regions_retained':True,'top_k_cutoff_used':False,'truth_used':False,'postcutoff_ART_family_used':False,'anchor_rule':pol['confirmation']['anchor_rule'],'background_per_window_length':24,'empirical_q_minimum':0.04,'informative_features_by_length':info,'glmy_seconds':gl_sec,'irregular_confirmation_seconds':ir_sec,'panel_sha256':sha256_file(pp),'confirmation_board_sha256':sha256_file(bp),'all_scores_sha256':sha256_file(allp),'authorization':'Blind confirmation board frozen. Do not reveal ART truth until cohort Stage03A final gate passes.'}
    write_json(out/'STAGE03A_CONFIRMATION_GATE.json',gate); print(json.dumps(gate,indent=2))
    if not a.keep_raw:
        for pth in [inp,raw,op]: Path(pth).unlink(missing_ok=True)
    if not ok: raise SystemExit(2)
if __name__=='__main__': main()
