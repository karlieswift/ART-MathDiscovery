#!/usr/bin/env python3
"""Resource-bounded replay. Reuses the frozen engine's scientific functions.
One *complete* window-length/strand stratum is computed at a time.
Never estimates medians, scales or ranks per processing chunk.
"""
from pathlib import Path
import argparse,ast,gc,hashlib,importlib.util,json,multiprocessing,os,sys,time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from io_utils import digest,js,read_json,runtime_info,manifest_tree

BLOCKS={'a0':np.arange(0,92),'a4':np.arange(92,348),
        'a5':np.arange(348,1372),'a6':np.arange(1372,5468)}
BOARDS=['STAGE03A_WINDOW_FAMILY_SCORES.tsv.gz',
        'STAGE03A_SELECTED_PROPOSAL_WINDOWS.tsv','STAGE03A_UNION_REGIONS.tsv']

def import_file(name,path):
    spec=importlib.util.spec_from_file_location(name,str(path))
    if spec is None or spec.loader is None:raise RuntimeError('Cannot import '+str(path))
    m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

def load_engine(eng):
    sys.path.insert(0,str(eng/'scripts'))
    common=import_file('common',eng/'scripts/common.py')
    primary=import_file('frozen_stage03a_primary',eng/'scripts/02_primary_scan.py')
    tree=ast.parse((eng/'scripts/02_primary_scan.py').read_text())
    # Read the immutable labels rather than rename or delete an engine branch.
    fun=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='merge_regions')
    fams=None
    for n in fun.body:
        if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='fams' for t in n.targets):
            fams=ast.literal_eval(n.value);break
    if fams is None or [x[1] for x in fams]!=['q_acnv_family','q_irr','q_recurrence']:
        raise RuntimeError('Frozen merge contract differs; refusing approximation')
    return common,primary,fams

def starts_for(n,L,fraction):
    if n<L:return []
    stride=max(1,round(L*fraction));x=list(range(0,n-L+1,stride))
    if x[-1]!=n-L:x.append(n-L)
    return x

def metadata(n,sid,lens,fraction):
    rows=[]
    for L in lens:
        for st in starts_for(n,L,fraction):
            for strand in ('+','-'):rows.append((sid,L,st,st+L,strand))
    return pd.DataFrame(rows,columns=['sequence_id','window_len','start0','end0','strand'])

def engine_signature(eng):
    pairs=[]
    for d in ('scripts','config','frozen_sources'):
        for p in sorted((eng/d).rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc':
                pairs.append((p.relative_to(eng).as_posix(),digest(p)))
    return hashlib.sha256(json.dumps(pairs,separators=(',',':')).encode()).hexdigest()

def checkpoint_signature(eng,fasta):
    return dict(engine_sha256=engine_signature(eng),fasta_file_sha256=digest(fasta),
                implementation_sha256=digest(Path(__file__)),runtime=runtime_info())

def extract_scores(common,primary,eng,seq,df,threads,work,signature):
    score={name:np.zeros(len(df),dtype=float) for name in list(BLOCKS)+['irr']}
    qc={name:[] for name in score}
    work.mkdir(parents=True,exist_ok=True)
    # Explicit fork is supported by the Linux execution target and preserves
    # the imported immutable worker module. No user installation is modified.
    ctx=multiprocessing.get_context('fork')
    for (L,strand),indices0 in df.groupby(['window_len','strand'],sort=True).groups.items():
        ii=np.asarray(list(indices0),dtype=int);n=len(ii)
        token=f'L{L}_{"plus" if strand=="+" else "minus"}'
        cp=work/(token+'.scores.npz');jp=work/(token+'.json')
        sig=dict(signature,window_length=int(L),strand=strand,n_rows=n)
        if cp.exists() and jp.exists():
            old=read_json(jp)
            if old.get('signature')==sig and old.get('scores_sha256')==digest(cp):
                with np.load(cp,allow_pickle=False) as z:
                    for name in score:score[name][ii]=z[name]
                for name in score:qc[name].append(old['qc'][name])
                print(f'[RESUME] {token} n={n}',flush=True);continue
        group={};qgroup={};t=time.time()
        for kind,width in (('acnv',5468),('aux',3570)):
            path=work/(token+'.'+kind+'.npy')
            mm=np.lib.format.open_memmap(path,mode='w+',dtype='float32',shape=(n,width))
            if kind=='acnv':
                init=primary.init_acnv
                args=(str(eng/'frozen_sources/ACNV_SOURCE__legacy_acnv_extractor_frozen.py'),)
                fn=primary.acnv_one
            else:
                init=primary.init_irr
                args=(str(eng/'frozen_sources/IRR_LEGACY_V3__hu_irregular_legacy_v3_features.py'),
                      str(eng/'frozen_sources/DNA_VOCAB_1TO3__vocab_dna_1to3mer.tsv'))
                fn=primary.irr_one
            with ProcessPoolExecutor(max_workers=threads,initializer=init,initargs=args,mp_context=ctx) as ex:
                # Bounded worker submissions; processing batches do not define
                # statistical calibration groups.
                for base in range(0,n,256):
                    glob=ii[base:base+256];items=[]
                    for idx in glob:
                        st=int(df.at[int(idx),'start0']);s=seq[st:st+int(L)]
                        if strand=='-':s=common.revcomp(s)
                        items.append(s if kind=='acnv' else (int(idx),s))
                    for j,value in enumerate(ex.map(fn,items,chunksize=4)):
                        if kind=='aux':
                            idx,x=value
                            if idx!=int(glob[j]):raise RuntimeError('Worker row order changed')
                        else:x=value
                        if np.shape(x)!=(width,) or not np.isfinite(x).all():raise RuntimeError('Invalid feature row')
                        mm[base+j]=x
                    if base==0 or (base//256)%20==0 or base+256>=n:
                        print(f'[{token} {kind}] {min(base+256,n)}/{n} elapsed={time.time()-t:.1f}s',flush=True)
            mm.flush();ri=np.arange(n,dtype=int)
            if kind=='acnv':
                for name,cols in BLOCKS.items():
                    # Same advanced-indexing shape and row order as frozen
                    # score_block; do not switch to feature-column tiles.
                    s,k=common.robust_score(mm[np.ix_(ri,cols)])
                    group[name]=s
                    qgroup[name]=dict(window_len=int(L),strand=strand,n=n,
                                      informative_features=int(k),total_features=len(cols))
            else:
                s,k=common.robust_score(mm[ri]);group['irr']=s
                qgroup['irr']=dict(window_len=int(L),strand=strand,n=n,
                                   informative_features=int(k),total_features=width)
            del mm;gc.collect();path.unlink()
        with cp.open('wb') as f:np.savez(f,**group)
        js(jp,dict(signature=sig,scores_sha256=digest(cp),qc=qgroup))
        for name in score:score[name][ii]=group[name];qc[name].append(qgroup[name])
    return score,qc

def assemble_boards(common,primary,fams,seq,df,score,pol):
    genomic=None
    for name in BLOCKS:
        c=common.collapse_oriented_q(df,score[name],name)
        genomic=c if genomic is None else genomic.merge(c,on=['sequence_id','window_len','start0','end0'],validate='one_to_one')
    other=common.collapse_oriented_q(df,score['irr'],'irr')
    genomic=genomic.merge(other,on=['sequence_id','window_len','start0','end0'],validate='one_to_one')
    arr=genomic[['q_a0','q_a4','q_a5','q_a6']].to_numpy(float)
    surpr=-np.log10(np.clip(arr,1e-300,1));genomic['acnv_family_stat']=surpr.max(axis=1)
    bi=surpr.argmax(axis=1);members=np.asarray(['ACNV92','ACNV4MER_BLOCK','ACNV5MER_BLOCK','ACNV6MER_BLOCK'])
    genomic['acnv_family_best_member']=members[bi];genomic['acnv_family_best_member_q']=arr[np.arange(len(arr)),bi]
    pct,q=common.empirical_upper_q(genomic.acnv_family_stat,genomic.window_len)
    genomic['acnv_family_pct']=pct;genomic['q_acnv_family']=q
    rp=pol['recurrence']
    recstat,bk,bl,bm=common.recurrence_scores(seq,genomic,tuple(rp['k_values']),rp['min_lag_bp'],rp['max_lag_bp'],rp['min_valid_pairs'],rp['min_recurrence_match_count'])
    genomic['recurrence_stat']=recstat;genomic['recurrence_best_k']=bk;genomic['recurrence_best_lag_bp']=bl;genomic['recurrence_best_match_count']=bm
    pct,q=common.empirical_upper_q(recstat,genomic.window_len);genomic['recurrence_pct']=pct;genomic['q_recurrence']=q
    chunks=[];thr=float(pol['proposal_q_threshold'])
    for fam,qc in fams:
        z=genomic.loc[genomic[qc]<=thr,['sequence_id','window_len','start0','end0',qc]].copy()
        z['family']=fam;z=z.rename(columns={qc:'family_q'});chunks.append(z)
    selected=pd.concat(chunks,ignore_index=True).sort_values(['start0','end0','family'])
    regions=primary.merge_regions(selected,genomic) if len(selected) else pd.DataFrame(columns=['region_rank','region_id','start0','end0','span_bp','selected_window_rows','support_family_n','support_families','sum_family_surprisal','best_single_family_q'])
    return genomic,selected,regions

def run(eng,fasta,accession,out,work,threads):
    common,primary,fams=load_engine(eng);pol=common.read_json(eng/'config/ENGINE_POLICY.json')
    header,seq=common.parse_fasta(fasta)
    if not seq or set(seq)-set('ACGT'):raise RuntimeError('Empty or non-ACGT input; no imputation authorized')
    df=metadata(len(seq),header,pol['window_lengths_bp'],pol['stride_fraction'])
    if df.empty:raise RuntimeError('No supported windows; do not fabricate a PASS')
    signature=checkpoint_signature(eng,fasta)
    score,qc=extract_scores(common,primary,eng,seq,df,threads,work,signature)
    out.mkdir(parents=True,exist_ok=True)
    genomic,selected,regions=assemble_boards(common,primary,fams,seq,df,score,pol)
    genomic.to_csv(out/BOARDS[0],sep='\t',index=False,compression={'method':'gzip','mtime':0})
    selected.to_csv(out/BOARDS[1],sep='\t',index=False);regions.to_csv(out/BOARDS[2],sep='\t',index=False)
    gate=dict(stage='Stage04C2-R2C',gate='PASS_COMPLETE_STRATUM_COMPUTATION',gate_pass=True,
              accession=accession,signature=signature,length_bp=len(seq),n_oriented_windows=len(df),
              n_genomic_windows=len(genomic),n_union_regions=len(regions),
              boards={n:digest(out/n) for n in BOARDS},qc=qc,truth_used=False,retuning_authorized=False)
    js(out/'GROUPWISE_GATE.json',gate)
    return gate

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--engine',type=Path,required=True);ap.add_argument('--fasta',type=Path,required=True)
    ap.add_argument('--accession',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--work',type=Path,required=True);ap.add_argument('--threads',type=int,default=8)
    a=ap.parse_args()
    if a.threads<1:raise SystemExit('threads must be positive')
    run(a.engine.resolve(),a.fasta.resolve(),a.accession,a.output.resolve(),a.work.resolve(),a.threads)
