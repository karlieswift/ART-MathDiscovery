#!/usr/bin/env python3
"""Exact formal ACNV+recurrence orchestration using unmodified frozen helpers.

The legacy IRR feature branch and legacy three-family union are not executed.
Neither is a dependency of the formal clean-two-arm board. This adapter keeps
every native genomic/oriented window and all original calibration operations.
It reads no labels, loci, source identities, or truth files.
"""
from __future__ import annotations
import argparse, importlib.util, json, sys, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd

KEY=['sequence_id','window_len','start0','end0']
MODE='EXACT_FORMAL_TWO_ARM_NO_UNUSED_IRR'

def load_primary(vendor):
    src=vendor/'arm_m/scripts'
    sys.path.insert(0,str(src));sys.modules.pop('common',None)
    name='j2_unmodified_primary_scan'
    spec=importlib.util.spec_from_file_location(name,src/'02_primary_scan.py')
    mod=importlib.util.module_from_spec(spec);sys.modules[name]=mod
    spec.loader.exec_module(mod)
    return mod

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',required=True);p.add_argument('--accession',required=True)
    p.add_argument('--vendor',required=True);p.add_argument('--threads',type=int,default=8)
    p.add_argument('--keep-raw',type=int,choices=[0,1],default=0)
    a=p.parse_args();root=Path(a.root).resolve();vendor=Path(a.vendor).resolve()
    m=load_primary(vendor);pol=m.read_json(root/'config/ENGINE_POLICY.json')
    if a.threads<1:raise ValueError('threads must be positive')
    acc=a.accession;work=root/'work'/acc;out=root/'outputs'/acc
    work.mkdir(parents=True,exist_ok=True);out.mkdir(parents=True,exist_ok=True)
    if (out/'STAGE03A_WINDOW_FAMILY_SCORES.tsv.gz').exists():
        raise ValueError('refusing to overwrite a completed window score matrix')
    header,seq=m.parse_fasta(root/'inputs'/f'{acc}.frozen.fasta')
    if set(seq)-set('ACGT'):raise ValueError('non-ACGT: frozen engine has no N policy')
    df=m.make_windows(seq,header,pol['window_lengths_bp'],pol['stride_fraction'])
    if not len(df):raise ValueError('no original windows')
    df.to_csv(work/'windows.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
    df.drop(columns=['sequence']).to_csv(out/'ORIENTED_WINDOW_KEYS.tsv.gz',sep='\t',index=False,
                                        compression={'method':'gzip','mtime':0})
    mm=np.lib.format.open_memmap(work/'acnv_5468.npy',mode='w+',dtype='float32',shape=(len(df),5468))
    tic=time.monotonic();src=str(root/'frozen_sources/ACNV_SOURCE__legacy_acnv_extractor_frozen.py')
    with ProcessPoolExecutor(max_workers=a.threads,initializer=m.init_acnv,initargs=(src,)) as ex:
        for i,x in enumerate(ex.map(m.acnv_one,df.sequence.astype(str).tolist(),chunksize=4)):
            mm[i]=x
            if (i+1)%1000==0:print(f'[{acc} ACNV exact] {i+1}/{len(df)} elapsed={time.monotonic()-tic:.1f}s',flush=True)
    mm.flush();acsec=time.monotonic()-tic;tic=time.monotonic()
    # Byte-identical block definitions and original score_block implementation.
    blocks={'a0':np.arange(0,92),'a4':np.arange(92,348),'a5':np.arange(348,1372),'a6':np.arange(1372,5468)}
    genomic=None;qc={}
    for name,cols in blocks.items():
        c,q=m.score_block(df,mm,cols,name);qc[name]=q
        genomic=c if genomic is None else genomic.merge(c,on=KEY,validate='one_to_one')
    # Original collapse_oriented_q has exactly one output row per complete key,
    # in sorted groupby-key order. Its legacy IRR merge is one-to-one and does
    # not filter/reorder the already complete ACNV keys. Verify this relation.
    allkeys=df.groupby(KEY,sort=True).size().index.to_frame(index=False).reset_index(drop=True)
    if not genomic[KEY].reset_index(drop=True).equals(allkeys):
        raise ValueError('full genomic key set/order differs from original collapse contract')
    if not (df.groupby(KEY,sort=True).size()==2).all():
        raise ValueError('each genomic window must retain both original orientations')
    # These inline calibration statements reproduce 02_primary_scan.py exactly;
    # they were not separate callable helpers in the frozen implementation.
    qcols=['q_a0','q_a4','q_a5','q_a6'];arr=genomic[qcols].to_numpy(float)
    surpr=-np.log10(np.clip(arr,1e-300,1));genomic['acnv_family_stat']=surpr.max(axis=1)
    bi=surpr.argmax(axis=1);members=np.asarray(['ACNV92','ACNV4MER_BLOCK','ACNV5MER_BLOCK','ACNV6MER_BLOCK'])
    genomic['acnv_family_best_member']=members[bi]
    genomic['acnv_family_best_member_q']=arr[np.arange(len(arr)),bi]
    pct,q=m.empirical_upper_q(genomic.acnv_family_stat,genomic.window_len)
    genomic['acnv_family_pct']=pct;genomic['q_acnv_family']=q
    acscoresec=time.monotonic()-tic;tic=time.monotonic()
    rec=pol['recurrence']
    rs,bk,bl,bm=m.recurrence_scores(seq,genomic,tuple(rec['k_values']),rec['min_lag_bp'],
                                 rec['max_lag_bp'],rec['min_valid_pairs'],rec['min_recurrence_match_count'])
    genomic['recurrence_stat']=rs;genomic['recurrence_best_k']=bk
    genomic['recurrence_best_lag_bp']=bl;genomic['recurrence_best_match_count']=bm
    pct,q=m.empirical_upper_q(rs,genomic.window_len)
    genomic['recurrence_pct']=pct;genomic['q_recurrence']=q
    recsec=time.monotonic()-tic
    fp=out/'STAGE03A_WINDOW_FAMILY_SCORES.tsv.gz'
    genomic.to_csv(fp,sep='\t',index=False,compression={'method':'gzip','mtime':0})
    kp=out/'GENOMIC_WINDOW_KEYS.tsv';genomic[KEY].to_csv(kp,sep='\t',index=False)
    report={'status':'EXACT_FORMAL_TWO_ARM_WINDOW_MATRIX_COMPLETE',
        'execution_mode':MODE,'unused_legacy_IRR_executed':False,
        'legacy_three_family_union_created':False,'formal_board_is_separate_original_clean_two_arm_call':True,
        'all_native_windows_preserved':True,'genomic_keys_equal_original_collapse_order':True,
        'native_length_bp':len(seq),'n_oriented_windows':len(df),'n_genomic_windows':len(genomic),
        'feature_dimensions':5468,'raw_memmap_bytes':int(len(df)*5468*4),
        'acnv_seconds':acsec,'acnv_calibration_seconds':acscoresec,'recurrence_seconds':recsec,
        'window_score_sha256':m.sha256_file(fp),'genomic_key_sha256':m.sha256_file(kp),
        'adapter_sha256':m.sha256_file(Path(__file__).resolve()),
        'original_primary_sha256':m.sha256_file(vendor/'arm_m/scripts/02_primary_scan.py'),
        'original_common_sha256':m.sha256_file(vendor/'arm_m/scripts/common.py'),
        'acnv_qc':qc,'truth_read':False,'new_external_scientific_validation_claimed':False}
    m.write_json(out/'EXACT_TWO_ARM_EXECUTION.json',report)
    print(json.dumps(report,indent=2))
    del mm
    if not a.keep_raw:(work/'acnv_5468.npy').unlink()
    return 0

if __name__=='__main__':sys.exit(main())
