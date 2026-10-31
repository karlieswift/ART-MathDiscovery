#!/usr/bin/env python3
"""Complete historical SA1 equivalence audit; never a new external benchmark."""
from __future__ import annotations
import importlib.util, math
from pathlib import Path
import numpy as np
import pandas as pd
from j2_pipeline import read_tsv,write_tsv,write_json,sha,Hold
from exact_two_arm_scan import load_primary,KEY

def compare(out,ref,vendor):
    cd=out/'scores/HIST_SA1';raw=cd/'arm/outputs/HIST_SA1'
    expected_path=ref/'EXPECTED_FULL_WINDOW_FAMILY_SCORES.tsv.gz'
    if not expected_path.is_file():raise Hold('complete historical SA1 window matrix unavailable; full equivalence remains pending')
    old=pd.read_csv(expected_path,sep='\t');new=pd.read_csv(raw/'STAGE03A_WINDOW_FAMILY_SCORES.tsv.gz',sep='\t')
    m=load_primary(vendor);header,seq=m.parse_fasta(ref/'SA1__MW218148.1.frozen.fasta')
    if set(old.sequence_id.astype(str))!={header}:raise Hold('reference window identity does not match the original native FASTA header')
    if set(new.sequence_id.astype(str))!={'HIST_SA1'}:raise Hold('unexpected recomputed sequence identity')
    old.sequence_id='HIST_SA1';mismatch=[];checks=[]
    keys_equal=old[KEY].reset_index(drop=True).equals(new[KEY].reset_index(drop=True))
    if not keys_equal:mismatch.append({'scope':'all_window_keys','kind':'genomic_key_count_or_order_mismatch','detail':f'{len(old)} vs {len(new)}'})
    columns=[c for c in old.columns if not(c.startswith('irr_') or c=='q_irr')]
    if set(columns)!=set(new.columns):mismatch.append({'scope':'all_window_columns','kind':'column_set_difference','detail':str(set(columns)^set(new.columns))})
    for c in columns:
        if c not in new or len(old)!=len(new):continue
        x=old[c];y=new[c]
        if pd.api.types.is_numeric_dtype(x):
            a=x.to_numpy(float);b=y.to_numpy(float)
            equal=np.isclose(a,b,rtol=1e-10,atol=1e-15,equal_nan=False)
            diff=float(np.max(np.abs(a-b))) if len(a) else 0
            exact=bool(np.array_equal(a,b))
        else:
            equal=x.astype(str).to_numpy()==y.astype(str).to_numpy();diff='';exact=bool(equal.all())
        nbad=int((~equal).sum())
        checks.append({'column':c,'rows_compared':len(x),'all_values_exact':exact,
          'maximum_absolute_difference':diff,'mismatch_rows_at_rtol1e_10_atol1e_15':nbad})
        if nbad:mismatch.append({'scope':'all_window_values','kind':'column_value_difference','detail':c+': '+str(nbad)})
    # Independently regenerate the original complete oriented-key table from the
    # unchanged native DNA and unchanged make_windows, including orientation hash.
    p=m.read_json(vendor/'arm_m/config/ENGINE_POLICY.json')
    expected_oriented=m.make_windows(seq,'HIST_SA1',p['window_lengths_bp'],p['stride_fraction']).drop(columns=['sequence'])
    actual_oriented=pd.read_csv(raw/'ORIENTED_WINDOW_KEYS.tsv.gz',sep='\t') if (raw/'ORIENTED_WINDOW_KEYS.tsv.gz').is_file() else pd.read_csv(cd/'arm/work/HIST_SA1/windows.tsv.gz',sep='\t').drop(columns=['sequence'])
    oriented_equal=expected_oriented.reset_index(drop=True).equals(actual_oriented.reset_index(drop=True))
    if not oriented_equal:mismatch.append({'scope':'all_window_keys','kind':'oriented_key_or_sequence_hash_mismatch','detail':'full native original make_windows comparison failed'})
    spec=importlib.util.spec_from_file_location('historical_clean_two_arm',vendor/'clean_two_arm_replay/scripts/common.py')
    clean=importlib.util.module_from_spec(spec);spec.loader.exec_module(clean)
    expected_sel,expected_board=clean.build_two_arm_board(old,thr=.005,prefix='AR')
    actual_board=pd.read_csv(cd/'CLEAN_TWO_ARM_UNION_REGIONS.tsv',sep='\t')
    board_comparison=[]
    if len(expected_board)!=len(actual_board):mismatch.append({'scope':'full_board','kind':'count_difference','detail':f'{len(expected_board)} vs {len(actual_board)}'})
    for idx in range(min(len(expected_board),len(actual_board))):
        for col in expected_board.columns:
            x=expected_board.iloc[idx][col];y=actual_board.iloc[idx][col]
            isnum=isinstance(x,(int,float,np.number))
            same=math.isclose(float(x),float(y),rel_tol=1e-10,abs_tol=1e-15) if isnum else str(x)==str(y)
            board_comparison.append({'region_rank':idx+1,'column':col,'historical_value':x,'recomputed_value':y,'match':same})
            if not same:mismatch.append({'scope':'full_board','kind':'field_difference','detail':str(idx+1)+' '+col})
    # The independent canonical source is also checked, not merely a board
    # freshly reconstructed with the current wrapper.
    canonical=read_tsv(ref/'EXPECTED_ARM_M_CANONICAL_UNION_REGIONS.tsv')
    canonical_signature=[(int(r['source_rank']),int(r['start0']),int(r['end0']),int(r['has_acnv']),int(r['has_recurrence'])) for r in canonical]
    board_signature=[(int(r.region_rank),int(r.start0),int(r.end0),int('ACNV_FAMILY' in r.support_families),int('GENERIC_RECURRENCE' in r.support_families)) for r in actual_board.itertuples()]
    canonical_equal=canonical_signature==board_signature
    if not canonical_equal:mismatch.append({'scope':'canonical_board','kind':'rank_interval_support_mismatch','detail':'all canonical rows checked'})
    expected_props=read_tsv(ref/'EXPECTED_CANONICAL_PROPOSALS.tsv')
    new_props=read_tsv(cd/'baseline_context/frozen/canonical/HIST_SA1/CANONICAL_PROPOSALS.tsv')
    alias={'ACNV':'ACNV_FAMILY','ACNV_FAMILY':'ACNV_FAMILY','RECURRENCE':'GENERIC_RECURRENCE','GENERIC_RECURRENCE':'GENERIC_RECURRENCE'}
    props=[]
    if len(expected_props)!=len(new_props):mismatch.append({'scope':'canonical_proposals','kind':'count_difference','detail':f'{len(expected_props)} vs {len(new_props)}'})
    for idx,(x,y) in enumerate(zip(expected_props,new_props),1):
        k0=(int(x['start0']),int(x['end0']),alias[x['family']]);k1=(int(y['start0']),int(y['end0']),alias[y['family']])
        same=k0==k1 and math.isclose(float(x['qvalue']),float(y['qvalue']),rel_tol=1e-10,abs_tol=1e-15)
        props.append({'row':idx,'historical_key':str(k0),'recomputed_key':str(k1),'historical_q':x['qvalue'],'recomputed_q':y['qvalue'],'match':same})
        if not same:mismatch.append({'scope':'canonical_proposals','kind':'key_order_or_q_difference','detail':str(idx)})
    write_tsv(out/'ALL_WINDOW_COLUMN_COMPARISON.tsv',checks)
    write_tsv(out/'ALL_BOARD_FIELD_COMPARISON.tsv',board_comparison)
    write_tsv(out/'ALL_CANONICAL_PROPOSAL_COMPARISON.tsv',props)
    result={'status':'COMPLETE_FORMAL_TWO_ARM_HISTORICAL_EQUIVALENCE_MATCH' if not mismatch else 'FORMAL_TWO_ARM_HISTORICAL_EQUIVALENCE_DIFFERENCES',
      'all_native_oriented_windows':len(actual_oriented),'oriented_keys_sequence_hashes_and_order_equal':oriented_equal,
      'all_native_genomic_windows':len(new),'genomic_keys_and_order_equal':keys_equal,
      'non_IRR_window_columns_compared':len(checks),'all_window_cells_compared':sum(r['rows_compared'] for r in checks),
      'numeric_tolerance':{'relative':1e-10,'absolute':1e-15},'window_columns_bitwise_numeric_or_exact_string_equal':sum(r['all_values_exact'] for r in checks),
      'canonical_region_n':len(canonical),'canonical_regions_ranks_coordinates_and_support_equal':canonical_equal,
      'full_board_fields_compared':len(board_comparison),'canonical_proposal_n':len(props),
      'canonical_proposals_matched':sum(r['match'] for r in props),'difference_n':len(mismatch),
      'original_full_window_matrix_sha256':sha(expected_path),'new_external_benchmark':False,
      'opaque_header_mapping':{header:'HIST_SA1'},
      'all_formal_ACNV_and_recurrence_calculations_recomputed_from_complete_native_DNA':True,
      'legacy_IRR_not_a_formal_method_dependency':True}
    write_json(out/'FULL_MATRIX_EQUIVALENCE.json',result)
    return result,mismatch
