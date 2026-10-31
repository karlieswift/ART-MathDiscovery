#!/usr/bin/env python3
from pathlib import Path
import sys
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from common import strict_overlap, overlap_bp, interval_gap_bp, classify_direct

assert strict_overlap(10,20,19,30)
assert not strict_overlap(10,20,20,30)
assert overlap_bp(10,20,15,25)==5
assert overlap_bp(10,20,20,25)==0
assert interval_gap_bp(10,20,20,25)==0
assert interval_gap_bp(10,20,23,25)==3
assert interval_gap_bp(23,25,10,20)==3

def df(fams):
    return pd.DataFrame({'support_families':fams})
assert classify_direct(df([]))=='NO_DIRECT_CLEAN_REGION'
assert classify_direct(df(['ACNV_FAMILY']))=='ACNV_ONLY_DIRECT'
assert classify_direct(df(['GENERIC_RECURRENCE']))=='RECURRENCE_ONLY_DIRECT'
assert classify_direct(df(['ACNV_FAMILY','GENERIC_RECURRENCE']))=='BOTH_SINGLE_ARMS_DIRECT_NONCONVERGENT'
assert classify_direct(df(['ACNV_FAMILY;GENERIC_RECURRENCE']))=='DUAL_ARM_DIRECT'
assert classify_direct(df(['ACNV_FAMILY','ACNV_FAMILY;GENERIC_RECURRENCE']))=='DUAL_ARM_DIRECT'
print('[PASS] interval and cassette-level classification semantics')
