from __future__ import annotations
import hashlib
from pathlib import Path

def sha256_file(path):
    p=Path(path)
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20), b''):
            h.update(b)
    return h.hexdigest()

def strict_overlap(a0,a1,b0,b1):
    return int(a0) < int(b1) and int(a1) > int(b0)

def overlap_bp(a0,a1,b0,b1):
    return max(0, min(int(a1),int(b1)) - max(int(a0),int(b0)))

def interval_gap_bp(a0,a1,b0,b1):
    a0,a1,b0,b1=map(int,(a0,a1,b0,b1))
    if strict_overlap(a0,a1,b0,b1):
        return 0
    if a1 <= b0:
        return b0-a1
    return a0-b1

def classify_direct(direct_all):
    if len(direct_all)==0:
        return 'NO_DIRECT_CLEAN_REGION'
    fams=[str(x) for x in direct_all['support_families'].tolist()]
    if any(x=='ACNV_FAMILY;GENERIC_RECURRENCE' for x in fams):
        return 'DUAL_ARM_DIRECT'
    has_a=any(x=='ACNV_FAMILY' for x in fams)
    has_r=any(x=='GENERIC_RECURRENCE' for x in fams)
    if has_a and has_r:
        return 'BOTH_SINGLE_ARMS_DIRECT_NONCONVERGENT'
    if has_a:
        return 'ACNV_ONLY_DIRECT'
    if has_r:
        return 'RECURRENCE_ONLY_DIRECT'
    raise ValueError('unexpected support_families: '+repr(sorted(set(fams))))
