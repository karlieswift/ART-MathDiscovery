#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from common import strict_overlap,overlap_bp,gap_bp,revcomp
assert strict_overlap(10,20,19,30)
assert not strict_overlap(10,20,20,30)
assert overlap_bp(10,20,15,17)==2
assert gap_bp(10,20,25,30)==5
assert revcomp('ACGTN')=='NACGT'
print('[PASS] R4 interval/sequence logic selftest')
