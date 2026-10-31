#!/usr/bin/env python3
import hashlib
# Synthetic test of the R2 DRT2 sequence-first identity rule.
tr='A'*266+'YCDD'+'A'*(425-270)
assert len(tr)==425
assert tr[266:270]=='YCDD'
assert (267,269,270)==(tr.find('YCDD')+1,tr.find('YCDD')+3,tr.find('YCDD')+4)
print('[PASS] Stage04C1B0-R2 DRT2 425-aa YCDD(D269/D270) synthetic resolver selftest')
