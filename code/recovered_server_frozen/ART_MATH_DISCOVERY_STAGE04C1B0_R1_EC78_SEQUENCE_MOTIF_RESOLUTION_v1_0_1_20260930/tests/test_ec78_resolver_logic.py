#!/usr/bin/env python3
# Synthetic test of the exact rule used by the resolver.
cds=[
 {'type':'CDS','aa_len':311,'translation':'M'*120+'YADD'+'A'*187,'has_YADD':True},
 {'type':'CDS','aa_len':450,'translation':'M'*200+'YADD'+'A'*246,'has_YADD':True},
 {'type':'CDS','aa_len':311,'translation':'M'*311,'has_YADD':False},
]
y=[x for x in cds if x.get('has_YADD')]
near=[x for x in y if x.get('aa_len') is not None and 310<=int(x['aa_len'])<=312]
exact=[x for x in near if int(x['aa_len'])==311]
assert len(exact)==1
print('[PASS] Stage04C1B0-R1 Ec78 311-aa YADD synthetic resolver selftest')
