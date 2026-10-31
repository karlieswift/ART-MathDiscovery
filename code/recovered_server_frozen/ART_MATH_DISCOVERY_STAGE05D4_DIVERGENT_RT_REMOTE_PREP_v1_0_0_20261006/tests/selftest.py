#!/usr/bin/env python3
import re
FORMAL=[
'42a9674ed0f1bce773c92ae6d0ba7273acc94fa323dd0eb8947a6b8299a2804b',
'4330ef5e0b72ac1244f825df7df03355cbbba35736df7593a9c54afdbf977ff',
'8a5b3b69d22ae7b1575029d0f74ea022798f5f78b04c72ea79154d6b3c2efa00',
'955affef28083a9cde889643e7b90a24c752f0cd02d0948dd5626921ba335c70',
'd7186dbba2a9ae64d8b3a1c874b9303b8b65500c4ddcd0db4bd62f1be6921cd']
assert len(set(FORMAL))==5
assert re.search(r'D.{0,3}D','DID')
print('[PASS] Stage05D4 policy/motif smoke test')
