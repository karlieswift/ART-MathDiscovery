#!/usr/bin/env python3
import sys, pandas as pd
from pathlib import Path
D=Path(sys.argv[1])
fl=pd.read_csv(Path(__file__).resolve().parents[1]/'config/FROZEN_LOCI.tsv',sep='\t')
cand=pd.read_csv(D/'config/NOVELTY_FIRST_999_CANDIDATES.tsv',sep='\t')
fan=pd.read_csv(D/'config/NOVELTY_FIRST_PROTEIN_FANOUT.tsv',sep='\t')
links=pd.read_csv(D/'config/NOVELTY_FIRST_CANDIDATE_ORF_LINKS.tsv',sep='\t')
assert set(fl.candidate_key).issubset(set(cand.candidate_key))
assert set(fl.broad_rt_target).issubset(set(fan.target_id))
for ck in fl.candidate_key:
    assert (links.candidate_key==ck).any(), ck
# formal Stage05D1A fixed-threshold output must contain exactly five reported target rows
rows=[x for x in open(D/'outputs/d1a_frozen/RVT_ALL_E1e4.tblout') if x.strip() and not x.startswith('#')]
assert len(rows)==5, f'expected 5 D1A rows, got {len(rows)}'
print('[PASS] Stage05D2 selftest: 3 frozen loci present; D1A fixed-threshold table has 5 rows')
