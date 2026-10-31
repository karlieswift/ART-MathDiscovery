#!/usr/bin/env python3
import csv, sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader(open(root/'config/FROZEN_REMAINING_BROAD_RT_TARGETS.tsv'),delimiter='\t'))
assert len(rows)==3
assert {r['accession'] for r in rows}=={'MT708547.1','NC_031049.1','PV251870.1'}
if len(sys.argv)>1 and sys.argv[1].strip():
    d1=Path(sys.argv[1])
    need=[
      d1/'config/NOVELTY_FIRST_24020_ORF_REGISTRY.tsv',
      d1/'config/NOVELTY_FIRST_PROTEIN_FANOUT.tsv',
      d1/'inputs/NOVELTY_FIRST_UNIQUE_PROTEINS.faa',
      d1/'outputs/d1a_frozen/RVT_ALL_E1e4.tblout',
      d1/'outputs/d1a_frozen/RVT_ALL_E1e4.domtblout',
    ]
    miss=[str(p) for p in need if not p.exists()]
    assert not miss, 'missing D1 inputs: '+', '.join(miss)
print('[PASS] Stage05D3 selftest: 3 remaining formal broad-RVT targets frozen; no DNA rescoring')
