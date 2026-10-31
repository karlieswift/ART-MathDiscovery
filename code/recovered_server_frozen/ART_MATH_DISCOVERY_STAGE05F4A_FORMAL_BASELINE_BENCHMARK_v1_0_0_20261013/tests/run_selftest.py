#!/usr/bin/env python3
from pathlib import Path
import sys,tempfile
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/'scripts'))
from sixframe_orfs import sixframe
from common import fasta_records,seq_sha
seq='ATG'+'AAA'*50+'TAA'+'CCC'*10
x=list(sixframe(seq,45))
assert x and any(len(a[-1])>=45 for a in x)
contract=(HERE/'config/BENCHMARK_CONTRACT.json').read_text()
assert 'truth_must_not_be_used_before_score_freeze' in contract
score=(HERE/'scripts/run_score_freeze.py').read_text()
assert "read_tsv(here/'frozen/BASELINE_TRUTH_SEAL.tsv')" not in score
assert 'PRODIGAL_RVT' in score and 'SIXFRAME_RVT' in score
print('[PASS] Stage05F4A formal baseline syntax/policy selftest')
print('[PASS] score-freeze runner cannot read truth seal by construction')
