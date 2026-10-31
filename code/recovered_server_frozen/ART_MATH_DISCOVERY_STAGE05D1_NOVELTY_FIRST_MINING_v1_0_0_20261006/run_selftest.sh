#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
python - "$ROOT" <<'PY'
import sys,pandas as pd,json
from Bio import SeqIO
from pathlib import Path
r=Path(sys.argv[1]); c=pd.read_csv(r/'config/NOVELTY_FIRST_999_CANDIDATES.tsv',sep='\t')
assert len(c)==999 and c.sequence_id.nunique()==402
assert (c.evidence_label=='REPEAT_CONTEXT_RT_NOT_DETECTED_BY_THIS_PANEL').sum()==987
assert (c.evidence_label=='RT_LIKE_AUXILIARY_PROFILE_CONTEXT').sum()==12
p=list(SeqIO.parse(r/'inputs/NOVELTY_FIRST_UNIQUE_PROTEINS.faa','fasta')); assert len(p)==17655
fan=pd.read_csv(r/'config/NOVELTY_FIRST_PROTEIN_FANOUT.tsv',sep='\t'); assert fan.record_key.nunique()==24020
s=json.load(open(r/'config/PREFREEZE_SUMMARY.json')); assert s['n_candidates']==999 and not s['selection_after_broad_RT_results']
print('[PASS] Stage05D1 prefreeze selftest: 999 candidates / 402 genomes / 17,655 unique AA targets')
PY
