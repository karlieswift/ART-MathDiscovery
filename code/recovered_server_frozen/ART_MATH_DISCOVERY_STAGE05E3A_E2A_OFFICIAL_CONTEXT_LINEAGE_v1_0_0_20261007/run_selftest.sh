#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
python - <<'PY'
import json, pandas as pd
from pathlib import Path
import Bio
r=Path('frozen')
fam=pd.read_csv(r/'STAGE05E3A_E2A_FAMILIES.tsv',sep='\t')
can=pd.read_csv(r/'STAGE05E3A_E2A_CANDIDATES.tsv',sep='\t')
wit=pd.read_csv(r/'STAGE05E3A_E2A_WITNESS_PAIRS.tsv',sep='\t')
ids=[x for x in (r/'STAGE05E3A_SEQUENCE_IDS.txt').read_text().splitlines() if x]
assert len(fam)==12
assert len(can)==119 and can.candidate_key.nunique()==119
assert len(ids)==61 and can.sequence_id.nunique()==61
assert len(wit)==80 and wit.candidate_key.nunique()==80 and wit.sequence_id.nunique()==44
assert set(fam.stage05e2_stratum)=={'E2A_EXACT_PRIMARY_JOINT_SIGNATURE'}
print('[PASS] Stage05E3A freeze/policy selftest; Biopython='+Bio.__version__)
PY
python -m py_compile scripts/fetch_genbank.py scripts/analyze_official_context.py
