#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
python - <<'PY'
import importlib.util,hashlib,json,tempfile,io,csv
from pathlib import Path
p=Path('scripts/select_logan_metadata.py');spec=importlib.util.spec_from_file_location('sel',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
assert m.classify('marine metagenome','')=='marine'
assert m.classify('viral metagenome','')=='virome'
assert m.classify('soil metagenome','')=='soil'
p=Path('scripts/analyze_frozen_query_expansion.py');spec=importlib.util.spec_from_file_location('ana',p);a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
assert len(a.CODON)==64 and a.CODON['AAG']=='K' and a.CODON['GAG']=='E'
assert a.tr('AAGAAGGAG')=='KKE'
seq='A'*12+'ATG'+'GCT'*50+'TAA'+'C'*12
segs=list(a.six_segments(seq,45));assert any(len(x[-1])>=51 for x in segs)
print('[PASS] Stage05F3A source-selection/translation/ORF selftest')
PY
python - <<'PY'
import json
x=json.load(open('config/F3_MASTER_POLICY.json'))
assert x['expected_query_n']==30 and x['expected_component_n']==26
assert x['source_min_bp']==5000000000 and x['source_max_bp']==10000000000
assert x['frozen_query_fasta_sha256']=='bfe8bfcf249b28e5a13e04cb52af26a58ee18b708c012b2d154a2056e095fac5'
assert x['translation_table_sha256']=='1cd328cab4f124096685dc86bb5dd6badefe8d417c4af97ce68ba789ff560745'
print('[PASS] Stage05F3A frozen policy selftest')
PY
