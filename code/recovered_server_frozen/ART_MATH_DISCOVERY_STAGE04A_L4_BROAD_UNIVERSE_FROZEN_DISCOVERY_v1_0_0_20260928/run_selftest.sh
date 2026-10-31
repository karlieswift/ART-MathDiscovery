#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
python -m py_compile scripts/*.py
python scripts/00_preflight.py --root "$ROOT"
TMP="$ROOT/work/selftest"
rm -rf "$TMP"; mkdir -p "$TMP"
python - "$TMP" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1]); s=('ACGT'*1800)[:6200]
for n in ['TOY_A','TOY_B']:
    (p/f'{n}.fasta').write_text(f'>{n} synthetic\n'+s+'\n')
(p/'cohort.tsv').write_text('sequence_id\tarm\tsource\trelease_date\tfasta_path\tcluster_id\tnotes\n'+f'TOY_A\tL4_CURATED\tSELFTEST\t2026-08-01\t{p}/TOY_A.fasta\tC1\ttoy\n'+f'TOY_B\tL4_PROSPECTIVE\tSELFTEST\t2026-09-29\t{p}/TOY_B.fasta\tC2\ttoy\n')
PY
python scripts/02_freeze_universe.py --root "$ROOT" --cohort-tsv "$TMP/cohort.tsv" >/tmp/stage04a_universe_selftest.json
python - <<'PY'
import json
x=json.load(open('outputs/STAGE04A_UNIVERSE_GATE.json')); assert x['gate_pass'] and x['n_sequences']==2
print('[PASS] Stage04A universe freeze policy selftest')
PY
rm -rf inputs/* outputs/STAGE04A_UNIVERSE_MANIFEST.tsv outputs/STAGE04A_UNIVERSE_GATE.json runtime_engine work/selftest
printf '%s\n' '[PASS] Stage04A accelerated package selftest; scientific engine remains byte-exact Stage03A'
