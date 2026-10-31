#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HERE"
rm -rf tests/_work tests/_out tests/_selftest_stage05e3c_stdout.json
python scripts/stage05e3c.py \
  --e3b tests/fixtures/e3b \
  --e3a tests/fixtures/e3a \
  --project-root "$HERE" \
  --outdir tests/_out \
  --workdir tests/_work \
  --expected-families 2 > tests/_selftest_stage05e3c_stdout.json
python - <<'PY'
import csv,json
s=json.load(open('tests/_out/STAGE05E3C_SUMMARY.json'))
assert s['gate_pass'] is True and s['n_families']==2
rows=list(csv.DictReader(open('tests/_out/STAGE05E3C_REVIEW_QUEUE.tsv'),delimiter='\t'))
by={r['family_id']:r for r in rows}
assert by['ARCH250_TEST1']['e3c_review_tier']=='E3C1_PRIORITY_BIOLOGICAL_CLOSURE'
assert by['ARCH250_TEST2']['e3c_review_tier']=='E3C5_KNOWN_SYSTEM_OR_SPECIAL_CONFLICT_REVIEW'
a=json.load(open('tests/_out/STAGE05E3C_ART95_AUDIT_STATUS.json'))
assert a['author_reported_family_member_n']==95 and a['exact_manifest_ready'] is False
print('[PASS] Stage05E3C policy/script selftest')
PY
rm -rf tests/_work tests/_out tests/_selftest_stage05e3c_stdout.json
