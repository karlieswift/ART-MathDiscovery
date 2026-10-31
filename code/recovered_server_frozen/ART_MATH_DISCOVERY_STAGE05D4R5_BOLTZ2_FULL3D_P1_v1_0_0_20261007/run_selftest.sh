#!/usr/bin/env bash
set -euo pipefail
R3=${1:?R3_DIR}; R4=${2:?R4_DIR}
python scripts/freeze_r5.py --r3 "$R3" --r4 "$R4" --out /tmp/stage05d4r5_selftest_$$ >/tmp/stage05d4r5_selftest_$$.log
python - <<'PY' /tmp/stage05d4r5_selftest_$$/STAGE05D4R5_FREEZE_SUMMARY.json
import json,sys
x=json.load(open(sys.argv[1])); assert x['n_P1']==43 and x['n_P2']==40 and x['n_r4_nohit']==244
print('[PASS] Stage05D4-R5 selftest: P1=43, P2=40, R4 no-hit=244; compute-priority only')
PY
rm -rf /tmp/stage05d4r5_selftest_$$ /tmp/stage05d4r5_selftest_$$.log
