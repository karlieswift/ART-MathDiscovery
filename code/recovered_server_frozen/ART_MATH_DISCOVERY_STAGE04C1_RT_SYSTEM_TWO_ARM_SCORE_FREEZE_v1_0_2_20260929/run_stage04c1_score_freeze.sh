#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -lt 2 ] || [ "$#" -gt 4 ]; then
 echo "Usage: $0 C0_ROOT C0R3_SEQUENCE_FREEZE_ROOT [THREADS=32] [KEEP_RAW=0]" >&2; exit 2
fi
ROOT=$(cd "$(dirname "$0")" && pwd); C0=$(cd "$1" && pwd); SF=$(cd "$2" && pwd); THREADS=${3:-32}; KEEP=${4:-0}
"$ROOT/run_selftest.sh"
python "$ROOT/scripts/00_preflight.py" --root "$ROOT" --c0-root "$C0" --seqfreeze-root "$SF" | tee "$ROOT/logs/00_preflight.log"
python - "$ROOT" "$SF" "$THREADS" "$KEEP" <<'PY'
import csv,json,subprocess,sys
from pathlib import Path
root=Path(sys.argv[1]); sf=Path(sys.argv[2]); threads=sys.argv[3]; keep=sys.argv[4]
pol=json.loads((root/'config/C1_POLICY.json').read_text()); man={r['control_id']:r for r in csv.DictReader((sf/'outputs/STAGE04C0_R3_FIXED_SEQUENCE_MANIFEST.tsv').open(),delimiter='\t')}
for cid in pol['eligible_control_ids']:
 m=man[cid]; acc=m['accession']; fa=sf/'inputs/sequences'/f'{cid}__{acc}.fasta'; log=root/'logs'/f'01_scan_{cid}.log'
 cmd=[sys.executable,str(root/'scripts/01_scan_two_arm.py'),'--root',str(root),'--control-id',cid,'--accession',acc,'--fasta',str(fa),'--threads',threads,'--keep-raw',keep]
 print('[RUN]',' '.join(cmd),flush=True)
 with log.open('w') as f: subprocess.run(cmd,check=True,stdout=f,stderr=subprocess.STDOUT)
 print(f'[PASS] {cid}',flush=True)
PY
python "$ROOT/scripts/02_finalize.py" "$ROOT" | tee "$ROOT/logs/02_finalize.log"
echo '[SEALED] C1 score boards frozen. Do not inspect known locus coordinates until separate post-freeze unblinding.'
