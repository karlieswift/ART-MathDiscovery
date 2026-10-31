#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
python scripts/verify.py --root "$ROOT"
OUT="$(mktemp -d "$ROOT/work_reproduce_XXXXXX")"
python scripts/prepare.py --root "$ROOT" --out "$OUT"
python - "$ROOT/precomputed" "$OUT" <<'PY'
from pathlib import Path
import sys,hashlib,json
p,q=map(Path,sys.argv[1:]);n=0
for a in p.iterdir():
 if not a.is_file():continue
 b=q/a.name
 if a.suffix in ('.tsv','.fasta'):
  if a.read_bytes()!=b.read_bytes():raise SystemExit('Reproduction difference: '+a.name)
  n+=1
 elif a.suffix=='.json':
  x,y=json.loads(a.read_text()),json.loads(b.read_text());x.pop('Bio_version',None);y.pop('Bio_version',None)
  if x!=y:raise SystemExit('Reproduction difference: '+a.name)
  n+=1
print('[PASS] Reproduced files:',n)
PY
