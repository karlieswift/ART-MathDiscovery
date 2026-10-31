#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 1 ]; then echo "usage: $0 <E3A_ROOT>" >&2; exit 2; fi
E3A="$(cd "$1" && pwd)"; FIX="$(cd "$(dirname "$0")" && pwd)"
mkdir -p /tmp/e3a_v103_probe/out
printf 'NC_005083.2\nOK905446.1\n' > /tmp/e3a_v103_probe/ids.txt
python "$FIX/scripts/recover_missing_v103.py" \
  --ids /tmp/e3a_v103_probe/ids.txt \
  --outdir /tmp/e3a_v103_probe/out \
  --manifest /tmp/e3a_v103_probe/manifest.tsv \
  --diagnostics /tmp/e3a_v103_probe/diagnostics.tsv
column -t -s $'\t' /tmp/e3a_v103_probe/manifest.tsv
column -t -s $'\t' /tmp/e3a_v103_probe/diagnostics.tsv | head -30
