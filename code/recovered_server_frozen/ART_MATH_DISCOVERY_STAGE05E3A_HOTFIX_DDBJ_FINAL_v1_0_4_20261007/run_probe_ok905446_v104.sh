#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 1 ]; then echo "usage: $0 <E3A_ROOT>" >&2; exit 2; fi
E3A="$(cd "$1" && pwd)"
FIX="$(cd "$(dirname "$0")" && pwd)"
"$FIX/run_selftest.sh"
rm -rf /tmp/e3a_ddbj_probe_v104
mkdir -p /tmp/e3a_ddbj_probe_v104/out
python "$FIX/scripts/recover_ok905446_ddbj_v104.py" \
  --outdir /tmp/e3a_ddbj_probe_v104/out \
  --manifest /tmp/e3a_ddbj_probe_v104/manifest.tsv \
  --diagnostics /tmp/e3a_ddbj_probe_v104/diagnostics.tsv
column -t -s $'\t' /tmp/e3a_ddbj_probe_v104/diagnostics.tsv
cat /tmp/e3a_ddbj_probe_v104/manifest.tsv
