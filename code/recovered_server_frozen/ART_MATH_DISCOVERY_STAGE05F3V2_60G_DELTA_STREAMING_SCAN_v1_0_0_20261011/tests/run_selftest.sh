#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")/.." && pwd); T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cat > "$T/x.fa" <<'EOF'
>ctg1
ATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAAATGAAATAA
EOF
python "$S/scripts/translate_sixframe.py" --accession TEST --min-aa 45 --stats "$T/stats.json" < "$T/x.fa" > "$T/o.faa"
grep -q '^>TEST|ctg1|+|f0|' "$T/o.faa"
python -m py_compile "$S/scripts/"*.py
for x in "$S"/*.sh "$S/scripts/"*.sh; do bash -n "$x"; done
echo '[PASS] F3V2 translation + syntax selftest'
