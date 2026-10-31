#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <ART_STAGE05E3A...R1.tar.gz> <project_root>" >&2; exit 2
fi
E3A_TAR=$(readlink -f "$1")
ROOT=$(readlink -f "$2")
HERE=$(cd "$(dirname "$0")" && pwd)
EXP=31c70ee1c7769baec58da11dfcc40fc3ec6ab7c17d6608391f737bc404041bca
ACT=$(sha256sum "$E3A_TAR" | awk '{print $1}')
[ "$ACT" = "$EXP" ] || { echo "[FAIL] E3A R1 SHA mismatch expected=$EXP actual=$ACT"; exit 3; }
WORK="$HERE/work/e3a_r1"
OUT="$HERE/outputs"
rm -rf "$WORK" "$OUT"
mkdir -p "$WORK" "$OUT" "$HERE/logs"
tar -xzf "$E3A_TAR" -C "$WORK"
python "$HERE/scripts/analyze_stage05e3b.py" --e3a-root "$WORK" --outdir "$OUT" 2>&1 | tee "$HERE/logs/analyze_stage05e3b.log"
cat "$OUT/GATE_STAGE05E3B.json"
