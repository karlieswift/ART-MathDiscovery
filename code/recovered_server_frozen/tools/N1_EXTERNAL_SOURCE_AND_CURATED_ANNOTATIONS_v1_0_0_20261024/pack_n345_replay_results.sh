#!/usr/bin/env bash
set -euo pipefail
P=${1:?Usage: bash pack_n345_replay_results.sh /absolute/project/root}
E="$P/exports"
test -f "$E/NCS_N1345_20261024/N5_VERIFY.json" || { echo 'Missing N5 verify results';exit 2; }
test -f "$E/NCS_N3_FULL_MATRIX_20261024/N3_SENSITIVITY_6_NEW_MATRICES.tsv" || { echo 'Missing full N3 replay';exit 2; }
NOW="$(date +%Y%m%d_%H%M%S)_$$"
A="$E/NCS_N345_VERIFIED_RESULTS_${NOW}.tar.gz"
tar -czf "$A" -C "$E" NCS_N1345_20261024 NCS_N3_FULL_MATRIX_20261024
sha256sum "$A" > "$A.sha256"
sha256sum -c "$A.sha256"
ls -lh "$A" "$A.sha256"
printf 'SEND_ARCHIVE=%s\nSEND_SHA256=%s\n' "$A" "$A.sha256"
