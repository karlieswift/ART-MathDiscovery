#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
for f in "$S/outputs/STAGE05F3B4E_G1_GATE.json" "$S/outputs/STAGE05F3B4E_G2_GATE.json" "$S/outputs/F3B4E_GLOBAL_PUBLIC_AUDIT_5.tsv" "$S/outputs/F3B4E_PUBLIC_ART8_FAMILY_SUMMARY_5.tsv"; do [[ -s "$f" ]] || { echo "[FATAL] missing $f" >&2; exit 2; }; done
python "$S/scripts/finalize.py" "$S/outputs/F3B4E_TARGETS_5.tsv" "$S/outputs/F3B4E_GLOBAL_PUBLIC_AUDIT_5.tsv" "$S/outputs/F3B4E_PUBLIC_ART8_FAMILY_SUMMARY_5.tsv" "$S/outputs" "$S/frozen/F3B4E_IMPLEMENTATION_LOCK.json" | tee "$S/logs/finalize.log"
(cd "$S" && sha256sum frozen/F3B4E_IMPLEMENTATION_LOCK.json inputs/* outputs/* > outputs/STAGE05F3B4E_SHA256.tsv)
echo PASS_STAGE05F3B4E_FINALIZE_COMPLETE
