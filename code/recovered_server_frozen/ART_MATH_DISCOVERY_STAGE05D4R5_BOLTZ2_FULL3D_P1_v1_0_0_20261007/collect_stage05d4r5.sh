#!/usr/bin/env bash
set -euo pipefail
[ -s outputs/GATE_STAGE05D4R5.json ] || { echo '[FAIL] no gate'; exit 2; }
ts=$(date +%Y%m%d_%H%M%S)
out="ART_STAGE05D4R5_BOLTZ2_FULL3D_P1_RESULTS_${ts}.tar.gz"
tar -czf "$out" README.md config docs scripts run_*.sh freeze_stage05d4r5.sh frozen outputs raw logs $(find boltz_out/predictions -type f \( -name '*_model_0.cif' -o -name 'confidence_*_model_0.json' \) | sort)
sha256sum "$out" | tee "$out.sha256"
echo "$PWD/$out"
