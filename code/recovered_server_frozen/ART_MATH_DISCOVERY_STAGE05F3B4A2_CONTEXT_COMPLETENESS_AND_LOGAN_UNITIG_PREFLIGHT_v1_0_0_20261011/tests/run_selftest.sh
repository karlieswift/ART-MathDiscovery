#!/usr/bin/env bash
set -euo pipefail
python -m py_compile scripts/audit_context_completeness.py scripts/scan_boundary_anchors.py
bash -n scripts/check_logan_unitigs.sh scripts/download_logan_unitigs.sh run_stage05f3b4a2.sh
echo '[PASS] F3B4A2 syntax selftest'
