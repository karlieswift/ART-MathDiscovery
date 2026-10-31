#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); cd "$HERE"
stamp=$(date +%Y%m%d_%H%M%S)
out="ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B7_RESULTS_${stamp}.tar.gz"
tar -czf "$out" README.md config docs outputs logs/*.log 2>/dev/null || tar -czf "$out" README.md config docs outputs
sha256sum "$out" | tee "$out.sha256"
echo "RESULT=$out"
