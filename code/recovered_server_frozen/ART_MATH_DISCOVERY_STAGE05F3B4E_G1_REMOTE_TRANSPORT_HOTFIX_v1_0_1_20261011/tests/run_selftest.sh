#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")/.." && pwd)
bash -n "$S/run_stage05f3b4e_global_remote_r1.sh"
grep -q 'one query per request' "$S/run_stage05f3b4e_global_remote_r1.sh"
grep -q 'technical-censored, not a biological no-hit' "$S/run_stage05f3b4e_global_remote_r1.sh"
echo '[PASS] F3B4E G1-R1 transport-hotfix syntax/policy selftest'
