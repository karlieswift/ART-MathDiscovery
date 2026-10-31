#!/usr/bin/env bash
set -euo pipefail
python -m py_compile scripts/build_architecture_families.py
grep -q 'construction_independent_of_R5_structure_results' config/POLICY.json
grep -q '250' config/POLICY.json
grep -q '1000' config/POLICY.json
echo '[PASS] Stage05E1 policy/script selftest'
