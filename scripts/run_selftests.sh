#!/usr/bin/env bash
set -euo pipefail
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$R"
python3 scripts/verify_release.py
python3 code/arm_m/tests/test_engine.py
python3 code/clean_two_arm_replay/tests/test_two_arm.py
python3 code/formal_baseline/tests/run_selftest.py
python3 code/ug27_exact_truth/tests/test_logic.py
bash code/ug27_preunblind/run_selftest.sh
python3 -m unittest discover -s code/j1_j2_supplement/j1/tests -p 'test_*.py' -v
OUT_TEST="$(mktemp -d)"
trap 'rm -rf "$OUT_TEST"' EXIT
bash code/j1_j2_supplement/j2/run_j2.sh test --output "$OUT_TEST/J2_GUARDS.json"
printf '%s\n' '[SMOKE_TESTS] COMPLETE: source guards only, not full scientific benchmarks'
