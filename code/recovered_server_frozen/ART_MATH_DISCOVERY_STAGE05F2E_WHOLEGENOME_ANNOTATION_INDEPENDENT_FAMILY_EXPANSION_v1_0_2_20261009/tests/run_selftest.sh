#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
python -m py_compile scripts/stage05f2e.py
python scripts/stage05f2e.py --selftest
python - <<'PY'
import json
p=json.load(open("config/FROZEN_POLICY.json"))
assert p["expected_query_n"]==30
assert p["expected_genome_n"]==919
assert p["annotation_used_to_generate_peptides"] is False
assert p["dna_score_recomputed"] is False
assert p["translation_table_codon_n"]==64
assert p["translation_table_sha256"]=="1cd328cab4f124096685dc86bb5dd6badefe8d417c4af97ce68ba789ff560745"
print("[PASS] Stage05F2E frozen policy selftest")
PY
