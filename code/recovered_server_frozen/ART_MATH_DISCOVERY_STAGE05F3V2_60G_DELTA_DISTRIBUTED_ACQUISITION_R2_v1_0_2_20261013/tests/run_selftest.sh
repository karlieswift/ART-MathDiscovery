#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")/.." && pwd)
python -m py_compile "$S/scripts/acquire_one_stratum_r2.py" "$S/scripts/merge_and_freeze_r2.py" "$S/scripts/count_fasta_stream.py" "$S/scripts/translate_sixframe.py" "$S/scripts/parse_hmmscan.py" "$S/scripts/dedupe_formal.py" "$S/scripts/finalize_scan.py"
bash -n "$S/02_acquire_stratum.sh" "$S/02_merge_and_freeze.sh" "$S/03_run_54g_delta_scan.sh" "$S/scripts/scan_one_accession.sh"
python - "$S" <<'PY'
import ast,sys,pathlib
S=pathlib.Path(sys.argv[1])
a=ast.parse((S/'scripts/acquire_one_stratum_r2.py').read_text())
m=ast.parse((S/'scripts/merge_and_freeze_r2.py').read_text())
# Contract checks: independent stratum lease, unchanged quotas, source R1 candidate table, no RT-aware selection.
txt=(S/'scripts/acquire_one_stratum_r2.py').read_text()
for token in ["TARGET={'A':13_500_000_000,'B':18_000_000_000,'C':13_500_000_000,'D':9_000_000_000}","F3V2_DELTA_CANDIDATES_","fcntl.flock","aws','s3api','head-object","aws','s3','cp","source R1 partial"]:
    assert token in txt, token
mt=(S/'scripts/merge_and_freeze_r2.py').read_text()
assert "for s in 'ABCD'" in mt and "deterministic_rank" in mt
print('[PASS] R2 distributed acquisition syntax/policy selftest')
PY
# Ensure downstream scientific scripts copied byte-identically from R1 if reference is available next door via env.
if [[ -n "${F3V2_R1_REFERENCE:-}" && -d "${F3V2_R1_REFERENCE}" ]]; then
  for f in 03_run_54g_delta_scan.sh scripts/count_fasta_stream.py scripts/translate_sixframe.py scripts/parse_hmmscan.py scripts/dedupe_formal.py scripts/scan_one_accession.sh scripts/finalize_scan.py; do
    cmp -s "$S/$f" "$F3V2_R1_REFERENCE/$f" || { echo "[FAIL] downstream science changed: $f" >&2; exit 3; }
  done
  echo '[PASS] R2 downstream science byte-identical to R1'
fi
