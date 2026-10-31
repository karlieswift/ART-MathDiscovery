#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")/.." && pwd)
bash -n "$S/run_stage05f3b4e_g1_dual_rid_resume_and_finalize.sh"
python -m py_compile "$S/scripts/ncbi_blast_urlapi_singlefile.py" "$S/scripts/summarize_global.py" "$S/scripts/common.py"
grep -q -- '--resume-rid "$PART_RID"' "$S/run_stage05f3b4e_g1_dual_rid_resume_and_finalize.sh"
grep -q -- '--resume-rid "$ARCH_RID"' "$S/run_stage05f3b4e_g1_dual_rid_resume_and_finalize.sh"
! grep -q -- '--resubmit-if-expired' "$S/run_stage05f3b4e_g1_dual_rid_resume_and_finalize.sh"
echo '[PASS] R5 dual-existing-RID resume/finalize syntax/policy selftest'
