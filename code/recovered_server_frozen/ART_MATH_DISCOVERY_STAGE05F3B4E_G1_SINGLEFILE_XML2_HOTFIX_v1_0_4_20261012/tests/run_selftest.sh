#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")/.." && pwd)
bash -n "$S/run_stage05f3b4e_global_urlapi_r4_singlefile.sh"
bash -n "$S/check_ready_rid_xml2s.sh"
python -m py_compile "$S/scripts/ncbi_blast_urlapi_singlefile.py" "$S/scripts/summarize_global.py"
grep -q "FORMAT_TYPE':'XML2_S'" "$S/scripts/ncbi_blast_urlapi_singlefile.py"
! grep -q "FORMAT_TYPE':'XML2'" "$S/scripts/ncbi_blast_urlapi_singlefile.py"
echo '[PASS] R4 single-file XML2 retrieval syntax/policy selftest'
