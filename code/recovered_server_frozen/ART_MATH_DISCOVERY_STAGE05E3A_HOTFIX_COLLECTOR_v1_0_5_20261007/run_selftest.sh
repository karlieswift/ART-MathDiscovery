#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
test -x "$ROOT/collect_stage05e3a_v105.sh"
grep -q 'STAGE05E3A_COMBINED_PROVIDER_MANIFEST_v104.tsv' "$ROOT/collect_stage05e3a_v105.sh"
if grep -q 'tar -czf.*raw/STAGE05E3A_GENBANK_RETRIEVAL_MANIFEST.tsv' "$ROOT/collect_stage05e3a_v105.sh"; then
  echo '[FAIL] obsolete single-provider manifest still hard-coded' >&2
  exit 1
fi
echo '[PASS] Stage05E3A v1.0.5 collector-only hotfix selftest'
