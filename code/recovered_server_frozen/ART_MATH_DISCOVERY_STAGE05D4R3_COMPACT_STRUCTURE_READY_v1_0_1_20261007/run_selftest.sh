#!/usr/bin/env bash
set -euo pipefail
D1="${1:?D1 path required}"
D4="${2:?D4 path required}"
R1="${3:?R1 path required}"
R2="${4:?R2 path required}"
python tests/selftest.py
for f in \
  "$D1/config/NOVELTY_FIRST_999_CANDIDATES.tsv" \
  "$D1/config/NOVELTY_FIRST_24020_ORF_REGISTRY.tsv" \
  "$D1/config/NOVELTY_FIRST_CANDIDATE_ORF_LINKS.tsv" \
  "$D4/outputs/STAGE05D4_N1N2_PROTEIN_UNIVERSE.tsv" \
  "$D4/outputs/STAGE05D4_N1N2_NONFORMAL_BROAD_RT_PROTEINS.faa" \
  "$R1/outputs/STAGE05D4R1_NR90_CLUSTERS.tsv" \
  "$R2/outputs/STAGE05D4R2_REMOTE_HOMOLOGY_TRIAGE.tsv"; do
  test -s "$f" || { echo "[FAIL] missing $f" >&2; exit 2; }
done
python - "$D1" "$D4" "$R1" "$R2" <<'PY'
import sys,pandas as pd,json
D1,D4,R1,R2=sys.argv[1:]
c=pd.read_csv(f"{D1}/config/NOVELTY_FIRST_999_CANDIDATES.tsv",sep="\t")
n12=c[c.novelty_architecture_stratum.astype(str).str.startswith(("N1_","N2_"))]
assert len(n12)==602, len(n12)
r=pd.read_csv(f"{R2}/outputs/STAGE05D4R2_REMOTE_HOMOLOGY_TRIAGE.tsv",sep="\t")
assert len(r)==9245, len(r)
assert (r.triage_label=="NO_SWISSPROT_HIT_NOT_NEGATIVE").sum()==8062
print("[PASS] source invariants: N1/N2=602, NR90=9245, SwissProt-no-hit=8062")
PY
