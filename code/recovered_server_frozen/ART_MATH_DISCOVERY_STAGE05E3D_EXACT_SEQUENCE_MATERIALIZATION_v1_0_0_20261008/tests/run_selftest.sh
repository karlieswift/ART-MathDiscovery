#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
mkdir -p "$T/e3c/outputs" "$T/e3a/frozen" "$T/e3a/outputs"
cat > "$T/e3c/outputs/STAGE05E3C_REVIEW_QUEUE.tsv" <<'TSV'
family_id	e3c_review_tier
ARCH250_0010	E3C1_PRIORITY_BIOLOGICAL_CLOSURE
ARCH250_0013	E3C1_PRIORITY_BIOLOGICAL_CLOSURE
ARCH250_0033	E3C1_PRIORITY_BIOLOGICAL_CLOSURE
TSV
cat > "$T/e3a/frozen/STAGE05E3A_E2A_CANDIDATES.tsv" <<'TSV'
family_id	candidate_id	sequence_id	repeat_sequence
ARCH250_0010	C1	S1	ACGTACGTACGT
ARCH250_0013	C2	S2	GCGTGCGTGCGT
ARCH250_0033	C3	S3	TATCTATCTATC
TSV
cat > "$T/e3a/frozen/STAGE05E3A_E2A_WITNESS_PAIRS.tsv" <<'TSV'
family_id	candidate_id	sequence_id	repeat_core
ARCH250_0010	C1	S1	ACGTACGTACGT
ARCH250_0013	C2	S2	GCGTGCGTGCGT
ARCH250_0033	C3	S3	TATCTATCTATC
TSV
cat > "$T/e3a/outputs/STAGE05E3A_CANDIDATE_OFFICIAL_CONTEXT.tsv" <<'TSV'
candidate_id	sequence_id	official_protein_id
C1	S1	P1.1
C2	S2	P2.1
C3	S3	P3.1
TSV
cat > "$T/e3a/outputs/STAGE05E3A_NEIGHBORHOOD_OFFICIAL_CDS.tsv" <<'TSV'
sequence_id	protein_id
S1	P1.1
S2	P2.1
S3	P3.1
TSV
cat > "$T/e3a/outputs/STAGE05E3A_NEIGHBORHOOD_OFFICIAL_PROTEINS.faa" <<'FAA'
>P1.1 protein one
MKKACDEFGHIKLMNPQRSTVWYACDEFGHIKLMN
>P2.1 protein two
MNNACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQ
>P3.1 protein three
MQQACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRST
FAA
python "$HERE/scripts/stage05e3d.py" --e3c "$T/e3c" --e3a "$T/e3a" --project-root "$T" --outdir "$T/out" --workdir "$T/work" --expected-families 3 >/dev/null
python - <<PY
import json
x=json.load(open('$T/out/STAGE05E3D_GATE.json'))
assert x['gate_pass'] is True, x
assert not x['priority_missing_families'], x
PY
echo '[PASS] Stage05E3D policy/schema selftest'
