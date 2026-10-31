#!/usr/bin/env bash
set -euo pipefail

S="$(cd "$(dirname "$0")" && pwd)"

THREADS="${1:-32}"

STATUS="$S/outputs/F3V2_STAGE04B3_BOLTZ_STATUS.tsv"

RTDB="/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE05D4R5C_CLEAN_RT_CALIBRATION_v1_0_0_20261007/references/clean_rt_afsp"

R1="$S/technical_repair_probe_gpu0/boltz_out"
R2="$S/technical_repair_r1/gpu0_2624_retry/boltz_out"

OUT="$S/raw/foldseek_repair9"
TMP="$S/tmp/foldseek_repair9"
LOG="$S/logs/stage5_repair9"

mkdir -p "$OUT" "$TMP" "$LOG"

echo "=== RESOLVE THE 9 TECHNICALLY RECOVERED CIFs ==="

python - "$STATUS" "$R1" "$R2" <<'PY' \
  > "$S/outputs/F3V2_STAGE04B3_REPAIRED9_CIF_RESOLVED.tsv"

from pathlib import Path
import csv
import sys

status = Path(sys.argv[1])
roots = [Path(sys.argv[2]), Path(sys.argv[3])]

expected_ranks = {2,6,9,10,14,18,22,26,30}

rows = []

with status.open(newline="") as f:
    r = csv.DictReader(f, delimiter="\t")

    for row in r:
        # remove possible CR from all values
        row = {
            k: (v.rstrip("\r") if isinstance(v, str) else v)
            for k, v in row.items()
        }

        if row["boltz_status"].startswith("PASS"):
            continue

        rank = int(row["shortlist_rank"])
        cid = row["candidate_id"]

        hits = []

        for root in roots:
            if not root.exists():
                continue

            candidate_root = root / cid

            if candidate_root.exists():
                hits.extend(
                    p.resolve()
                    for p in candidate_root.rglob("*_model_0.cif")
                    if p.is_file()
                )

        # unique paths only
        hits = sorted(set(hits))

        if len(hits) != 1:
            raise SystemExit(
                f"FAIL rank={rank} cid={cid}: "
                f"expected exactly one repaired model_0.cif, found {len(hits)}"
            )

        rows.append((
            rank,
            row["family_rep_id"],
            cid,
            row["accession"],
            row["aa_length"],
            str(hits[0]),
        ))

got = {x[0] for x in rows}

if got != expected_ranks:
    raise SystemExit(
        f"FAIL repaired rank set mismatch: "
        f"got={sorted(got)} expected={sorted(expected_ranks)}"
    )

if len(rows) != 9:
    raise SystemExit(f"FAIL expected 9 repaired candidates, got {len(rows)}")

print("\t".join([
    "shortlist_rank",
    "family_rep_id",
    "candidate_id",
    "accession",
    "aa_length",
    "repaired_cif",
]))

for row in sorted(rows):
    print("\t".join(map(str, row)))
PY

column -t -s $'\t' \
  "$S/outputs/F3V2_STAGE04B3_REPAIRED9_CIF_RESOLVED.tsv"

echo
echo "=== RUN EXACT HISTORICAL R5C FOLDSEEK ON REPAIRED 9 ==="

run_foldseek() {
    conda run -n art_novelty foldseek "$@"
}

tail -n +2 "$S/outputs/F3V2_STAGE04B3_REPAIRED9_CIF_RESOLVED.tsv" |
while IFS=$'\t' read -r rank fam cid acc aa_len cif; do

    out="$OUT/${cid}.m8"
    tmp="$TMP/${cid}"
    log="$LOG/foldseek_${cid}.log"
    marker="$LOG/foldseek_${cid}.PASS"

    if [[ -f "$marker" ]]; then
        echo "[SKIP COMPLETE] rank=$rank cid=$cid"
        continue
    fi

    echo
    echo "============================================================"
    echo "[RUN] repaired rank=$rank aa=$aa_len"
    echo "CID=$cid"
    echo "CIF=$cif"
    echo "============================================================"

    test -s "$cif"

    rm -rf "$tmp"
    rm -f "$out"

    set +e

    run_foldseek easy-search \
      "$cif" \
      "$RTDB" \
      "$out" \
      "$tmp" \
      -e 1.0 \
      --max-seqs 100 \
      --threads "$THREADS" \
      --exact-tmscore 1 \
      --format-output \
      'query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore' \
      > "$log" 2>&1

    rc=$?

    set -e

    if [[ "$rc" -ne 0 ]]; then
        echo "[FAIL] Foldseek rc=$rc rank=$rank cid=$cid"
        tail -100 "$log"
        exit "$rc"
    fi

    touch "$marker"

    echo "[PASS] repaired9 foldseek $cid"
done

echo
echo "=== REPAIRED9 FOLDSEEK CLOSEOUT ==="

NMARK=$(
    find "$LOG" \
      -maxdepth 1 \
      -type f \
      -name '*.PASS' \
      | wc -l
)

NOUT=$(
    find "$OUT" \
      -maxdepth 1 \
      -type f \
      -name '*.m8' \
      | wc -l
)

echo "PASS markers = $NMARK / 9"
echo "m8 outputs   = $NOUT / 9"

[[ "$NMARK" -eq 9 ]]
[[ "$NOUT"  -eq 9 ]]

find "$OUT" \
  -maxdepth 1 \
  -type f \
  -name '*.m8' \
  -print0 \
| sort -z \
| xargs -0 sha256sum \
> "$S/outputs/STAGE04B3_REPAIRED9_FOLDSEEK_SHA256.tsv"

echo
echo "PASS_STAGE04B3_REPAIRED9_FOLDSEEK_TECHNICAL_CLOSURE"
