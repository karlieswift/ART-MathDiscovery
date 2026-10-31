#!/usr/bin/env bash
set -euo pipefail

S="$(cd "$(dirname "$0")" && pwd)"

THREADS="${1:-32}"

RTDB="/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1/ART_MATH_DISCOVERY_STAGE05D4R5C_CLEAN_RT_CALIBRATION_v1_0_0_20261007/references/clean_rt_afsp"

STATUS="$S/outputs/F3V2_STAGE04B3_BOLTZ_STATUS.tsv"

mkdir -p \
  "$S/raw/foldseek" \
  "$S/tmp/foldseek" \
  "$S/logs/stage5_resume"

run_foldseek() {
    conda run -n art_novelty foldseek "$@"
}

echo "=== Stage04B3 Stage5 technical resume ==="
echo "science_change=false"
echo "threshold_change=false"
echo "Boltz_rerun=false"
echo "historical_R5C_rule_change=false"
echo

python - "$STATUS" "$S" <<'PY' > "$S/outputs/F3V2_STAGE04B3_CIF_RESOLVED.tsv"
from pathlib import Path
import csv
import sys

status = Path(sys.argv[1])
root = Path(sys.argv[2])

print("\t".join([
    "rank",
    "family",
    "candidate_id",
    "aa_length",
    "recorded_cif",
    "resolved_cif",
    "recorded_cif_exists",
    "resolution_status",
]))

with status.open() as f:
    r = csv.reader(f, delimiter="\t")
    header = next(r)

    for row in r:

        if len(row) < 10:
            continue

        rank, fam, cid, acc, aa_len, gpu, bst, rc, elapsed, recorded = row[:10]

        if not bst.startswith("PASS"):
            continue

        rec = Path(recorded)

        actual = sorted(
            p.resolve()
            for p in (root / "boltz_out" / cid).rglob("*_model_0.cif")
            if p.is_file()
        )

        if len(actual) != 1:
            raise SystemExit(
                f"FAIL: {cid}: expected exactly one model_0 CIF, found {len(actual)}"
            )

        print("\t".join([
            rank,
            fam,
            cid,
            aa_len,
            recorded,
            str(actual[0]),
            str(rec.is_file()),
            "UNIQUE_ACTUAL_CIF_RESOLVED",
        ]))
PY

echo "[PASS] resolved all Boltz-PASS structures to unique actual CIFs"

tail -n +2 "$S/outputs/F3V2_STAGE04B3_CIF_RESOLVED.tsv" |
while IFS=$'\t' read -r rank fam cid aa_len recorded_cif cif recorded_exists resolution; do

    out="$S/raw/foldseek/${cid}.m8"
    tmp="$S/tmp/foldseek/${cid}"
    log="$S/logs/stage5_resume/foldseek_${cid}.log"
    done_marker="$S/logs/stage5_resume/foldseek_${cid}.PASS"

    if [[ -f "$done_marker" ]]; then
        echo "[SKIP COMPLETE] rank=$rank cid=$cid"
        continue
    fi

    echo
    echo "===================================================="
    echo "[RUN] rank=$rank aa=$aa_len"
    echo "CID=$cid"
    echo "CIF=$cif"
    echo "===================================================="

    test -s "$cif" || {
        echo "[FAIL] CIF missing: $cif"
        exit 20
    }

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
        echo "[FAIL] Foldseek rc=$rc cid=$cid"
        tail -80 "$log"
        exit "$rc"
    fi

    touch "$done_marker"

    echo "[PASS] foldseek $cid"
done

echo
echo "=== STAGE5 RESUME COMPLETE ==="

expected=$(
  tail -n +2 "$S/outputs/F3V2_STAGE04B3_CIF_RESOLVED.tsv" |
  wc -l
)

done_n=$(
  find "$S/logs/stage5_resume" \
    -maxdepth 1 \
    -type f \
    -name '*.PASS' |
  wc -l
)

echo "expected=$expected"
echo "completed=$done_n"

if [[ "$done_n" -ne "$expected" ]]; then
    echo "[FAIL] incomplete Stage5 resume"
    exit 30
fi

echo "PASS_STAGE04B3_STAGE5_FOLDSEEK_TECHNICAL_RESUME_COMPLETE"
