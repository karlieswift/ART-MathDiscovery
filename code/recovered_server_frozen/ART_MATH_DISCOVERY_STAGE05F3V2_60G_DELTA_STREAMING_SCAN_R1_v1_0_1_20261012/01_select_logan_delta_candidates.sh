#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
P=${ART_MATHDISCOVERY_ROOT:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
RESUME_FROM=${F3V2_RESUME_FROM:-$P/ART_MATH_DISCOVERY_STAGE05F3V2_60G_DELTA_STREAMING_SCAN_v1_0_0_20261011}
SEED=F3V2_60G_DELTA_20261011
META_URI=s3://logan-pub/stats/sra_taxid.csv.zst
mkdir -p "$S"/{inputs,logs,frozen,work/candidate_r1,tmp/sort}

# Reuse the exact old-universe exclusion list when available; otherwise recover it exactly as v1.0.0 did.
if [[ ! -s "$S/inputs/OLD_F3_ACCESSIONS.txt" ]]; then
  if [[ -s "$RESUME_FROM/inputs/OLD_F3_ACCESSIONS.txt" ]]; then
    cp "$RESUME_FROM/inputs/OLD_F3_ACCESSIONS.txt" "$S/inputs/OLD_F3_ACCESSIONS.txt"
    echo "[R1] imported old-universe exclusion list from $RESUME_FROM"
  else
    python "$S/scripts/discover_old_accessions.py" "$P" "$S/inputs/OLD_F3_ACCESSIONS.txt" | tee "$S/logs/discover_old_accessions.log"
  fi
fi
N=$(grep -cve '^$' "$S/inputs/OLD_F3_ACCESSIONS.txt" || true)
if (( N < 131 )); then
  echo "[FATAL] OLD_F3_ACCESSIONS.txt has $N accessions; require >=131." >&2
  exit 2
fi
sha256sum "$S/inputs/OLD_F3_ACCESSIONS.txt" > "$S/frozen/OLD_F3_ACCESSIONS.sha256"

# Resume if all four complete ranked candidate tables already exist and are nonempty.
complete=1
for x in A B C D; do [[ -s "$S/inputs/F3V2_DELTA_CANDIDATES_${x}.tsv" ]] || complete=0; done
if (( complete == 0 )); then
  chmod u+w "$S"/inputs/F3V2_DELTA_CANDIDATES_*.tsv 2>/dev/null || true
  rm -f "$S"/inputs/F3V2_DELTA_CANDIDATES_*.tsv
  rm -rf "$S/work/candidate_r1" "$S/tmp/sort"
  mkdir -p "$S/work/candidate_r1" "$S/tmp/sort"
  echo '[R1 STREAM] Logan metadata -> ALL source-blind eligible candidates; no fixed K ceiling'
  aws s3 cp "$META_URI" - --no-sign-request \
    | zstd -dc \
    | python "$S/scripts/stream_logan_candidates_r1.py" \
        --exclude "$S/inputs/OLD_F3_ACCESSIONS.txt" \
        --seed "$SEED" \
        --out-dir "$S/work/candidate_r1" \
        2> >(tee "$S/logs/select_logan_candidates_r1.stderr" >&2)

  for STR in A B C D; do
    raw="$S/work/candidate_r1/$STR.unsorted.tsv"
    sorted="$S/work/candidate_r1/$STR.sorted.tsv"
    out="$S/inputs/F3V2_DELTA_CANDIDATES_${STR}.tsv"
    [[ -s "$raw" ]] || { echo "[FATAL] no eligible candidates for stratum $STR" >&2; exit 2; }
    # Same v1.0.0 order: ascending 64-bit SHA prefix, accession tie-break.
    LC_ALL=C sort -T "$S/tmp/sort" -t $'\t' -k2,2n -k3,3 "$raw" > "$sorted"
    python "$S/scripts/assign_candidate_rank_r1.py" "$sorted" "$out"
    n=$(($(wc -l < "$out")-1))
    echo "[R1 CANDIDATES] $STR n=$n"
  done
fi

# A manifest of the full deterministic order replaces the old top-K pool hash.
: > "$S/frozen/F3V2_DELTA_CANDIDATE_POOLS_R1.sha256"
for STR in A B C D; do sha256sum "$S/inputs/F3V2_DELTA_CANDIDATES_${STR}.tsv" >> "$S/frozen/F3V2_DELTA_CANDIDATE_POOLS_R1.sha256"; done

python - "$S" "$SEED" "$META_URI" <<'PY'
import json,sys,hashlib,csv
S,seed,uri=sys.argv[1:]
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
counts={}
for s in 'ABCD':
 p=f'{S}/inputs/F3V2_DELTA_CANDIDATES_{s}.tsv'
 with open(p) as f: counts[s]=max(sum(1 for _ in f)-1,0)
g={'stage':'Stage05F3V2-R1-source-blind-candidate-order','gate':'PASS_F3V2_R1_UNBOUNDED_DETERMINISTIC_CANDIDATE_ORDER_FROZEN','gate_pass':all(counts[s]>0 for s in 'ABCD'),
'seed':seed,'metadata_uri':uri,'candidate_counts':counts,'old_accessions_sha256':sha(f'{S}/inputs/OLD_F3_ACCESSIONS.txt'),
'ordering':'ascending uint64(first16hex(SHA256(seed|stratum|accession))); accession tie-break; identical to v1.0.0 top-K ordering',
'acquisition_change_only':'fixed K=5000 ceiling removed; all eligible candidates are externally sorted and downstream acquisition stops only after frozen bp quota is met',
'source_selection_used_rt_or_candidate_information':False,
'selection_inputs':'SRA accession/assay/organism/taxonomy plus downstream Logan technical file size and realized nucleotide bp only'}
open(f'{S}/frozen/F3V2_R1_CANDIDATE_ORDER_FREEZE.json','w').write(json.dumps(g,indent=2)+'\n')
print(json.dumps(g,indent=2))
assert g['gate_pass']
PY

OLDPOOL="$RESUME_FROM/inputs/F3V2_DELTA_ACCESSION_CANDIDATES.tsv"
python "$S/scripts/verify_candidate_prefix_r1.py" "$OLDPOOL" "$S/inputs" "$S/frozen/F3V2_R1_PREFIX_COMPATIBILITY.json"
chmod a-w "$S"/inputs/F3V2_DELTA_CANDIDATES_*.tsv "$S/frozen/F3V2_R1_CANDIDATE_ORDER_FREEZE.json" "$S/frozen/F3V2_R1_PREFIX_COMPATIBILITY.json" 2>/dev/null || true
echo PASS_F3V2_R1_SOURCE_BLIND_UNBOUNDED_CANDIDATE_ORDER_FROZEN
