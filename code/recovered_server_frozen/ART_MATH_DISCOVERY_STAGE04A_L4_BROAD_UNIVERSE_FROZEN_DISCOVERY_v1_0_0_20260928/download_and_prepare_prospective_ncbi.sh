#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
TAXON="${STAGE04_TAXON:-Caudoviricetes}"
D="$ROOT/work/ncbi_prospective_20260929"
mkdir -p "$D"
command -v datasets >/dev/null || { echo 'ERROR: NCBI datasets CLI not found'; exit 2; }
echo "[1/4] metadata: $TAXON complete genomes first released >= 2026-09-29"
datasets summary virus genome taxon "$TAXON" --complete-only --released-after 2026-09-29 --as-json-lines > "$D/data_report.jsonl"
python - "$D/data_report.jsonl" "$D/accessions.txt" <<'PY'
import json,re,sys
inp,out=sys.argv[1:3]; ex={'MW248466','MW349129','MW349128','MW218148'}; acc=[]
def find(o):
 if isinstance(o,dict):
  for k,v in o.items():
   if re.sub('[^a-z]','',k.lower()) in {'accession','nucleotideaccession'} and isinstance(v,str): return v
  for v in o.values():
   x=find(v)
   if x:return x
 elif isinstance(o,list):
  for v in o:
   x=find(v)
   if x:return x
for line in open(inp):
 try:a=find(json.loads(line))
 except Exception:continue
 if a and a.split('.')[0] not in ex: acc.append(a)
acc=sorted(set(acc)); open(out,'w').write('\n'.join(acc)+('\n' if acc else '')); print('selected',len(acc))
PY
[[ -s "$D/accessions.txt" ]] || { echo '[INFO] No prospective genomes available yet; retry later without changing policy.'; exit 0; }
echo '[2/4] download sequence-only NCBI package'
datasets download virus genome accession --inputfile "$D/accessions.txt" --include genome --filename "$D/ncbi_prospective.zip"
rm -rf "$D/unzip"; unzip -q "$D/ncbi_prospective.zip" -d "$D/unzip"
echo '[3/4] prepare cohort'
python "$ROOT/scripts/01_prepare_ncbi_cohort.py" --genomic-fna "$D/unzip/ncbi_dataset/data/genomic.fna" --data-report-jsonl "$D/data_report.jsonl" --arm L4_PROSPECTIVE --source "NCBI_NEW_${TAXON}" --outdir "$D/prepared"
echo '[4/4] ready'
echo "RUN: $ROOT/run_stage04a_primary.sh $D/prepared/STAGE04A_L4_PROSPECTIVE_COHORT.tsv 4 8 0"
