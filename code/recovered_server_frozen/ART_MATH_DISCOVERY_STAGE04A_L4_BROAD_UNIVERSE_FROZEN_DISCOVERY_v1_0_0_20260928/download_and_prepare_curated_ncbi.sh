#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
TAXON="${STAGE04_TAXON:-Caudoviricetes}"
SAMPLE_N="${STAGE04_SAMPLE_N:-250}"
D="$ROOT/work/ncbi_curated"
mkdir -p "$D"
command -v datasets >/dev/null || { echo 'ERROR: NCBI datasets CLI not found'; exit 2; }
echo "[1/4] metadata: $TAXON complete RefSeq"
datasets summary virus genome taxon "$TAXON" --refseq --complete-only --as-json-lines > "$D/data_report.jsonl"
python - "$D/data_report.jsonl" "$D/accessions.txt" "$SAMPLE_N" <<'PY'
import json,hashlib,re,sys
inp,out,n=sys.argv[1],sys.argv[2],int(sys.argv[3]); ex={'MW248466','MW349129','MW349128','MW218148'}; acc=[]
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
acc=sorted(set(acc),key=lambda a:hashlib.sha256(('STAGE04A_CURATED_SAMPLE_V1|'+a).encode()).hexdigest())
if n>0: acc=acc[:n]
open(out,'w').write('\n'.join(acc)+'\n'); print('selected',len(acc))
PY
echo '[2/4] download sequence-only NCBI package for deterministic cohort'
datasets download virus genome accession --inputfile "$D/accessions.txt" --include genome --filename "$D/ncbi_curated.zip"
rm -rf "$D/unzip"; unzip -q "$D/ncbi_curated.zip" -d "$D/unzip"
echo '[3/4] prepare per-genome FASTAs + cohort TSV'
python "$ROOT/scripts/01_prepare_ncbi_cohort.py" --genomic-fna "$D/unzip/ncbi_dataset/data/genomic.fna" --data-report-jsonl "$D/data_report.jsonl" --arm L4_CURATED --source "NCBI_REFSEQ_${TAXON}" --sample-n 0 --outdir "$D/prepared"
echo '[4/4] ready'
echo "COHORT=$D/prepared/STAGE04A_L4_CURATED_COHORT.tsv"
echo "RUN: $ROOT/run_stage04a_primary.sh $D/prepared/STAGE04A_L4_CURATED_COHORT.tsv 4 8 0"
