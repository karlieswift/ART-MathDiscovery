#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
C="$S/inputs/F3V2_DELTA_ACCESSION_CANDIDATES.tsv"; OLD="$S/inputs/OLD_F3_ACCESSIONS.txt"
[[ -s "$C" && -s "$OLD" ]] || { echo '[FATAL] run 01_select_logan_delta_candidates.sh first' >&2; exit 2; }
mkdir -p "$S"/{data/contigs,work/download,logs,frozen}
OUT="$S/frozen/F3V2_DELTA_SELECTED_ACCESSIONS.tsv"; TMP="$S/work/download/selected.tmp.tsv"
# Frozen 54-Gb delta quotas scale the original 6-Gb 1.5/2.0/1.5/1.0 design by 9x.
declare -A TARGET=( [A]=13500000000 [B]=18000000000 [C]=13500000000 [D]=9000000000 )
MAX_ACC_BP=${MAX_ACCESSION_BP:-1500000000}
MIN_ZST_BYTES=${MIN_ZST_BYTES:-10000000}
MAX_ZST_BYTES=${MAX_ZST_BYTES:-2000000000}
printf 'stratum\taccession\tnucleotide_bp\tcontig_n\tcompressed_bytes\tzst_sha256\ts3_uri\n' > "$TMP"
for STR in A B C D; do
  got=0
  mapfile -t ACCS < <(awk -F'\t' -v s="$STR" 'NR>1 && $1==s {print $3}' "$C")
  echo "[STRATUM $STR] target_bp=${TARGET[$STR]} candidates=${#ACCS[@]}"
  for acc in "${ACCS[@]}"; do
    (( got >= TARGET[$STR] )) && break
    grep -qxF "$acc" "$OLD" && { echo "[FATAL] old accession leaked into delta: $acc" >&2; exit 2; }
    key="c/$acc/$acc.contigs.fa.zst"; uri="s3://logan-pub/$key"; z="$S/data/contigs/$acc.contigs.fa.zst"
    if [[ ! -s "$z" ]]; then
      bytes=$(aws s3api head-object --bucket logan-pub --key "$key" --no-sign-request --query ContentLength --output text 2>/dev/null || echo 0)
      [[ "$bytes" =~ ^[0-9]+$ ]] || bytes=0
      if (( bytes < MIN_ZST_BYTES || bytes > MAX_ZST_BYTES )); then continue; fi
      echo "[DOWNLOAD] $acc compressed=$bytes"
      aws s3 cp "$uri" "$z.part" --no-sign-request --only-show-errors || { rm -f "$z.part"; continue; }
      mv "$z.part" "$z"
    fi
    bytes=$(stat -c %s "$z")
    js=$(zstd -dc "$z" | python "$S/scripts/count_fasta_stream.py") || { echo "[WARN] decompression/count failed $acc"; rm -f "$z"; continue; }
    bp=$(python -c 'import json,sys;print(json.loads(sys.argv[1])["bp"])' "$js")
    nr=$(python -c 'import json,sys;print(json.loads(sys.argv[1])["records"])' "$js")
    if (( bp <= 0 || bp > MAX_ACC_BP )); then echo "[SKIP_SIZE] $acc bp=$bp"; rm -f "$z"; continue; fi
    sh=$(sha256sum "$z" | awk '{print $1}')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$STR" "$acc" "$bp" "$nr" "$bytes" "$sh" "$uri" >> "$TMP"
    got=$((got+bp)); echo "[ACCEPT] $STR $acc bp=$bp cumulative=$got/${TARGET[$STR]}"
  done
  if (( got < TARGET[$STR] )); then echo "[FATAL] stratum $STR insufficient: $got < ${TARGET[$STR]}; rerun selector with larger K" >&2; exit 2; fi
done
mv "$TMP" "$OUT"
python - "$OUT" "$OLD" "$S/frozen/F3V2_DELTA_SOURCE_FREEZE.json" <<'PY'
import csv,json,sys,hashlib
p,oldp,out=sys.argv[1:]
r=list(csv.DictReader(open(p),delimiter='\t')); old=set(x.strip() for x in open(oldp) if x.strip())
assert len({x['accession'] for x in r})==len(r)
assert not ({x['accession'] for x in r}&old)
bp=sum(int(x['nucleotide_bp']) for x in r); by={s:sum(int(x['nucleotide_bp']) for x in r if x['stratum']==s) for s in 'ABCD'}
sha=lambda p:hashlib.sha256(open(p,'rb').read()).hexdigest()
g={'stage':'Stage05F3V2-source-freeze','gate':'PASS_F3V2_54G_DELTA_SOURCE_BLIND_FREEZE','gate_pass':50000000000<=bp<=60000000000,
'delta_accession_n':len(r),'delta_nucleotide_bp':bp,'stratum_bp':by,'old_accession_overlap_n':0,
'source_selection_used_rt_or_candidate_information':False,'selection_inputs':'SRA accession/assay/organism/taxonomy plus Logan technical file size and realized nucleotide bp only',
'contigs_primary':True,'unitigs_used':False,'selected_manifest_sha256':sha(p),
'claim_boundary':'This freezes a post-discovery external delta universe before any F3V2 HMM result is inspected; it does not modify the original 6.12-Gb universe.'}
open(out,'w').write(json.dumps(g,indent=2)+'\n'); print(json.dumps(g,indent=2)); assert g['gate_pass']
PY
chmod a-w "$OUT" "$S/frozen/F3V2_DELTA_SOURCE_FREEZE.json"
echo PASS_F3V2_54G_DELTA_DOWNLOAD_AND_SOURCE_FREEZE_COMPLETE
