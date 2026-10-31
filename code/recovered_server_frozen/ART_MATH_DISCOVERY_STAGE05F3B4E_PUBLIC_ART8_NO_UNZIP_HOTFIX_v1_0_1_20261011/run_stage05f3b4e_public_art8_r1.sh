#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
B="$P/ART_MATH_DISCOVERY_STAGE05F3B4E_GLOBAL_NOVELTY_AND_PUBLIC_ART_CROSSWALK_v1_0_0_20261011"
for x in datasets makeblastdb tblastn blastn python sha256sum; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x" >&2; exit 2; }; done
for f in "$B/inputs/PUBLIC_ART8_ACCESSIONS.txt" "$B/inputs/PUBLIC_ART8_MANIFEST.tsv" "$B/inputs/F3B4E_RT5.faa" "$B/inputs/F3B4E_PARTNER5.faa" "$B/inputs/F3B4E_ARCHITECTURE5.fna" "$B/scripts/summarize_art8.py" "$B/scripts/common.py"; do [[ -s "$f" ]] || { echo "[FATAL] missing $f; run base prepare first" >&2; exit 2; }; done
# Freeze implementation-only transport/extraction hotfix before ART8 results.
python - "$B" "$S/frozen/F3B4E_PUBLIC_ART8_R1_LOCK.json" "$S/run_stage05f3b4e_public_art8_r1.sh" <<'PY'
import sys,json,hashlib,os
B,out,run=sys.argv[1:]
def h(p):
 x=hashlib.sha256();
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):x.update(b)
 return x.hexdigest()
g={'stage':'Stage05F3B4E-G2-R1','gate':'PASS_F3B4E_PUBLIC_ART8_NO_UNZIP_HOTFIX_FROZEN_PRE_RESULTS','gate_pass':True,
'hotfix_role':'IMPLEMENTATION_EXTRACTION_DEPENDENCY_ONLY','base_stage':B,
'previous_failure':'base public-ART8 runner aborted at preflight because external unzip executable was absent; no ART8 BLAST result was generated',
'change':'replace external unzip with Python stdlib zipfile extraction; datasets request, ART8 accession list, BLAST programs/parameters and summarizer remain unchanged',
'scientific_parameters_changed':False,'target_membership_changed':False,'public_art8_membership_changed':False,
'base_manifest_sha256':h(B+'/inputs/PUBLIC_ART8_MANIFEST.tsv'),'rt5_sha256':h(B+'/inputs/F3B4E_RT5.faa'),'partner5_sha256':h(B+'/inputs/F3B4E_PARTNER5.faa'),'architecture5_sha256':h(B+'/inputs/F3B4E_ARCHITECTURE5.fna'),'runner_sha256':h(run),
'claim_boundary':'Implementation-only extraction hotfix. Public ART8 remains a named eight-genome subset, not the full ART95 catalogue.'}
open(out,'w').write(json.dumps(g,indent=2)+'\n');print(json.dumps(g,indent=2))
PY
mkdir -p "$B"/{outputs,logs,work/art8,tmp/art8db}
ZIP="$B/work/art8/public_art8_ncbi_dataset.zip"; U="$B/work/art8/unpacked_r1"
if [[ ! -s "$ZIP" ]]; then
 datasets download virus genome accession --inputfile "$B/inputs/PUBLIC_ART8_ACCESSIONS.txt" --include genome,protein,annotation --filename "$ZIP" --no-progressbar 2>&1 | tee "$B/logs/art8_datasets_download_r1.log"
fi
rm -rf "$U"; mkdir -p "$U"
python - "$ZIP" "$U" <<'PY'
import zipfile,sys,os
z,u=sys.argv[1:]
with zipfile.ZipFile(z) as f:
 f.extractall(u)
print('[PASS] Python zipfile extraction')
PY
G=$(find "$U" -type f -path '*/ncbi_dataset/data/genomic.fna' -print | head -1)
[[ -n "$G" && -s "$G" ]] || { echo '[FATAL] genomic.fna not found in ART8 package' >&2; exit 2; }
# Virus data package contains one combined genomic.fna for the requested accession set. Assert all 8 requested accession headers are present.
python - "$G" "$B/inputs/PUBLIC_ART8_ACCESSIONS.txt" <<'PY'
import re,sys
p,a=sys.argv[1:]; want=[x.strip() for x in open(a) if x.strip()]; heads=[]
with open(p,errors='ignore') as f:
 for line in f:
  if line.startswith('>'): heads.append(line[1:].split()[0])
miss=[x for x in want if not any(h==x or h.split('.')[0]==x.split('.')[0] for h in heads)]
print('requested_accession_n =',len(want)); print('genomic_fasta_record_n =',len(heads)); print('missing_accession_n =',len(miss))
if miss: raise SystemExit('FATAL missing ART8 accessions in genomic.fna: '+','.join(miss))
PY
DB="$B/tmp/art8db/ART8_R1"; rm -f "$DB"* 2>/dev/null || true
makeblastdb -in "$G" -dbtype nucl -parse_seqids -out "$DB" > "$B/logs/art8_makeblastdb_r1.log" 2>&1
OF='6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen qcovs'
tblastn -query "$B/inputs/F3B4E_RT5.faa" -db "$DB" -evalue 1e-5 -max_target_seqs 1000 -seg yes -outfmt "$OF" > "$B/outputs/F3B4E_RT5_VS_PUBLIC_ART8_TBLASTN.tsv"
tblastn -query "$B/inputs/F3B4E_PARTNER5.faa" -db "$DB" -evalue 1e-3 -max_target_seqs 1000 -seg yes -outfmt "$OF" > "$B/outputs/F3B4E_PARTNER5_VS_PUBLIC_ART8_TBLASTN.tsv"
blastn -task dc-megablast -query "$B/inputs/F3B4E_ARCHITECTURE5.fna" -db "$DB" -evalue 1e-10 -max_target_seqs 1000 -dust yes -outfmt "$OF" > "$B/outputs/F3B4E_ARCHITECTURE5_VS_PUBLIC_ART8_BLASTN.tsv"
PYTHONPATH="$B/scripts${PYTHONPATH:+:$PYTHONPATH}" python "$B/scripts/summarize_art8.py" "$B/outputs/F3B4E_TARGETS_5.tsv" "$B/inputs/PUBLIC_ART8_MANIFEST.tsv" "$B/outputs/F3B4E_RT5_VS_PUBLIC_ART8_TBLASTN.tsv" "$B/outputs/F3B4E_PARTNER5_VS_PUBLIC_ART8_TBLASTN.tsv" "$B/outputs/F3B4E_ARCHITECTURE5_VS_PUBLIC_ART8_BLASTN.tsv" "$G" "$B/outputs/F3B4E_PUBLIC_ART8_DETAIL_40.tsv" "$B/frozen/F3B4E_IMPLEMENTATION_LOCK.json" | tee "$B/logs/summarize_art8_r1.log"
echo PASS_STAGE05F3B4E_PUBLIC_ART8_R1_NO_UNZIP_COMPLETE
