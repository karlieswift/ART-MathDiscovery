#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$S"/{outputs,logs,work/art8,tmp/art8db}
for x in datasets unzip makeblastdb tblastn blastn python; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x. Suggested: conda install -n art_novelty -y -c conda-forge -c bioconda ncbi-datasets-cli blast" >&2; exit 2; }; done
for f in "$S/inputs/PUBLIC_ART8_ACCESSIONS.txt" "$S/inputs/PUBLIC_ART8_MANIFEST.tsv" "$S/inputs/F3B4E_RT5.faa" "$S/inputs/F3B4E_PARTNER5.faa" "$S/inputs/F3B4E_ARCHITECTURE5.fna"; do [[ -s "$f" ]] || { echo "[FATAL] run prepare first; missing $f" >&2; exit 2; }; done
ZIP="$S/work/art8/public_art8_ncbi_dataset.zip"; U="$S/work/art8/unpacked"
if [[ ! -s "$ZIP" ]]; then
  datasets download virus genome accession --inputfile "$S/inputs/PUBLIC_ART8_ACCESSIONS.txt" --include genome,protein,annotation --filename "$ZIP" --no-progressbar 2>&1 | tee "$S/logs/art8_datasets_download.log"
fi
rm -rf "$U"; mkdir -p "$U"; unzip -q "$ZIP" -d "$U"
G=$(find "$U" -type f -name genomic.fna -print | head -1)
[[ -n "$G" && -s "$G" ]] || { echo '[FATAL] genomic.fna not found in ART8 data package' >&2; exit 2; }
DB="$S/tmp/art8db/ART8"
rm -f "$DB"* 2>/dev/null || true
makeblastdb -in "$G" -dbtype nucl -parse_seqids -out "$DB" > "$S/logs/art8_makeblastdb.log" 2>&1
OF='6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen qcovs'
tblastn -query "$S/inputs/F3B4E_RT5.faa" -db "$DB" -evalue 1e-5 -max_target_seqs 1000 -seg yes -outfmt "$OF" > "$S/outputs/F3B4E_RT5_VS_PUBLIC_ART8_TBLASTN.tsv"
tblastn -query "$S/inputs/F3B4E_PARTNER5.faa" -db "$DB" -evalue 1e-3 -max_target_seqs 1000 -seg yes -outfmt "$OF" > "$S/outputs/F3B4E_PARTNER5_VS_PUBLIC_ART8_TBLASTN.tsv"
blastn -task dc-megablast -query "$S/inputs/F3B4E_ARCHITECTURE5.fna" -db "$DB" -evalue 1e-10 -max_target_seqs 1000 -dust yes -outfmt "$OF" > "$S/outputs/F3B4E_ARCHITECTURE5_VS_PUBLIC_ART8_BLASTN.tsv"
python "$S/scripts/summarize_art8.py" "$S/outputs/F3B4E_TARGETS_5.tsv" "$S/inputs/PUBLIC_ART8_MANIFEST.tsv" "$S/outputs/F3B4E_RT5_VS_PUBLIC_ART8_TBLASTN.tsv" "$S/outputs/F3B4E_PARTNER5_VS_PUBLIC_ART8_TBLASTN.tsv" "$S/outputs/F3B4E_ARCHITECTURE5_VS_PUBLIC_ART8_BLASTN.tsv" "$G" "$S/outputs/F3B4E_PUBLIC_ART8_DETAIL_40.tsv" "$S/frozen/F3B4E_IMPLEMENTATION_LOCK.json" | tee "$S/logs/summarize_art8.log"
echo PASS_STAGE05F3B4E_PUBLIC_ART8_RUN_COMPLETE
