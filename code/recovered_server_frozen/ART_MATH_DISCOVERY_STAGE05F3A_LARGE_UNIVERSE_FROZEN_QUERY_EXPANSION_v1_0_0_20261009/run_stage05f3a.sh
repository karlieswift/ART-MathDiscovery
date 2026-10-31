#!/usr/bin/env bash
set -euo pipefail
# Usage: ./run_stage05f3a.sh F2A_ROOT F2E_ROOT PROJECT_ROOT [THREADS]
F2A=${1:?F2A_ROOT}
F2E=${2:?F2E_ROOT}
P=${3:?PROJECT_ROOT}
THREADS=${4:-32}
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
mkdir -p outputs logs work raw external
for x in python mmseqs hmmscan zstdcat sha256sum; do command -v "$x" >/dev/null || { echo "[FATAL] missing tool: $x" >&2; exit 2; }; done
if ! command -v aws >/dev/null && ! command -v wget >/dev/null; then echo '[FATAL] need aws CLI or wget for Logan downloads' >&2; exit 2; fi
Q="$F2A/outputs/STAGE05F2A_NOVELTY_FIRST_STRICT.faa"
QT="$F2A/outputs/STAGE05F2A_NOVELTY_FIRST_STRICT_RESCUES.tsv"
COMP="$F2E/outputs/STAGE05F2E_QUERY_COMPONENT_SUMMARY.tsv"
test -s "$Q";test -s "$QT";test -s "$COMP"
[[ "$(sha256sum "$Q"|awk '{print $1}')" == bfe8bfcf249b28e5a13e04cb52af26a58ee18b708c012b2d154a2056e095fac5 ]]
[[ "$(sha256sum "$QT"|awk '{print $1}')" == cacc0efcb4233b77c0f397f673c7aee107f267472e598cb802b194e2dfe436bf ]]
[[ "$(sha256sum "$COMP"|awk '{print $1}')" == 819493afaf96022b8f9740dbc6fc849e2406f3220f679328e64cfd35f7a9fdfb ]]
# Resolve frozen RVT HMM without using candidate outcomes.
HMM=${STAGE05F3_RVT_HMM:-}
if [[ -z "$HMM" ]]; then
  HMM=$(find "$P" -type f -name '*.hmm' 2>/dev/null | while read -r f; do n=$(grep -c '^HMMER3/f' "$f" 2>/dev/null || true); if [[ "$n" -ge 40 && "$n" -le 60 ]]; then printf '%05d\t%s\n' "$(( n>45 ? n-45 : 45-n ))" "$f"; fi; done | sort | head -1 | cut -f2-)
fi
[[ -n "$HMM" && -s "$HMM" ]] || { echo '[FATAL] cannot resolve frozen ~45-profile RVT HMM; set STAGE05F3_RVT_HMM' >&2; exit 2; }
# ---- F3A0 QUERY-BLIND SOURCE FREEZE ----
META=external/sra_taxid.csv.zst
CANDS=outputs/STAGE05F3A_SOURCE_CANDIDATE_POOL.tsv
if [[ ! -s "$META" ]]; then
  if command -v aws >/dev/null; then aws s3 cp s3://logan-pub/stats/sra_taxid.csv.zst "$META" --no-sign-request --only-show-errors; else wget -c -O "$META" https://s3.amazonaws.com/logan-pub/stats/sra_taxid.csv.zst; fi
fi
sha256sum "$META" | tee outputs/STAGE05F3A_SOURCE_METADATA_SHA256.txt
zstdcat "$META" | python scripts/select_logan_metadata.py --out "$CANDS" --salt ART_STAGE05F3A_20261009_SOURCE_FREEZE_V1 --per-stratum 1000 2> logs/source_metadata_selection.log
EXTROOT="$P/external_db/F3_LARGE_UNIVERSE_20261009"
mkdir -p "$EXTROOT"
EXTFA="$EXTROOT/F3A_EXTERNAL_UNIVERSE.fna"
EXTMAN="$EXTROOT/F3A_EXTERNAL_SOURCE_MANIFEST.tsv"
EXTGATE="$EXTROOT/F3A0_SOURCE_FREEZE_GATE.json"
python scripts/materialize_logan_universe.py --candidates "$CANDS" --external-root "$EXTROOT" --out-fasta "$EXTFA" --out-manifest "$EXTMAN" --out-gate "$EXTGATE" --target-bp 6000000000 --min-bp 5000000000 --max-bp 10000000000 --min-contig-nt 300 --max-accession-bp 250000000 | tee logs/source_materialization.log
cp "$EXTGATE" outputs/STAGE05F3A0_SOURCE_FREEZE_GATE.json
cp "$EXTMAN" outputs/STAGE05F3A0_EXTERNAL_SOURCE_MANIFEST.tsv
sha256sum "$EXTFA" "$EXTMAN" > outputs/STAGE05F3A0_SOURCE_FREEZE_SHA256.tsv
# Critical invariant: source freeze exists before query database is created.
python - <<'PY'
import json
x=json.load(open('outputs/STAGE05F3A0_SOURCE_FREEZE_GATE.json'))
assert x['gate_pass'] and 5_000_000_000 <= x['total_bp'] <= 10_000_000_000,x
print('[PASS] F3A0 source universe frozen before query search')
PY
# ---- F3A1 FROZEN QUERY EXPANSION ----
rm -rf work/queryDB* work/targetDB* work/mmseqs_tmp raw/F3A_QUERY30_VS_EXTERNAL*
mmseqs createdb "$Q" work/queryDB > logs/mmseqs_createdb_query.log 2>&1
mmseqs createdb "$EXTFA" work/targetDB > logs/mmseqs_createdb_target.log 2>&1
mmseqs search work/queryDB work/targetDB raw/F3A_QUERY30_VS_EXTERNAL.db work/mmseqs_tmp --threads "$THREADS" -s 7.5 -e 1e-2 --max-seqs 50000 --mask 0 --orf-start-mode 1 --min-length 45 --translation-table 11 -a > logs/mmseqs_search.log 2>&1
mmseqs convertalis work/queryDB work/targetDB raw/F3A_QUERY30_VS_EXTERNAL.db raw/STAGE05F3A_QUERY30_VS_EXTERNAL.m8 --format-output query,target,evalue,bits,qcov,alnlen,qlen,tlen,fident,qstart,qend,tstart,tend,taln > logs/mmseqs_convertalis.log 2>&1
python scripts/analyze_frozen_query_expansion.py --queries "$Q" --query-table "$QT" --components "$COMP" --external-fasta "$EXTFA" --m8 raw/STAGE05F3A_QUERY30_VS_EXTERNAL.m8 --hmm "$HMM" --threads "$THREADS" --outdir outputs | tee logs/f3a_analysis.log
python - <<'PY'
import json,hashlib
from pathlib import Path
g=json.load(open('outputs/STAGE05F3A_GATE.json'));s=json.load(open('outputs/STAGE05F3A0_SOURCE_FREEZE_GATE.json'))
g.update({'source_universe_bp':s['total_bp'],'source_accession_n':s['accepted_accession_n'],'source_fasta_sha256':s['external_fasta_sha256'],'source_manifest_sha256':s['source_manifest_sha256'],'source_metadata_sha256':Path('outputs/STAGE05F3A_SOURCE_METADATA_SHA256.txt').read_text().split()[0],'f3b_source_lock_ready':True,'f3b_source_fasta_sha256':s['external_fasta_sha256']})
Path('outputs/STAGE05F3A_GATE.json').write_text(json.dumps(g,indent=2)+'\n')
print(json.dumps(g,indent=2))
PY
sha256sum config/F3_MASTER_POLICY.json docs/SCIENTIFIC_BOUNDARY.md docs/F3B_FROZEN_PROTOCOL.md outputs/STAGE05F3A*.json outputs/STAGE05F3A*.tsv raw/STAGE05F3A_QUERY30_VS_EXTERNAL.m8 > outputs/STAGE05F3A_SHA256.tsv
