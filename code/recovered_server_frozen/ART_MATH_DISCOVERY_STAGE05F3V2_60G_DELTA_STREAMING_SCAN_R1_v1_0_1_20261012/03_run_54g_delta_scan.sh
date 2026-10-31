#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd); P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
CPU=${1:-12}; JOBS=${2:-4}; SEL="$S/frozen/F3V2_DELTA_SELECTED_ACCESSIONS.tsv"
[[ -s "$SEL" && -s "$S/frozen/F3V2_DELTA_SOURCE_FREEZE.json" ]] || { echo '[FATAL] source universe not frozen' >&2; exit 2; }
for x in zstd hmmscan hmmpress mmseqs python; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x" >&2; exit 2; }; done
EXPECTED=74556bd21fd92b947af2759878c3433361b442f2499a62492166a5416971e25e
if [[ -n "${RVT_HMM:-}" ]]; then HMM="$RVT_HMM"; else
  mapfile -t HITS < <(find "$P" -type f -name 'RVT-All.hmm' -print 2>/dev/null)
  HMM=''
  for h in "${HITS[@]}"; do [[ $(sha256sum "$h"|awk '{print $1}') == "$EXPECTED" ]] && { HMM="$h"; break; }; done
fi
[[ -n "$HMM" && -s "$HMM" ]] || { echo '[FATAL] frozen RVT-All.hmm with expected SHA not found; export RVT_HMM=/path/RVT-All.hmm' >&2; exit 2; }
[[ $(sha256sum "$HMM"|awk '{print $1}') == "$EXPECTED" ]] || { echo '[FATAL] RVT-All.hmm SHA mismatch' >&2; exit 2; }
# Freeze execution/science contract before results.
python - "$SEL" "$HMM" "$S/frozen/F3V2_SCAN_IMPLEMENTATION_LOCK.json" "$CPU" "$JOBS" <<'PY'
import sys,json,hashlib
sel,hmm,out,cpu,jobs=sys.argv[1:]
sha=lambda p:hashlib.sha256(open(p,'rb').read()).hexdigest()
g={'stage':'Stage05F3V2-scan-freeze','gate':'PASS_F3V2_DELTA_SCAN_RULES_FROZEN_PRE_RESULTS','gate_pass':True,
'selected_manifest_sha256':sha(sel),'rvt_hmm_sha256':sha(hmm),'rvt_hmm_expected_sha256':'74556bd21fd92b947af2759878c3433361b442f2499a62492166a5416971e25e',
'translation':'six-frame maximal stop-to-stop; genetic code 11; no start-codon requirement; min 45 aa','formal_full_evalue':1e-4,'formal_domain_i_evalue':1e-4,'exploratory_full_evalue':1e-2,'exploratory_domain_i_evalue':1e-2,
'family_min_seq_id':0.30,'family_coverage':0.70,'family_cov_mode':0,'family_cluster_mode':0,'family_sensitivity':7.5,
'cpu_per_hmmscan_job':int(cpu),'parallel_accession_jobs':int(jobs),'annotation_used_for_candidate_generation':False,'old_6p12gb_recomputed':False,
'claim_boundary':'Resource sharding may change runtime only; source universe and frozen F3 scientific rules are fixed before delta scan results.'}
open(out,'w').write(json.dumps(g,indent=2)+'\n');print(json.dumps(g,indent=2))
PY
mkdir -p "$S"/{work/per_accession,work/done,outputs/formal_faa,logs/per_accession,tmp/mmseqs}
cut -f2 "$SEL" | tail -n +2 > "$S/work/accessions_to_scan.txt"
export KEEP_ORF_FASTA=${KEEP_ORF_FASTA:-0}
export S HMM CPU
cat "$S/work/accessions_to_scan.txt" | xargs -n1 -P "$JOBS" -I{} bash "$S/scripts/scan_one_accession.sh" "$S" {} "$HMM" "$CPU"
# Combine candidate tables and formal FASTA.
python - "$S" <<'PY'
import csv,glob,os,sys
S=sys.argv[1]; files=sorted(glob.glob(S+'/work/per_accession/*/candidates.tsv'))
out=S+'/outputs/F3V2_FORMAL_INSTANCES.tsv'; rows=[]
for p in files:
 acc=os.path.basename(os.path.dirname(p))
 for r in csv.DictReader(open(p),delimiter='\t'):
  r={'source_accession':acc,**r}; rows.append(r)
fields=['source_accession','query_id','best_rvt_target','full_evalue','best_domain_i_evalue','domain_score','formal','exploratory']
with open(out,'w',newline='') as f:
 w=csv.DictWriter(f,delimiter='\t',fieldnames=fields);w.writeheader();w.writerows(rows)
print('exploratory_query_n =',len(rows));print('formal_query_n =',sum(x['formal']=='True' for x in rows))
PY
cat "$S"/outputs/formal_faa/*.formal.faa > "$S/outputs/F3V2_FORMAL_ALL.faa" || true
python "$S/scripts/dedupe_formal.py" "$S/outputs/F3V2_FORMAL_ALL.faa" "$S/outputs/F3V2_FORMAL_EXACT_UNIQUE.faa" "$S/outputs/F3V2_FORMAL_EXACT_MAP.tsv"
rm -rf "$S/tmp/mmseqs" "$S/outputs/F3V2_FAMILY_"* 2>/dev/null || true
if [[ -s "$S/outputs/F3V2_FORMAL_EXACT_UNIQUE.faa" ]]; then
  mmseqs easy-cluster "$S/outputs/F3V2_FORMAL_EXACT_UNIQUE.faa" "$S/outputs/F3V2_FAMILY" "$S/tmp/mmseqs" --min-seq-id 0.30 -c 0.70 --cov-mode 0 --cluster-mode 0 -s 7.5 --threads "$((CPU*JOBS))" > "$S/logs/mmseqs_family.log" 2>&1
else
  : > "$S/outputs/F3V2_FAMILY_cluster.tsv"
fi
python "$S/scripts/finalize_scan.py" "$S" | tee "$S/logs/finalize_scan.log"
echo PASS_STAGE05F3V2_54G_DELTA_SCAN_COMPLETE
