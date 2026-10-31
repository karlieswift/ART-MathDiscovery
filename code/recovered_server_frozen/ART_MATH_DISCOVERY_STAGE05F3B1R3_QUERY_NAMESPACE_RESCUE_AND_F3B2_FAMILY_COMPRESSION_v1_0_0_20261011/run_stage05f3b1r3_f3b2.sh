#!/usr/bin/env bash
set -euo pipefail
THREADS=${1:-32}
S=$(cd "$(dirname "$0")" && pwd)
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1

S0="$P/ART_MATH_DISCOVERY_STAGE05F3B_DE_NOVO_ANNOTATION_INDEPENDENT_RT_DISCOVERY_v1_0_0_20261010"
S1="$P/ART_MATH_DISCOVERY_STAGE05F3B_DE_NOVO_ANNOTATION_INDEPENDENT_RT_DISCOVERY_v1_0_1_20261011"
OLD="$P/ART_MATH_DISCOVERY_STAGE05F3B2_FORMAL_CANDIDATE_FREEZE_AND_FAMILY_COMPRESSION_v1_0_0_20261011"

ORFDB="$S0/work/F3B1_orfs_aa"
ORFIDX="$ORFDB.index"
ORFDBH="${ORFDB}_h"
ORFHDX="${ORFDB}_h.index"
HMM_INPUT="$S1/work/F3B1_ALL_ORFS_NUMERIC.faa"
HMM="$S1/references/RVT-All.hmm"
FROZEN_QNAMES="$S1/outputs/STAGE05F3B1_FORMAL_IDS.txt"
F3B1_GATE="$S1/outputs/STAGE05F3B1_GATE.json"
F3B1_AUDIT="$S1/outputs/STAGE05F3B1_SHARDING_POSTRUN_AUDIT.json"
RUNTIME_LOCK="$S1/frozen/F3B0_V101_IMPLEMENTATION_LOCK.json"
OLD_LOCK="$OLD/frozen/F3B2_IMPLEMENTATION_LOCK.json"
PROTOCOL="$S1/docs/F3B_FROZEN_PROTOCOL.md"

EXTRACT="$S/scripts/extract_namespace_instances.py"
COLLECT="$S/scripts/collect_instance_hmm.py"
ANALYSIS="$S/scripts/stage05f3b2_dynamic.py"
LOCK="$S/frozen/F3B1R3_F3B2_IMPLEMENTATION_LOCK.json"

mkdir -p "$S"/{frozen,raw,work,tmp,outputs,logs,docs}
for x in python hmmscan mmseqs sha256sum; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x" >&2; exit 2; }; done
for x in "$ORFDB" "$ORFIDX" "$ORFDBH" "$ORFHDX" "$HMM_INPUT" "$HMM" "$FROZEN_QNAMES" "$F3B1_GATE" "$F3B1_AUDIT" "$RUNTIME_LOCK" "$OLD_LOCK" "$PROTOCOL" "$EXTRACT" "$COLLECT" "$ANALYSIS"; do [[ -e "$x" ]] || { echo "[FATAL] missing input $x" >&2; exit 2; }; done

HMM_SHA=$(sha256sum "$HMM"|awk '{print $1}')
PROTOCOL_SHA=$(sha256sum "$PROTOCOL"|awk '{print $1}')

if [[ ! -s "$LOCK" ]]; then
  echo "[1/8] freeze namespace-rescue implementation before corrected instance HMM results"
  python - "$FROZEN_QNAMES" "$F3B1_GATE" "$F3B1_AUDIT" "$RUNTIME_LOCK" "$OLD_LOCK" "$HMM_INPUT" "$HMM" "$PROTOCOL" "$EXTRACT" "$COLLECT" "$ANALYSIS" "$LOCK" <<'PY'
import hashlib,json,sys
(qnames,g1,a1,r1,old,hmm_input,hmm,protocol,extract,collect,analysis,out)=sys.argv[1:13]
def h(p):
    x=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''):x.update(b)
    return x.hexdigest()
ids=[x.strip() for x in open(qnames) if x.strip()]
g=json.load(open(g1));a=json.load(open(a1));r=json.load(open(r1));o=json.load(open(old))
assert len(ids)==10464 and len(set(ids))==10464
assert g['formal_candidate_n']==10464 and g['gate_pass'] is True
assert a['sum_shard_record_n']==65315646 and a['gate_pass'] is True
assert r['orf_instance_n']==65315646
for k,v in {'family_min_seq_id':0.30,'family_coverage':0.70,'family_cov_mode':0,'family_cluster_mode':0,'family_sensitivity':7.5}.items(): assert o[k]==v,(k,o[k],v)
lock={'stage':'Stage05F3B1R3+F3B2','gate':'PASS_F3B1R3_NAMESPACE_RESCUE_AND_F3B2_IMPLEMENTATION_FROZEN_PRE_RERUN','gate_pass':True,'failure_discovered':'F3B1 HMM-input FASTA query-name token is not unique; duplicate formal query name 54 proves old 10,464 count is a namespace-collapsed name count rather than a guaranteed peptide-instance count','rescue_rule':'expand every frozen formal query name to every matching HMM-input FASTA record, assign a unique ordinal-derived instance ID, validate against original ORF DB, rerun only these instances with the unchanged frozen RVT-All HMM and thresholds','full_65m_hmm_rerun_required':False,'old_formal_query_name_n':10464,'old_formal_query_names_sha256':h(qnames),'hmm_input_fasta_sha256':h(hmm_input),'rvt_hmm_sha256':h(hmm),'protocol_sha256':h(protocol),'formal_full_E':1e-4,'formal_best_domain_i_E':1e-4,'exploratory_full_E':1e-2,'exploratory_best_domain_i_E':1e-2,'scientific_threshold_changed':False,'annotation_filter_added':False,'family_min_seq_id':0.30,'family_coverage':0.70,'family_cov_mode':0,'family_cluster_mode':0,'family_sensitivity':7.5,'family_parameters_changed':False,'extract_script_sha256':h(extract),'collector_script_sha256':h(collect),'analysis_script_sha256':h(analysis),'old_f3b1_gate_sha256':h(g1),'old_f3b1_sharding_audit_sha256':h(a1),'old_f3b1_runtime_lock_sha256':h(r1),'old_f3b2_lock_sha256':h(old),'claim_boundary':'This is an identifier-namespace rescue. Candidate-generation HMM, thresholds, source universe, and F3B2 family parameters are unchanged.'}
with open(out,'w') as f:json.dump(lock,f,indent=2);f.write('\n')
print(json.dumps(lock,indent=2))
PY
  sha256sum "$LOCK" "$EXTRACT" "$COLLECT" "$ANALYSIS" > "$S/frozen/F3B1R3_F3B2_FREEZE_SHA256.tsv"
  chmod a-w "$LOCK" "$EXTRACT" "$COLLECT" "$ANALYSIS"
else
  echo "[1/8] [SKIP] implementation lock already present"
fi

echo "[2/8] expand the frozen 10,464 query-name namespace to unique peptide instances"
python "$EXTRACT" "$FROZEN_QNAMES" "$HMM_INPUT" "$ORFDB" "$ORFIDX" "$ORFDBH" "$ORFHDX" \
  "$S/raw/F3B1R3_NAMESPACE_INSTANCES.faa" \
  "$S/raw/F3B1R3_INSTANCE_MAP.tsv" \
  "$S/outputs/F3B1R3_NAMESPACE_EXTRACTION_AUDIT.json" \
  "$S/outputs/F3B1R3_FORMAL_QNAME_MULTIPLICITY.tsv" \
  | tee "$S/logs/f3b1r3_namespace_extract.log"

echo "[3/8] rerun ONLY the expanded formal-name instances with unchanged frozen RVT-All screen"
hmmscan --cpu "$THREADS" --noali -o /dev/null -E 1e-2 --domE 1000 \
  --tblout "$S/raw/F3B1R3_RVT_ALL.tblout" \
  --domtblout "$S/raw/F3B1R3_RVT_ALL.domtblout" \
  "$HMM" "$S/raw/F3B1R3_NAMESPACE_INSTANCES.faa" \
  2> "$S/logs/f3b1r3_hmmscan.stderr"

echo "[4/8] corrected instance-level formal candidate freeze"
python "$COLLECT" \
  "$S/raw/F3B1R3_RVT_ALL.tblout" "$S/raw/F3B1R3_RVT_ALL.domtblout" \
  "$S/raw/F3B1R3_INSTANCE_MAP.tsv" "$FROZEN_QNAMES" \
  "$S/outputs/F3B1R3_NAMESPACE_EXTRACTION_AUDIT.json" "$S/outputs" "$HMM_SHA" "$PROTOCOL_SHA" \
  | tee "$S/logs/f3b1r3_collect.log"

echo "[5/8] exact-AA compression of corrected formal instances"
python "$ANALYSIS" prepare \
  "$S/outputs/STAGE05F3B1R3_FORMAL_INSTANCE_IDS.txt" \
  "$S/raw/F3B1R3_NAMESPACE_INSTANCES.faa" \
  "$S/raw/F3B1R3_INSTANCE_MAP.tsv" \
  "$S/outputs/STAGE05F3B1R3_INSTANCE_RVT_PAIR_HITS.tsv.gz" \
  "$S/outputs" | tee "$S/logs/f3b2_prepare.log"

echo "[6/8] frozen operational sequence-family clustering"
rm -rf "$S/work/F3B2_exactDB"* "$S/work/F3B2_familyDB"* "$S/tmp/F3B2_family"
mmseqs createdb "$S/outputs/F3B2_EXACT_AA_REPRESENTATIVES.faa" "$S/work/F3B2_exactDB" --shuffle 0 > "$S/logs/f3b2_createdb.log" 2>&1
mmseqs cluster "$S/work/F3B2_exactDB" "$S/work/F3B2_familyDB" "$S/tmp/F3B2_family" \
  --min-seq-id 0.30 -c 0.70 --cov-mode 0 --cluster-mode 0 -s 7.5 --threads "$THREADS" > "$S/logs/f3b2_cluster.log" 2>&1
mmseqs createtsv "$S/work/F3B2_exactDB" "$S/work/F3B2_exactDB" "$S/work/F3B2_familyDB" "$S/raw/F3B2_FAMILY_CLUSTER.tsv" > "$S/logs/f3b2_createtsv.log" 2>&1
python "$ANALYSIS" finalize "$S/outputs/F3B2_CANDIDATE_BASE.tsv" "$S/outputs/F3B2_EXACT_AA_MAP.tsv" "$S/raw/F3B2_FAMILY_CLUSTER.tsv" "$LOCK" "$S/outputs" | tee "$S/logs/f3b2_finalize.log"

echo "[7/8] hashes"
sha256sum "$LOCK" "$S/outputs/F3B1R3_NAMESPACE_EXTRACTION_AUDIT.json" "$S/outputs/STAGE05F3B1R3_GATE.json" "$S/outputs/STAGE05F3B1R3_FORMAL_INSTANCE_IDS.txt" "$S/outputs/F3B2_EXACT_AA_MAP.tsv" "$S/outputs/F3B2_CANDIDATE_BASE.tsv" "$S/raw/F3B2_FAMILY_CLUSTER.tsv" "$S/outputs/F3B2_FAMILY_SUMMARY.tsv" "$S/outputs/STAGE05F3B2_R3_GATE.json" > "$S/outputs/STAGE05F3B1R3_F3B2_SHA256.tsv"

echo "[8/8] compact result"
cat "$S/outputs/STAGE05F3B1R3_GATE.json"
echo
cat "$S/outputs/STAGE05F3B2_R3_GATE.json"
echo
echo "=== TOP 20 FAMILIES ==="
column -t -s $'\t' "$S/outputs/F3B2_FAMILY_TOP100.tsv" | head -21 || true
echo
echo "PASS_STAGE05F3B1R3_QUERY_NAMESPACE_RESCUE_AND_F3B2_FAMILY_COMPRESSION_COMPLETE"
