#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
R3="$P/ART_MATH_DISCOVERY_STAGE05F3B1R3_QUERY_NAMESPACE_RESCUE_AND_F3B2_FAMILY_COMPRESSION_v1_0_0_20261011"
R4="$P/ART_MATH_DISCOVERY_STAGE05F3B2_R4_SOURCE_ACCESSION_METADATA_RESCUE_v1_0_0_20261011"
G4="$R4/outputs/STAGE05F3B2_R4_GATE.json"
FMAP="$R4/outputs/F3B2_R4_CANDIDATE_FAMILY_MAP_WITH_SOURCE.tsv"
FSUM="$R4/outputs/F3B2_R4_FAMILY_SUMMARY_WITH_SOURCE.tsv"
FORMAL="$R3/outputs/STAGE05F3B1R3_FORMAL_INSTANCE_IDS.txt"
NAMESPACE="$R3/raw/F3B1R3_NAMESPACE_INSTANCES.faa"
PAIR="$R3/outputs/STAGE05F3B1R3_INSTANCE_RVT_PAIR_HITS.tsv.gz"
DOM="$R3/raw/F3B1R3_RVT_ALL.domtblout"
SCRIPT="$S/scripts/stage05f3b3a.py"
LOCK="$S/frozen/F3B3A_IMPLEMENTATION_LOCK.json"
SHORTLIST_N=32
mkdir -p "$S"/{frozen,outputs,logs}
for x in python sha256sum; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x" >&2; exit 2; }; done
for x in "$G4" "$FMAP" "$FSUM" "$FORMAL" "$NAMESPACE" "$PAIR" "$DOM" "$SCRIPT"; do [[ -s "$x" ]] || { echo "[FATAL] missing input: $x" >&2; exit 2; }; done

if [[ ! -s "$LOCK" ]]; then
  echo "[1/4] freeze F3B3A post-freeze characterization/priority rule before motif scan"
  python - "$G4" "$FMAP" "$FSUM" "$FORMAL" "$NAMESPACE" "$PAIR" "$DOM" "$SCRIPT" "$LOCK" <<'PY'
import hashlib,json,sys
(g4,fmap,fsum,formal,ns,pair,dom,script,out)=sys.argv[1:10]
def h(p):
    x=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''):x.update(b)
    return x.hexdigest()
g=json.load(open(g4))
assert g['gate_pass'] is True
assert g['formal_candidate_instance_n']==10598
assert g['exact_aa_unique_n']==10329
assert g['sequence_family_n']==4912
assert g['formal_accession_n']==98
assert g['multirun_family_n']==733
assert g['max_accession_replication_n']==48
assert g['formal_unknown_accession_n']==0
lock={
 'stage':'Stage05F3B3A-implementation-freeze',
 'gate':'PASS_F3B3A_PRIORITY_RULE_FROZEN_BEFORE_MOTIF_TOPOLOGY_RESULTS',
 'gate_pass':True,
 'role':'post-freeze descriptive characterization and expensive-compute prioritization only',
 'r4_gate_sha256':h(g4),'r4_candidate_family_map_sha256':h(fmap),'r4_family_summary_sha256':h(fsum),
 'formal_ids_sha256':h(formal),'namespace_fasta_sha256':h(ns),'pair_hits_sha256':h(pair),'domtblout_sha256':h(dom),
 'analysis_script_sha256':h(script),
 'motif_c_like_regex':'[YF]xDD',
 'motif_scope':'best formal HMM-domain envelope; whole-sequence presence also reported descriptively',
 'motif_is_identity_threshold':False,
 'structure_priority_shortlist_family_n':32,
 'shortlist_quota_per_axis':8,
 'shortlist_axis_order':['replication','max_orf_length','formal_model_coverage','motif_c_topology'],
 'shortlist_deduplication':'preserve axis order then deterministic backfill to 32 unique families',
 'known_rt_annotation_used_for_ranking':False,
 'novelty_used_for_ranking':False,
 'candidate_membership_changed':False,'family_membership_changed':False,'hmm_results_changed':False,'scientific_threshold_changed':False,
 'chronology':'POST_F3B2_COMPUTE_PRIORITY_IMPLEMENTATION_FREEZE; authorized priority attributes were prespecified in F3B_FROZEN_PROTOCOL.md, but this exact quota/ranking implementation is not claimed as pre-F3A1',
 'claim_boundary':'Priority attributes and motif topology are descriptive/computational triage only and do not establish RT/ART identity or novelty.'
}
with open(out,'w') as f:json.dump(lock,f,indent=2);f.write('\n')
print(json.dumps(lock,indent=2))
PY
  sha256sum "$LOCK" "$SCRIPT" > "$S/frozen/F3B3A_FREEZE_SHA256.tsv"
  chmod a-w "$LOCK" "$SCRIPT"
else
  echo "[1/4] [SKIP] F3B3A implementation lock already present"
fi

echo "[2/4] characterize 10,598 formal instances and 4,912 frozen families"
python "$SCRIPT" "$G4" "$FMAP" "$FSUM" "$FORMAL" "$NAMESPACE" "$PAIR" "$DOM" "$LOCK" "$S/outputs" "$SHORTLIST_N" | tee "$S/logs/f3b3a.log"

echo "[3/4] hashes"
sha256sum "$LOCK" "$S/outputs/F3B3A_CANDIDATE_CHARACTERIZATION.tsv" "$S/outputs/F3B3A_FAMILY_PRIORITY.tsv" "$S/outputs/F3B3A_STRUCTURE_PRIORITY_SHORTLIST_32.tsv" "$S/outputs/F3B3A_STRUCTURE_PRIORITY_SHORTLIST_32.faa" "$S/outputs/STAGE05F3B3A_GATE.json" > "$S/outputs/STAGE05F3B3A_SHA256.tsv"

echo "[4/4] compact result"
cat "$S/outputs/STAGE05F3B3A_GATE.json"
echo
echo "=== STRUCTURE PRIORITY SHORTLIST 32 ==="
column -t -s $'\t' "$S/outputs/F3B3A_STRUCTURE_PRIORITY_SHORTLIST_32.tsv" | head -33 || true
echo
echo "PASS_STAGE05F3B3A_POSTFREEZE_PRIORITY_AND_CATALYTIC_TOPOLOGY_COMPLETE"
