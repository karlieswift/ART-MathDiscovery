#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1

R3="$P/ART_MATH_DISCOVERY_STAGE05F3B1R3_QUERY_NAMESPACE_RESCUE_AND_F3B2_FAMILY_COMPRESSION_v1_0_0_20261011"
S0="$P/ART_MATH_DISCOVERY_STAGE05F3B_DE_NOVO_ANNOTATION_INDEPENDENT_RT_DISCOVERY_v1_0_0_20261010"
MAN="$P/external_db/F3_LARGE_UNIVERSE_20261009/F3A_EXTERNAL_SOURCE_MANIFEST.tsv"

INSTANCE_MAP="$R3/raw/F3B1R3_INSTANCE_MAP.tsv"
FORMAL_IDS="$R3/outputs/STAGE05F3B1R3_FORMAL_INSTANCE_IDS.txt"
BASE="$R3/outputs/F3B2_CANDIDATE_BASE.tsv"
FAMILY_MAP="$R3/outputs/F3B2_CANDIDATE_FAMILY_MAP.tsv"
G1="$R3/outputs/STAGE05F3B1R3_GATE.json"
G2="$R3/outputs/STAGE05F3B2_R3_GATE.json"
GENOME_H="$S0/work/genomeDB_h"
GENOME_H_IDX="$S0/work/genomeDB_h.index"
SCRIPT="$S/scripts/rescue_source_accession.py"
LOCK="$S/frozen/F3B2_R4_IMPLEMENTATION_LOCK.json"

mkdir -p "$S"/{frozen,outputs,logs}

for x in python sha256sum; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x" >&2; exit 2; }; done
for x in "$INSTANCE_MAP" "$FORMAL_IDS" "$BASE" "$FAMILY_MAP" "$G1" "$G2" "$GENOME_H" "$GENOME_H_IDX" "$MAN" "$SCRIPT"; do
  [[ -e "$x" ]] || { echo "[FATAL] missing input: $x" >&2; exit 2; }
done

if [[ ! -s "$LOCK" ]]; then
  echo "[1/4] freeze source-accession metadata rescue implementation"
  python - "$G1" "$G2" "$INSTANCE_MAP" "$FORMAL_IDS" "$BASE" "$FAMILY_MAP" "$GENOME_H" "$GENOME_H_IDX" "$MAN" "$SCRIPT" "$LOCK" <<'PY'
import hashlib,json,sys
(g1,g2,imap,ids,base,fmap,gh,ghi,man,script,out)=sys.argv[1:12]
def h(p):
    x=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''): x.update(b)
    return x.hexdigest()
a=json.load(open(g1)); b=json.load(open(g2))
assert a['gate_pass'] is True and a['corrected_formal_instance_n']==10598
assert b['gate_pass'] is True and b['formal_candidate_instance_n']==10598
assert b['exact_aa_unique_n']==10329 and b['sequence_family_n']==4912
lock={
  'stage':'Stage05F3B2-R4-implementation-freeze',
  'gate':'PASS_F3B2_R4_SOURCE_ACCESSION_METADATA_RESCUE_IMPLEMENTATION_FROZEN',
  'gate_pass':True,
  'failure_discovered':'R3 ORF headers encode parent contig database keys and coordinates, not SRA/ENA run accessions; therefore accession_n was spuriously zero for every family',
  'repair_rule':'parse parent contig key from each R3 ORF header, fetch the original contig header from the canonical F3B1 genomeDB_h, require membership in the frozen 131-accession F3A manifest, patch accession metadata, and recompute descriptive replication summaries without reclustering',
  'candidate_membership_changed':False,
  'hmm_results_changed':False,
  'exact_aa_membership_changed':False,
  'family_membership_changed':False,
  'family_parameters_changed':False,
  'known_rt_annotation_used':False,
  'novelty_screen_performed':False,
  'f3b1r3_gate_sha256':h(g1),
  'f3b2_r3_gate_sha256':h(g2),
  'r3_instance_map_sha256':h(imap),
  'r3_formal_ids_sha256':h(ids),
  'r3_candidate_base_sha256':h(base),
  'r3_candidate_family_map_sha256':h(fmap),
  'genomeDB_h_sha256':h(gh),
  'genomeDB_h_index_sha256':h(ghi),
  'canonical_manifest_sha256':h(man),
  'repair_script_sha256':h(script),
  'claim_boundary':'Metadata/provenance rescue only; no candidate/HMM/family/novelty inference is changed.'
}
with open(out,'w') as f: json.dump(lock,f,indent=2); f.write('\n')
print(json.dumps(lock,indent=2))
PY
  sha256sum "$LOCK" "$SCRIPT" > "$S/frozen/F3B2_R4_FREEZE_SHA256.tsv"
  chmod a-w "$LOCK" "$SCRIPT"
else
  echo "[1/4] [SKIP] R4 implementation lock already present"
fi

echo "[2/4] rescue canonical source accession for all R3 namespace instances and formal candidates"
python "$SCRIPT" \
  "$INSTANCE_MAP" "$FORMAL_IDS" "$BASE" "$FAMILY_MAP" \
  "$G1" "$G2" "$GENOME_H" "$GENOME_H_IDX" "$MAN" \
  "$S/outputs" "$LOCK" \
  | tee "$S/logs/f3b2_r4_rescue.log"

echo "[3/4] hashes"
sha256sum \
  "$LOCK" \
  "$S/outputs/F3B1R3_INSTANCE_MAP_WITH_CANONICAL_SOURCE.tsv" \
  "$S/outputs/F3B2_R4_CANDIDATE_BASE_WITH_SOURCE.tsv" \
  "$S/outputs/F3B2_R4_CANDIDATE_FAMILY_MAP_WITH_SOURCE.tsv" \
  "$S/outputs/F3B2_R4_FAMILY_SUMMARY_WITH_SOURCE.tsv" \
  "$S/outputs/STAGE05F3B2_R4_GATE.json" \
  > "$S/outputs/STAGE05F3B2_R4_SHA256.tsv"

echo "[4/4] compact result"
cat "$S/outputs/STAGE05F3B2_R4_GATE.json"
echo
echo "=== TOP 20 FAMILIES AFTER ACCESSION RESCUE ==="
column -t -s $'\t' "$S/outputs/F3B2_R4_FAMILY_TOP100_WITH_SOURCE.tsv" | head -21 || true
echo
echo "PASS_STAGE05F3B2_R4_SOURCE_ACCESSION_METADATA_RESCUE_COMPLETE"
