#!/usr/bin/env bash
set -euo pipefail
THREADS=${1:-32}
P=${ART_MATH_DISCOVERY_ROOT:-/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$HERE"/{outputs,logs,frozen,work,tmp}
R4="$P/ART_MATH_DISCOVERY_STAGE05F3B2_R4_SOURCE_ACCESSION_METADATA_RESCUE_v1_0_0_20261011/outputs"
FM="$R4/F3B2_R4_CANDIDATE_FAMILY_MAP_WITH_SOURCE.tsv"
IM="$R4/F3B1R3_INSTANCE_MAP_WITH_CANONICAL_SOURCE.tsv"
FS="$R4/F3B2_R4_FAMILY_SUMMARY_WITH_SOURCE.tsv"
R4G="$R4/STAGE05F3B2_R4_GATE.json"
S0="$P/ART_MATH_DISCOVERY_STAGE05F3B_DE_NOVO_ANNOTATION_INDEPENDENT_RT_DISCOVERY_v1_0_0_20261010/work"
ORFDB="$S0/F3B1_orfs_aa"
ORFI="$S0/F3B1_orfs_aa.index"
GDB="$S0/genomeDB"
GIDX="$S0/genomeDB.index"
VREP="$HERE/scripts/vendor_repeat_array_scan.py"
VPLC="$HERE/scripts/vendor_family_wide_context_rescue.py"
LOCK="$HERE/frozen/F3B4C_IMPLEMENTATION_LOCK.json"
for x in python mmseqs sha256sum; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x" >&2; exit 2; }; done
for f in "$FM" "$IM" "$FS" "$R4G" "$ORFDB" "$ORFI" "$GDB" "$GIDX" "$VREP" "$VPLC"; do [[ -e "$f" ]] || { echo "[FATAL] missing $f" >&2; exit 2; }; done

"$HERE/tests/run_selftest.sh"

echo '[1/6] freeze global 10,598-locus repeat census AND partner rules before any new census result'
python - "$HERE" "$FM" "$IM" "$FS" "$R4G" "$LOCK" <<'PY'
import csv,json,hashlib,sys
from pathlib import Path
here,fm,im,fs,r4g,out=map(Path,sys.argv[1:])
def sha(p):
 h=hashlib.sha256();
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(16<<20),b''): h.update(b)
 return h.hexdigest()
g=json.load(open(r4g)); assert g['gate_pass'] is True
assert g['formal_candidate_instance_n']==10598 and g['sequence_family_n']==4912
assert sha(fm)=='859ea2ef382668853e413837c98f235dada60810ea83ea335d30a316c5185406'
assert sha(im)=='6881ebb12c558b8713a0596a3e412ce6ed375691a17ee01002c88e18617538cf'
assert sha(here/'scripts/vendor_repeat_array_scan.py')=='9e3c5b79f807a273e715eff4d01ddd45fdce96cb8047aefdd077ba425e44174f'
assert sha(here/'scripts/vendor_family_wide_context_rescue.py')=='fe9dd9ca8190bdb3d066b33f5ba36de4c818662aec6592a7bd1471e081b6034f'
lock={
 'stage':'Stage05F3B4C-implementation-freeze','gate':'PASS_F3B4C_GLOBAL_REPEAT_AND_PARTNER_RULES_FROZEN_PRE_RESULTS','gate_pass':True,
 'formal_rt_locus_n':10598,'formal_rt_family_n':4912,
 'repeat_census_rule':'all and only frozen F3B2-R4 formal RVT-supported loci; no family/structure/repeat preselection',
 'repeat_detector':'byte-identical F3B4B1 detector','repeat_seed_k_values':[8,10,12,14,16],'repeat_min_copies':3,'repeat_start_gap_bp':[60,360],'repeat_max_exact_core_bp':80,
 'strict_published_art_geometry':'UPSTREAM; 3-21 copies; every inferred spacer 120-220 bp','broad_art_like_geometry':'UPSTREAM; 3-30 copies; median spacer 100-240 bp; >=75% spacers 80-280 bp',
 'placement_rule':'byte-identical F3B4B0 six-frame exact peptide replay; unique exact loci are architecture-scannable; no manual tie-break',
 'partner_scan_rule':'for every uniquely placed formal RT locus, scan up to 10 kb physically downstream after the RT terminal stop; enumerate complete code-11 ORFs >=30 aa on both strands; nearest physical ORF is the primary partner candidate',
 'partner_start_codons':['ATG','GTG','TTG'],'partner_stop_codons':['TAA','TAG','TGA'],'partner_min_complete_orf_aa':30,'partner_scan_bp':10000,
 'partner_cluster_rule':'MMseqs2 easy-cluster; min-seq-id 0.30; coverage 0.70; cov-mode 0; cluster-mode 0; sensitivity 7.5',
 'recurrent_partner_rule':'same nearest-downstream ORF cluster in >=2 independent accessions among loci carrying the indicated repeat geometry',
 'partner_function_confirmation_performed':False,'dedicated_partner_gene_confirmation_performed':False,'novelty_decision_performed':False,
 'r4_family_map_sha256':sha(fm),'r4_instance_map_sha256':sha(im),'r4_family_summary_sha256':sha(fs),'r4_gate_sha256':sha(r4g),
 'vendored_repeat_detector_sha256':sha(here/'scripts/vendor_repeat_array_scan.py'),'vendored_placement_engine_sha256':sha(here/'scripts/vendor_family_wide_context_rescue.py'),
 'global_census_script_sha256':sha(here/'scripts/global_formal_repeat_census.py'),'partner_orf_script_sha256':sha(here/'scripts/partner_orf_scan.py'),'partner_summary_script_sha256':sha(here/'scripts/summarize_partner_architecture.py'),
 'claim_boundary':'Exhaustive computational architecture census. Repeat geometry plus recurrent adjacent ORF cluster supports ART-like architecture but does not establish biochemical activity, dedicated interaction, function, or novelty.'}
with open(out,'w') as f: json.dump(lock,f,indent=2); f.write('\n')
print(json.dumps(lock,indent=2))
PY

echo '[2/6] exhaustive frozen repeat-array census across all 10,598 formal RT loci'
python "$HERE/scripts/global_formal_repeat_census.py" \
  "$FM" "$IM" "$ORFDB" "$ORFI" "$GDB" "$GIDX" \
  "$HERE/outputs" "$LOCK" "$VPLC" "$VREP" "$FS" \
  2>&1 | tee "$HERE/logs/f3b4c1_global_repeat_census.log"

echo '[3/6] enumerate downstream partner-ORF candidates for every uniquely placed formal RT locus'
python "$HERE/scripts/partner_orf_scan.py" \
  "$HERE/outputs/F3B4C1_GLOBAL_FORMAL_RT_LOCUS_REPEAT_CENSUS_10598.tsv" \
  "$GDB" "$GIDX" "$HERE/outputs" "$LOCK" "$VPLC" "$VREP" \
  2>&1 | tee "$HERE/logs/f3b4c2_partner_orf_scan.log"

echo '[4/6] operationally cluster nearest downstream ORF proteins; same 30% identity / 70% coverage compression rule'
FA="$HERE/outputs/F3B4C2_NEAREST_DOWNSTREAM_ORF_PROTEINS.faa"
PREFIX="$HERE/work/F3B4C2_NEAREST_PARTNER"
TMP="$HERE/tmp/mmseqs_partner"
rm -rf "$TMP" "$PREFIX"* 2>/dev/null || true
if [[ -s "$FA" ]]; then
  mmseqs easy-cluster "$FA" "$PREFIX" "$TMP" \
    --threads "$THREADS" \
    --min-seq-id 0.30 \
    -c 0.70 \
    --cov-mode 0 \
    --cluster-mode 0 \
    -s 7.5 \
    > "$HERE/logs/mmseqs_partner_cluster.log" 2>&1
  [[ -s "${PREFIX}_cluster.tsv" ]] || { echo '[FATAL] mmseqs partner cluster table missing' >&2; exit 2; }
else
  : > "${PREFIX}_cluster.tsv"
fi
cp -f "${PREFIX}_cluster.tsv" "$HERE/outputs/F3B4C2_NEAREST_PARTNER_OPERATIONAL_CLUSTER.tsv"

echo '[5/6] integrate repeat geometry + nearest downstream partner-cluster recurrence'
python "$HERE/scripts/summarize_partner_architecture.py" \
  "$HERE/outputs/F3B4C1_GLOBAL_FORMAL_RT_LOCUS_REPEAT_CENSUS_10598.tsv" \
  "$HERE/outputs/F3B4C1_ALL_DETECTED_REPEAT_ARRAYS.tsv" \
  "$HERE/outputs/F3B4C2_NEAREST_DOWNSTREAM_ORF_PER_RT_LOCUS.tsv" \
  "$HERE/outputs/F3B4C2_NEAREST_PARTNER_OPERATIONAL_CLUSTER.tsv" \
  "$HERE/outputs/F3B4C1_GLOBAL_FORMAL_RT_FAMILY_REPEAT_SUMMARY_4912.tsv" \
  "$HERE/outputs" "$LOCK" \
  2>&1 | tee "$HERE/logs/f3b4c2_partner_architecture_summary.log"

echo '[6/6] hashes and compact closeout'
(cd "$HERE" && sha256sum frozen/F3B4C_IMPLEMENTATION_LOCK.json outputs/* > outputs/STAGE05F3B4C_SHA256.tsv)
echo
cat "$HERE/outputs/STAGE05F3B4C1_GATE.json"
echo
cat "$HERE/outputs/STAGE05F3B4C2_GATE.json"
echo
printf '%s\n' '=== TOP ARRAY + RECURRENT PARTNER FAMILIES ==='
python - "$HERE/outputs/F3B4C2_ARRAY_ASSOCIATED_RT_PARTNER_CANDIDATE_FAMILIES.tsv" <<'PY'
import csv,sys
r=list(csv.DictReader(open(sys.argv[1]),delimiter='\t'))
for x in r[:40]:
 print('\t'.join(str(x[k]) for k in ['family_rep_id','strict_array_locus_n','strict_array_dominant_partner_cluster_accession_n','strict_array_dominant_partner_cluster_median_gap_bp','broad_array_locus_n','broad_array_dominant_partner_cluster_accession_n','upstream_array_locus_n','dominant_partner_cluster_accession_n']))
PY
echo PASS_STAGE05F3B4C_GLOBAL_FORMAL_RT_ARCHITECTURE_CENSUS_AND_PARTNER_SCAN_RUN_COMPLETE
