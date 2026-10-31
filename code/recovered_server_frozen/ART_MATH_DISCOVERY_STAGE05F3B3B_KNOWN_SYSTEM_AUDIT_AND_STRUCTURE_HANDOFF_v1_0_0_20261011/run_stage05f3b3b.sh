#!/usr/bin/env bash
set -euo pipefail
THREADS=${1:-32}
S=$(cd "$(dirname "$0")" && pwd)
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1
A="$P/ART_MATH_DISCOVERY_STAGE05F3B3A_POSTFREEZE_PRIORITY_AND_CATALYTIC_TOPOLOGY_v1_0_0_20261011"
G="$A/outputs/STAGE05F3B3A_GATE.json"
QTSV="$A/outputs/F3B3A_STRUCTURE_PRIORITY_SHORTLIST_32.tsv"
QFA="$A/outputs/F3B3A_STRUCTURE_PRIORITY_SHORTLIST_32.faa"
LOCK="$S/frozen/F3B3B_IMPLEMENTATION_LOCK.json"
mkdir -p "$S"/{frozen,outputs,logs,raw,references,work,tmp}
for x in python sha256sum mmseqs; do command -v "$x" >/dev/null || { echo "[FATAL] missing $x" >&2; exit 2; }; done
for x in "$G" "$QTSV" "$QFA" "$S/scripts/build_known_reference_panel.py" "$S/scripts/summarize_known_audit.py"; do [[ -s "$x" ]] || { echo "[FATAL] missing $x" >&2; exit 2; }; done

if [[ ! -s "$LOCK" ]]; then
  echo "[1/6] freeze F3B3B before known-system sequence audit"
  python - "$G" "$QTSV" "$QFA" "$S/scripts/build_known_reference_panel.py" "$S/scripts/summarize_known_audit.py" "$LOCK" <<'PY'
import hashlib,json,sys
(g,tsv,fa,builder,summarizer,out)=sys.argv[1:7]
def h(p):
 x=hashlib.sha256();
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(16<<20),b''):x.update(b)
 return x.hexdigest()
z=json.load(open(g))
assert z['gate_pass'] is True and z['structure_priority_shortlist_family_n']==32
assert z['shortlist_tsv_sha256']==h(tsv)
assert z['shortlist_fasta_sha256']==h(fa)
lock={
 'stage':'Stage05F3B3B-implementation-freeze','gate':'PASS_F3B3B_KNOWN_SYSTEM_AUDIT_RULE_FROZEN_BEFORE_SEQUENCE_SIMILARITY_RESULTS','gate_pass':True,
 'f3b3a_gate_sha256':h(g),'shortlist_tsv_sha256':h(tsv),'shortlist_fasta_sha256':h(fa),
 'reference_resolution':'explicit F3B3B_REFERENCE_FASTAS if set; otherwise deterministic scan of pre-F3B frozen ART project reference directories using path-role and protein-content guards',
 'swissprot_focus':'when prior D4R5B RVTALL_SWISSPROT.tblout exists, full SwissProt references are restricted to those pre-existing RVT-All positive-calibration IDs',
 'sequence_similarity_engine':'MMseqs2 easy-search','sequence_similarity_evalue':1e-3,'sequence_similarity_sensitivity':7.5,'sequence_similarity_max_seqs':100,
 'descriptive_exact_match_rule':'fident>=1.0, qcov>=0.99, tcov>=0.99','descriptive_high_identity_rule':'fident>=0.90, qcov>=0.80, tcov>=0.80',
 'reference_hit_used_to_remove_shortlist_family':False,'novelty_decision_authorized':False,
 'builder_sha256':h(builder),'summarizer_sha256':h(summarizer),
 'claim_boundary':'Post-freeze known-system audit only. Known matches are retained. No-hit is not novelty evidence.'
}
with open(out,'w') as f:json.dump(lock,f,indent=2);f.write('\n')
print(json.dumps(lock,indent=2))
PY
  sha256sum "$LOCK" "$S/scripts/build_known_reference_panel.py" "$S/scripts/summarize_known_audit.py" > "$S/frozen/F3B3B_FREEZE_SHA256.tsv"
  chmod a-w "$LOCK" "$S/scripts/build_known_reference_panel.py" "$S/scripts/summarize_known_audit.py"
else
  echo "[1/6] [SKIP] implementation lock exists"
fi

echo "[2/6] freeze byte-identical 32-sequence structure handoff"
cp -f "$QFA" "$S/outputs/F3B3B_STRUCTURE_HANDOFF_32.faa"
H1=$(sha256sum "$QFA" | awk '{print $1}')
H2=$(sha256sum "$S/outputs/F3B3B_STRUCTURE_HANDOFF_32.faa" | awk '{print $1}')
[[ "$H1" == "$H2" ]] || { echo '[FATAL] structure handoff FASTA hash mismatch' >&2; exit 2; }

echo "[3/6] build pre-existing known-system reference panel"
python "$S/scripts/build_known_reference_panel.py" "$P" "$S/references" | tee "$S/logs/reference_panel.log"
REFN=$(python - "$S/references/F3B3B_REFERENCE_PANEL_SUMMARY.json" <<'PY'
import json,sys
print(json.load(open(sys.argv[1]))['reference_unique_aa_n'])
PY
)

echo "reference_unique_aa_n=$REFN"

echo "[4/6] descriptive sequence-similarity audit"
HITS="$S/raw/F3B3B_SHORTLIST32_VS_KNOWN_REFERENCE.m8"
rm -f "$HITS"
if [[ "$REFN" -gt 0 ]]; then
  mmseqs easy-search \
    "$QFA" \
    "$S/references/F3B3B_KNOWN_SYSTEM_REFERENCE_PANEL.faa" \
    "$HITS" \
    "$S/tmp/mmseqs_known" \
    --threads "$THREADS" \
    -s 7.5 \
    -e 1e-3 \
    --max-seqs 100 \
    --mask 0 \
    --format-output 'query,target,fident,alnlen,qcov,tcov,evalue,bits' \
    > "$S/logs/mmseqs_known.log" 2>&1
  { echo -e 'query\ttarget\tfident\talnlen\tqcov\ttcov\tevalue\tbits'; cat "$HITS"; } > "$HITS.tmp"
  mv "$HITS.tmp" "$HITS"
else
  printf 'query\ttarget\tfident\talnlen\tqcov\ttcov\tevalue\tbits\n' > "$HITS"
fi

echo "[5/6] summarize; do not alter shortlist"
python "$S/scripts/summarize_known_audit.py" \
  "$QTSV" "$QFA" \
  "$S/references/F3B3B_REFERENCE_PANEL_SUMMARY.json" \
  "$S/references/F3B3B_KNOWN_SYSTEM_REFERENCE_MAP.tsv" \
  "$HITS" "$S/outputs" | tee "$S/logs/f3b3b_summary.log"

echo "[6/6] hashes and compact output"
sha256sum \
  "$LOCK" \
  "$S/outputs/F3B3B_STRUCTURE_HANDOFF_32.faa" \
  "$S/outputs/F3B3B_KNOWN_SYSTEM_AUDIT_32.tsv" \
  "$S/outputs/STAGE05F3B3B_GATE.json" \
  "$S/references/F3B3B_REFERENCE_INVENTORY.tsv" \
  "$S/references/F3B3B_REFERENCE_PANEL_SUMMARY.json" \
  > "$S/outputs/STAGE05F3B3B_SHA256.tsv"
cat "$S/outputs/STAGE05F3B3B_GATE.json"
echo
echo '=== KNOWN-SYSTEM AUDIT 32 ==='
column -t -s $'\t' "$S/outputs/F3B3B_KNOWN_SYSTEM_AUDIT_32.tsv" | head -33 || true
echo
echo 'PASS_STAGE05F3B3B_KNOWN_SYSTEM_AUDIT_AND_STRUCTURE_HANDOFF_COMPLETE'
