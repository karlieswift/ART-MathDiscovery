#!/usr/bin/env bash
set -euo pipefail

THREADS=${1:-32}
S=$(cd "$(dirname "$0")" && pwd)
P=/mnt/202001Main/BioData/hgq/ART_MathDiscovery_v1

S0="$P/ART_MATH_DISCOVERY_STAGE05F3B_DE_NOVO_ANNOTATION_INDEPENDENT_RT_DISCOVERY_v1_0_0_20261010"
S1="$P/ART_MATH_DISCOVERY_STAGE05F3B_DE_NOVO_ANNOTATION_INDEPENDENT_RT_DISCOVERY_v1_0_1_20261011"
OLD="$P/ART_MATH_DISCOVERY_STAGE05F3B2_FORMAL_CANDIDATE_FREEZE_AND_FAMILY_COMPRESSION_v1_0_0_20261011"
R1="$P/ART_MATH_DISCOVERY_STAGE05F3B2_R1_ORFDB_PATH_HOTFIX_v1_0_1_20261011"

ORFDB="$S0/work/F3B1_orfs_aa"
ORFIDX="$ORFDB.index"
ORFDBH="${ORFDB}_h"
ORFHDX="${ORFDB}_h.index"
HMM_INPUT="$S1/work/F3B1_ALL_ORFS_NUMERIC.faa"

OLD_FORMAL="$OLD/frozen/F3B2_FORMAL_CANDIDATE_IDS.txt"
OLD_LOCK="$OLD/frozen/F3B2_IMPLEMENTATION_LOCK.json"
F3B1_FORMAL="$S1/outputs/STAGE05F3B1_FORMAL_IDS.txt"
F3B1_PAIRS="$S1/outputs/STAGE05F3B1_RVT_PAIR_HITS.tsv.gz"
F3B1_GATE="$S1/outputs/STAGE05F3B1_GATE.json"
F3B1_AUDIT="$S1/outputs/STAGE05F3B1_SHARDING_POSTRUN_AUDIT.json"
F3B1_RUNTIME_LOCK="$S1/frozen/F3B0_V101_IMPLEMENTATION_LOCK.json"

FORMAL="$S/frozen/F3B2_FORMAL_CANDIDATE_IDS.txt"
LOCK="$S/frozen/F3B2_R2_IMPLEMENTATION_LOCK.json"
ANALYSIS="$S/scripts/stage05f3b2.py"
EXTRACTOR="$S/scripts/extract_formal_from_serialized_fasta.py"

mkdir -p "$S"/{frozen,raw,work,tmp,outputs,logs,docs}

for x in python mmseqs sha256sum; do
  command -v "$x" >/dev/null || { echo "[FATAL] missing tool: $x" >&2; exit 2; }
done

for x in \
  "$OLD_FORMAL" "$OLD_LOCK" "$F3B1_FORMAL" "$F3B1_PAIRS" \
  "$F3B1_GATE" "$F3B1_AUDIT" "$F3B1_RUNTIME_LOCK" \
  "$ANALYSIS" "$EXTRACTOR" "$HMM_INPUT" \
  "$ORFDB" "$ORFIDX" "$ORFDB.dbtype" "$ORFDBH" "$ORFHDX" "${ORFDB}_h.dbtype"
do
  [[ -e "$x" ]] || { echo "[FATAL] required input missing: $x" >&2; exit 2; }
done

echo "[1/7] freeze R2 query-ID/key-mapping hotfix before family results"

python - \
  "$OLD_FORMAL" "$OLD_LOCK" "$F3B1_FORMAL" "$F3B1_GATE" "$F3B1_AUDIT" \
  "$F3B1_RUNTIME_LOCK" "$ANALYSIS" "$EXTRACTOR" "$HMM_INPUT" \
  "$FORMAL" "$LOCK" "$ORFDB" "$R1" <<'PY'
import hashlib,json,os,shutil,sys
from pathlib import Path

(old_formal,old_lock,f3b1_formal,f3b1_gate,f3b1_audit,runtime_lock,
 analysis,extractor,hmm_input,formal_out,lock_out,orfdb,r1)=sys.argv[1:14]

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''):
            h.update(b)
    return h.hexdigest()

def ids(p): return [x.strip() for x in open(p) if x.strip()]

old_ids=ids(old_formal); f3_ids=ids(f3b1_formal)
assert len(old_ids)==10464 and len(set(old_ids))==10464
assert len(f3_ids)==10464 and set(old_ids)==set(f3_ids)
og=json.load(open(old_lock)); g=json.load(open(f3b1_gate)); a=json.load(open(f3b1_audit)); r=json.load(open(runtime_lock))
assert og['gate_pass'] is True and og['formal_candidate_n']==10464
assert g['gate_pass'] is True and g['formal_candidate_n']==10464
assert a['gate_pass'] is True and a['sum_shard_record_n']==65315646
assert r['orf_instance_n']==65315646
expected={'family_min_seq_id':0.30,'family_coverage':0.70,'family_cov_mode':0,'family_cluster_mode':0,'family_sensitivity':7.5}
for k,v in expected.items(): assert og[k]==v,(k,og[k],v)
Path(formal_out).parent.mkdir(parents=True,exist_ok=True)
shutil.copyfile(old_formal,formal_out)
assert sha(formal_out)==sha(old_formal)

lock={
 'stage':'Stage05F3B2-R2',
 'gate':'PASS_F3B2_R2_SERIALIZED_QUERYID_KEYMAP_HOTFIX_FROZEN_BEFORE_FAMILY_RESULTS',
 'gate_pass':True,
 'hotfix_type':'IMPLEMENTATION_IDENTIFIER_MAPPING_ONLY',
 'r1_failure':'createsubdb treated hmmscan serialized query IDs as MMseqs DB keys; 10232/10464 were recovered, proving the identifiers are not universally interchangeable',
 'canonical_selection_source':'exact F3B1 hmmscan input FASTA query identifiers',
 'canonical_hmm_input_fasta':hmm_input,
 'canonical_hmm_input_fasta_sha256':sha(hmm_input),
 'mapping_method':'pair each serialized FASTA record ordinal with the original ORFDB.index row at the same ordinal; validate all 10,464 mapped peptides against ORFDB bytes before use',
 'corrected_orfdb_path':orfdb,
 'scientific_protocol_changed':False,
 'formal_candidate_membership_changed':False,
 'family_parameters_changed':False,
 'known_rt_annotation_used':False,
 'novelty_screen_performed':False,
 'formal_candidate_n':10464,
 'formal_candidate_ids_sha256':sha(formal_out),
 'source_v100_candidate_ids_sha256':sha(old_formal),
 'source_v100_implementation_lock_sha256':sha(old_lock),
 'f3b1_gate_sha256':sha(f3b1_gate),
 'f3b1_sharding_audit_sha256':sha(f3b1_audit),
 'f3b1_runtime_lock_sha256':sha(runtime_lock),
 'f3b1_orf_instance_n':65315646,
 'f3b1_orf_fasta_sha256':r['orf_fasta_sha256'],
 'family_min_seq_id':0.30,'family_coverage':0.70,'family_cov_mode':0,'family_cluster_mode':0,'family_sensitivity':7.5,
 'analysis_script_sha256':sha(analysis),
 'extractor_script_sha256':sha(extractor),
 'claim_boundary':'R2 corrects identifier mapping only. Frozen F3B1 formal membership and F3B2 family-compression parameters are unchanged.'
}
with open(lock_out,'w') as f: json.dump(lock,f,indent=2); f.write('\n')
print(json.dumps(lock,indent=2))
PY

sha256sum "$FORMAL" "$LOCK" "$ANALYSIS" "$EXTRACTOR" > "$S/frozen/F3B2_R2_FREEZE_SHA256.tsv"
chmod a-w "$FORMAL" "$LOCK" "$ANALYSIS" "$EXTRACTOR"

echo "[2/7] extract all 10,464 formal peptides by exact hmmscan query IDs"
rm -f "$S/raw/F3B2_FORMAL_NUMERIC.faa" "$S/raw/F3B2_FORMAL_WITH_SOURCE_HEADERS.faa" "$S/raw/F3B2_QUERY_TO_ORFDB_KEY.tsv"

python "$EXTRACTOR" \
  "$FORMAL" "$HMM_INPUT" "$ORFDB" "$ORFIDX" "$ORFDBH" "$ORFHDX" \
  "$S/raw/F3B2_FORMAL_NUMERIC.faa" \
  "$S/raw/F3B2_FORMAL_WITH_SOURCE_HEADERS.faa" \
  "$S/raw/F3B2_QUERY_TO_ORFDB_KEY.tsv" \
  | tee "$S/logs/f3b2_r2_extract.log"

python - "$FORMAL" "$S/raw/F3B2_FORMAL_NUMERIC.faa" "$S/raw/F3B2_FORMAL_WITH_SOURCE_HEADERS.faa" <<'PY'
import sys

def fasta(path):
    h=None;s=[]
    with open(path,encoding='ascii') as f:
        for line in f:
            line=line.rstrip('\r\n')
            if line.startswith('>'):
                if h is not None: yield h,''.join(s).upper()
                h=line[1:];s=[]
            elif line:s.append(line.strip())
        if h is not None: yield h,''.join(s).upper()
formal=[x.strip() for x in open(sys.argv[1]) if x.strip()]
num=list(fasta(sys.argv[2])); hdr=list(fasta(sys.argv[3]))
assert len(formal)==10464
assert len(num)==10464,len(num)
assert len(hdr)==10464,len(hdr)
assert [h.split()[0] for h,_ in num]==formal
for (nh,ns),(hh,hs) in zip(num,hdr): assert ns==hs,nh
print('formal_candidate_n =',len(formal))
print('numeric_fasta_n =',len(num))
print('source_header_fasta_n =',len(hdr))
print('PASS_F3B2_R2_FORMAL_ORF_EXTRACTION')
PY

echo "[3/7] exact-AA compression and HMM descriptive features"
python "$ANALYSIS" prepare \
  "$FORMAL" \
  "$S/raw/F3B2_FORMAL_NUMERIC.faa" \
  "$S/raw/F3B2_FORMAL_WITH_SOURCE_HEADERS.faa" \
  "$F3B1_PAIRS" \
  "$S/outputs" \
  | tee "$S/logs/f3b2_prepare.log"

echo "[4/7] operational sequence-family clustering"
rm -rf "$S/work/F3B2_exactDB"* "$S/work/F3B2_familyDB"* "$S/tmp/F3B2_family"
mmseqs createdb "$S/outputs/F3B2_EXACT_AA_REPRESENTATIVES.faa" "$S/work/F3B2_exactDB" --shuffle 0 > "$S/logs/f3b2_createdb_exact.log" 2>&1
mmseqs cluster "$S/work/F3B2_exactDB" "$S/work/F3B2_familyDB" "$S/tmp/F3B2_family" \
  --min-seq-id 0.30 -c 0.70 --cov-mode 0 --cluster-mode 0 -s 7.5 --threads "$THREADS" \
  > "$S/logs/f3b2_family_cluster.log" 2>&1
mmseqs createtsv "$S/work/F3B2_exactDB" "$S/work/F3B2_exactDB" "$S/work/F3B2_familyDB" "$S/raw/F3B2_FAMILY_CLUSTER.tsv" \
  > "$S/logs/f3b2_family_tsv.log" 2>&1

echo "[5/7] finalize frozen F3B2 family summary"
python "$ANALYSIS" finalize \
  "$S/outputs/F3B2_CANDIDATE_BASE.tsv" \
  "$S/outputs/F3B2_EXACT_AA_MAP.tsv" \
  "$S/raw/F3B2_FAMILY_CLUSTER.tsv" \
  "$LOCK" "$S/outputs" \
  | tee "$S/logs/f3b2_finalize.log"

echo "[6/7] integrity hashes"
sha256sum \
  "$FORMAL" "$LOCK" \
  "$S/raw/F3B2_QUERY_TO_ORFDB_KEY.tsv" \
  "$S/raw/F3B2_FORMAL_NUMERIC.faa" \
  "$S/raw/F3B2_FORMAL_WITH_SOURCE_HEADERS.faa" \
  "$S/outputs/F3B2_EXACT_AA_MAP.tsv" \
  "$S/outputs/F3B2_CANDIDATE_BASE.tsv" \
  "$S/raw/F3B2_FAMILY_CLUSTER.tsv" \
  "$S/outputs/F3B2_CANDIDATE_FAMILY_MAP.tsv" \
  "$S/outputs/F3B2_FAMILY_SUMMARY.tsv" \
  "$S/outputs/STAGE05F3B2_GATE.json" \
  > "$S/outputs/STAGE05F3B2_R2_SHA256.tsv"

echo "[7/7] compact result"
cat "$S/outputs/STAGE05F3B2_GATE.json"
echo
echo "=== TOP 20 FAMILIES ==="
column -t -s $'\t' "$S/outputs/F3B2_FAMILY_TOP100.tsv" | head -21 || true

echo
echo "PASS_STAGE05F3B2_R2_SERIALIZED_QUERYID_KEYMAP_AND_FAMILY_COMPRESSION_COMPLETE"
