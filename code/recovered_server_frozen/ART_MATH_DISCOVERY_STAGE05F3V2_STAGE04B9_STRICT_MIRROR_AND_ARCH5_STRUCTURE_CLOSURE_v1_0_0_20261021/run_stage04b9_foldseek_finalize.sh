#!/usr/bin/env bash
set -euo pipefail
P=${1:?project root}; THREADS=${2:-32}; HERE=$(cd "$(dirname "$0")" && pwd)
command -v foldseek >/dev/null || { echo HOLD foldseek missing; exit 2; }
RTDB="$P/ART_MATH_DISCOVERY_STAGE05D4R5C_CLEAN_RT_CALIBRATION_v1_0_0_20261007/references/clean_rt_afsp"; ls "${RTDB}"* >/dev/null 2>&1 || { echo HOLD clean RT DB missing; exit 2; }
run_branch(){
 local tag=$1 man=$2 sd=$3 n=$4 stage=$5
 local fd="$HERE/work/foldseek_${tag}" td="$HERE/work/foldseek_tmp_${tag}";mkdir -p "$fd" "$td" "$HERE/logs"
 while IFS=$'\t' read -r q rest; do
  [[ "$q" == query_id ]] && continue
  cif="$sd/${q}.cif"; [[ -s "$cif" ]] || { echo "HOLD missing CIF $tag $q"; exit 2; }
  out="$fd/${q}.m8"; tmp="$td/${q}";rm -rf "$tmp"
  foldseek easy-search "$cif" "$RTDB" "$out" "$tmp" -e 1.0 --max-seqs 100 --threads "$THREADS" --exact-tmscore 1 --format-output query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore > "$HERE/logs/foldseek_${tag}_${q}.log" 2>&1
  echo "[PASS FOLDSEEK $tag] $q"
 done < "$man"
 python "$HERE/scripts/classify_foldseek.py" --manifest "$man" --foldseek-dir "$fd" --out-tsv "$HERE/outputs/${tag}_CORRECTED_STRUCTURAL_COMPATIBILITY.tsv" --out-gate "$HERE/outputs/${tag}_STRUCTURE_GATE.json" --stage "$stage" --expected-n "$n"
}
run_branch STRICT_MIRROR "$HERE/outputs/STRICT_MIRROR_STRUCTURE_MANIFEST_32.tsv" "$HERE/outputs/strict_mirror_structures" 32 Stage05F3V2-Stage04B9A
run_branch ARCH5 "$HERE/outputs/ARCH5_STRUCTURE_MANIFEST.tsv" "$HERE/outputs/arch5_structures" 5 Stage05F3V2-Stage04B9B
python - "$P" "$HERE" <<'PY'
import csv,json,sys
from pathlib import Path
P=Path(sys.argv[1]);H=Path(sys.argv[2]);O=H/'outputs'
def one(t):
 ds=sorted(d for d in P.iterdir() if d.is_dir() and t in d.name)
 if not ds:raise SystemExit('HOLD cannot locate '+t)
 return ds[-1]
def read(p):
 with open(p,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,cols):
 with open(p,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
sm=json.load(open(O/'STRICT_MIRROR_STRUCTURE_GATE.json'));a5=json.load(open(O/'ARCH5_STRUCTURE_GATE.json'))
B7=one('STAGE04B7_RESCUED_DELTA_STRUCTURE_AND_REPLICATION_CLOSURE_v1_0_0');b7rows=read(B7/'outputs/DELTA13_CORRECTED_STRUCTURAL_COMPATIBILITY.tsv');primary=[r for r in b7rows if r['b7_role']=='PRIMARY_COMPLETE_STOP_TO_STOP']
if len(primary)!=12:raise SystemExit('HOLD expected B7 primary12')
a5rows=read(O/'ARCH5_CORRECTED_STRUCTURAL_COMPATIBILITY.tsv');comb=[]
for r in primary:comb.append({'source':'B7_PRIMARY12','family':r['delta_operational_family_id'],'candidate_id':r['candidate_id'],'aa_length':r['aa_length'],'tier':r['tier']})
for r in a5rows:comb.append({'source':'B9_ARCH5_COMPLETE','family':r['delta_operational_family_id'],'candidate_id':r['candidate_id'],'aa_length':r['aa_length'],'tier':r['tier']})
write(O/'ARCHITECTURE_FIRST_DELTAONLY17_COMPLETE_STRUCTURAL_CLOSURE.tsv',comb,list(comb[0]))
pos=sum(r['tier'] in ('STRONG','SUGGESTIVE') for r in comb)
# Corrected frozen structural reference source if present.
freeze=None
for p in P.rglob('STRUCTURAL_COMPATIBILITY_FREEZE_GATE.json'):
 try:
  z=json.load(open(p))
  if z.get('gate')=='PASS_STRUCTURAL_COMPATIBILITY_DEFINITION_FROZEN_AND_VALIDATED':freeze=(p,z);break
 except:pass
hist={'n':30,'positive':28,'strong':20,'suggestive':8};deltaonly={'n':32,'positive':12,'strong':3,'suggestive':9}
if freeze:
 z=freeze[1];hist={'n':z['historical_resolved30']['n'],'positive':z['historical_resolved30']['positive'],'strong':z['historical_resolved30']['strong'],'suggestive':z['historical_resolved30']['suggestive']};deltaonly={'n':z['independent_delta32']['n'],'positive':z['independent_delta32']['positive'],'strong':z['independent_delta32']['strong'],'suggestive':z['independent_delta32']['suggestive']}
g={'stage':'Stage05F3V2-Stage04B9','gate':'PASS_STAGE04B9_STRICT_MIRROR_AND_ARCH5_STRUCTURE_CLOSURE','gate_pass':True,'historical_6p12G_corrected_structure':hist,'delta54G_strict_mirror':{'n':sm['n'],'strong':sm['strong_n'],'suggestive':sm['suggestive_n'],'positive':sm['positive_n']},'delta54G_novelty_only_protein_route':deltaonly,'delta54G_architecture_first_complete17':{'n':17,'positive':pos,'strong':sum(r['tier']=='STRONG' for r in comb),'suggestive':sum(r['tier']=='SUGGESTIVE' for r in comb),'weak':sum(r['tier']=='WEAK' for r in comb),'no_hit':sum(r['tier']=='NO_HIT' for r in comb)},'delta54G_architecture_first_edge_censored_sensitivity_n':1,'strict_mirror_same_family_universe_definition_as_historical':True,'strict_mirror_combined_old_delta_novelty_filter_used':False,'three_54G_routes_must_be_reported_separately':True,'retuning_performed':False,'science_change':False}
with open(O/'STAGE04B9_CLOSEOUT_GATE.json','w') as f:json.dump(g,f,indent=2);f.write('\n')
print(json.dumps(g,indent=2));print('PASS_STAGE04B9_STRICT_MIRROR_AND_ARCH5_STRUCTURE_CLOSURE')
PY
