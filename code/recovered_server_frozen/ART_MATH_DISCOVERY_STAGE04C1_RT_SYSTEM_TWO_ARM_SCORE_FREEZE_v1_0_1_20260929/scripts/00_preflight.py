#!/usr/bin/env python3
import argparse,csv,hashlib,json
from pathlib import Path

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); ap.add_argument('--c0-root',required=True); ap.add_argument('--seqfreeze-root',required=True); a=ap.parse_args()
 root=Path(a.root); c0=Path(a.c0_root); sf=Path(a.seqfreeze_root)
 pol=json.loads((root/'config/C1_POLICY.json').read_text())
 reg=list(csv.DictReader((c0/'registry/RT_SYSTEM_CONTROL_REGISTRY.tsv').open(),delimiter='\t')); regmap={r['control_id']:r for r in reg}
 gate_path=sf/'outputs/STAGE04C0_R3_FINAL_GATE.json'; man_path=sf/'outputs/STAGE04C0_R3_FIXED_SEQUENCE_MANIFEST.tsv'
 if not gate_path.is_file(): raise SystemExit(f'[FAIL] required R3 gate missing: {gate_path}')
 if not man_path.is_file(): raise SystemExit(f'[FAIL] required R3 manifest missing: {man_path}')
 gate=json.loads(gate_path.read_text())
 if not gate.get('gate_pass') or gate.get('gate')!='PASS_STAGE04C0_R3_SEQUENCE_LENGTH_CORRECTION_FIXED':
  raise SystemExit('[FAIL] Stage04C0-R3 sequence freeze gate not PASS')
 man=list(csv.DictReader(man_path.open(),delimiter='\t')); mmap={r['control_id']:r for r in man}
 rows=[]
 for cid in pol['eligible_control_ids']:
  if cid not in regmap: raise SystemExit(f'[FAIL] control absent registry: {cid}')
  if regmap[cid]['c1_scan_eligible_now']!='1': raise SystemExit(f'[FAIL] control not C1 eligible: {cid}')
  m=mmap.get(cid)
  if not m or m['status']!='PRESENT_VALIDATED': raise SystemExit(f'[FAIL] no validated R3 sequence: {cid}')
  acc=m['accession']; fa=sf/'inputs/sequences'/f'{cid}__{acc}.fasta'
  if not fa.is_file(): raise SystemExit(f'[FAIL] fasta missing {fa}')
  if sha(fa)!=m['fasta_sha256']: raise SystemExit(f'[FAIL] fasta hash mismatch {cid}')
  rows.append({'control_id':cid,'registry_accession':regmap[cid]['sequence_accession'],'accession':acc,'fasta':str(fa),'length_bp':int(m['length_bp']),'fasta_sha256':m['fasta_sha256']})
 out={'gate':'PASS_STAGE04C1_TWO_ARM_PREFLIGHT','gate_pass':True,'truth_used':False,'eligible_controls':rows,
      'sequence_freeze_gate':gate['gate'],'sequence_freeze_manifest_sha256':sha(man_path),
      'retuning_authorized':False}
 (root/'outputs/C1_PREFLIGHT_GATE.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
if __name__=='__main__': main()
