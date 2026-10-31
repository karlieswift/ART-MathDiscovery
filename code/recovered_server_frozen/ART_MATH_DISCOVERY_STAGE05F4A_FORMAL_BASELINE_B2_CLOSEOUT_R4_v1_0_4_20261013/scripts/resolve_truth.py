#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path

def sha256(p:Path):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def read_tsv(p):
    with open(p,newline='') as f: return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(p,rows,fields):
    p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore'); w.writeheader(); w.writerows(rows)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--project-root',required=True); a=ap.parse_args()
    here=Path(__file__).resolve().parents[1]; P=Path(a.project_root).resolve()
    imp=here/'frozen/IMPORTED_R2_B1_ARCHIVE.json'
    if not imp.is_file(): raise SystemExit('import frozen R2 B1 archive first')
    iobj=json.loads(imp.read_text())
    if not iobj.get('gate_pass'): raise SystemExit('R2 B1 import gate not PASS')
    cases=read_tsv(here/'config/CASES.tsv')

    repeat_file=P/'ART_MATH_DISCOVERY_STAGE05C2B_TEN_LOCUS_VARIATION_RNA_v1_0_0_20261006/inputs/evidence/CORE_CANDIDATE_QUALIFYING_REPEAT_DESCRIPTIONS.tsv'
    ctrl_file=P/'ART_MATH_DISCOVERY_STAGE04C1B0_R2_DRT2_SEQUENCE_IDENTITY_RESOLUTION_v1_0_2_20261002/outputs/C1B0_CONTROL_TRUTH_RESOLVED.tsv'
    for p in [repeat_file,ctrl_file]:
        if not p.is_file(): raise SystemExit(f'missing truth source: {p}')
    rep=read_tsv(repeat_file)
    artmap={'LY01':'OR836606.1','S6':'LC680885.1','PALS_2':'MN091626.1','UFV_DC4':'MZ779063.1'}
    ctrmap={
      'SEN2':('SEN2_14028S','RT_rrtT_STM14_4641'),
      'EC78_ECO7':('EC78_ECO7','RT_Ec78'),
      'DRT2':('DRT2_KPN','RT_DRT2'),
    }
    ctrl=[r for r in read_tsv(ctrl_file) if r.get('control_id') and r.get('control_id')!='control_id']
    rep_sha=sha256(repeat_file); ctrl_sha=sha256(ctrl_file)

    for c in cases:
        cid=c['case_id']
        c['truth_source_class']='PRELOCKED'
        c['truth_source_path']='config/CASES.tsv'
        c['truth_source_sha256']=sha256(here/'config/CASES.tsv')
        c['truth_resolution_rule']='prelocked interval carried unchanged'
        c['truth_independence_tier']='A_PRELOCKED_EXTERNAL_OR_RETROSPECTIVE'
        if cid in artmap:
            sid=artmap[cid]
            rows=[r for r in rep if r.get('sequence_id')==sid and r.get('assay')=='EXACT_CORE_PRIMARY_DESCRIPTION']
            if not rows: raise SystemExit(f'no EXACT_CORE_PRIMARY_DESCRIPTION rows for {cid}/{sid}')
            starts=[int(r['start0']) for r in rows]; ends=[int(r['end0']) for r in rows]
            c['truth_start0']=str(min(starts)); c['truth_end0']=str(max(ends))
            c['truth_mode']='PROJECT_POSTFREEZE_EXACT_REPEAT_ENVELOPE'
            c['truth_source_class']='PROJECT_RECONSTRUCTED_EXACT_REPEAT'
            c['truth_source_path']=str(repeat_file)
            c['truth_source_sha256']=rep_sha
            c['truth_resolution_rule']='min(start0)..max(end0) over all EXACT_CORE_PRIMARY_DESCRIPTION rows for the frozen sequence_id; no current baseline score/rank consulted'
            c['truth_independence_tier']='B_PROJECT_RECONSTRUCTED_PUBLIC_SOURCE'
            c['truth_exact_row_n']=str(len(rows))
        elif cid in ctrmap:
            co,comp=ctrmap[cid]
            rows=[r for r in ctrl if r.get('control_id')==co and r.get('component')==comp and r.get('truth_status')=='RESOLVED']
            if len(rows)!=1: raise SystemExit(f'expected one resolved control truth row for {cid}, found {len(rows)}')
            r=rows[0]
            c['truth_start0']=r['start0']; c['truth_end0']=r['end0']; c['truth_mode']='EXACT_RT_ANCHOR'
            c['truth_source_class']='HISTORICAL_POSTFREEZE_CONTROL_TRUTH'
            c['truth_source_path']=str(ctrl_file); c['truth_source_sha256']=ctrl_sha
            c['truth_resolution_rule']=f'exact resolved component {comp}; interval copied unchanged from C1B0 control truth resolution'
            c['truth_independence_tier']='A_HISTORICAL_CONTROL_EXACT_RT'
        elif cid=='PB50':
            c['truth_independence_tier']='C_RT_PROXY_PRELOCKED'
            c['truth_source_class']='PRELOCKED_RT_PROXY'
            c['truth_resolution_rule']='prelocked RT proxy interval carried unchanged; no exact array truth asserted'
    unresolved=[c['case_id'] for c in cases if c.get('primary_eval')=='1' and (c.get('truth_start0') in {'','.','NA'} or c.get('truth_end0') in {'','.','NA'})]
    if unresolved: raise SystemExit('unresolved truth remains: '+','.join(unresolved))
    fields=['case_id','reference_accession','sequence_accession','role','group','primary_eval','truth_mode','truth_start0','truth_end0','truth_independence_tier','truth_source_class','truth_source_path','truth_source_sha256','truth_resolution_rule','truth_exact_row_n','notes']
    out=here/'frozen/FULLY_RESOLVED_BASELINE_TRUTH_SEAL.tsv'; write_tsv(out,cases,fields)
    obj={
      'stage':'Stage05F4A-B2-R4-postfreeze-truth-resolution','version':'v1.0.4',
      'gate':'PASS_STAGE05F4A_B2_R4_FULLY_RESOLVED_TRUTH_SEAL', 'gate_pass':True,
      'truth_resolution_after_b1_score_freeze':True,'truth_used_before_score_freeze':False,'retuning_authorized':False,
      'primary_positive_n':sum(c['group']=='ART_POSITIVE' and c['role']=='PRIMARY_EXTERNAL' for c in cases),
      'primary_exact_array_n':sum(c['group']=='ART_POSITIVE' and c['role']=='PRIMARY_EXTERNAL' and c['truth_mode'] in {'STRICT_REPEAT_CORE','PROJECT_POSTFREEZE_EXACT_REPEAT_ENVELOPE'} for c in cases),
      'primary_rt_proxy_n':sum(c['group']=='ART_POSITIVE' and c['role']=='PRIMARY_EXTERNAL' and c['truth_mode']=='RT_PROXY' for c in cases),
      'mechanistic_control_n':sum(c['group']=='NONARRAY_RT_CONTROL' for c in cases),
      'truth_seal_path':str(out),'truth_seal_sha256':sha256(out),
      'note':'The pre-registered 8-case primary panel is preserved. A 7-case exact-array subset is labeled descriptive sensitivity; PB50 remains the prelocked RT-proxy challenge case.'
    }
    (here/'frozen/BASELINE_TRUTH_RESOLUTION_GATE.json').write_text(json.dumps(obj,indent=2)+"\n")
    print(json.dumps(obj,indent=2))
if __name__=='__main__': main()
