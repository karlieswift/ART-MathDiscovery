#!/usr/bin/env python3
from pathlib import Path
import argparse, json, sys
from common import read_tsv, write_tsv, dump_json, strict_overlap, overlap_bp, gap_bp, sha256

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--r4a-out',required=True); ap.add_argument('--r2d-root',required=True); ap.add_argument('--out',required=True)
    a=ap.parse_args(); r4=Path(a.r4a_out); r2d=Path(a.r2d_root); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    gate=json.loads((r4/'STAGE04C2_R4A_TRUTH_GATE.json').read_text())
    if not gate.get('gate_pass'): raise SystemExit('R4A truth gate is not PASS; exact-array evaluation forbidden')
    units=read_tsv(r4/'UG27_EXACT_NCRNA_UNITS.tsv'); spans=read_tsv(r4/'UG27_EXACT_NCRNA_ARRAY_SPANS.tsv')
    allr=read_tsv(r2d/'outputs/CLEAN_TWO_ARM_ALL_REGIONS.tsv'); dual=read_tsv(r2d/'outputs/CLEAN_TWO_ARM_DUAL_ARM_VIEW.tsv')
    byall={}; bydual={}
    for r in allr: byall.setdefault((r['assembly_id'],r['scoring_id']),[]).append(r)
    for r in dual: bydual.setdefault((r['assembly_id'],r['scoring_id']),[]).append(r)
    # map system to assembly/scoring from R4 spans: source contig -> find matching scoring row through R2D region tables impossible if no region; infer fixed EXP order
    sysmap={'Alistipes-putredinis':('UG27_EXP1','UG27_EXP1__C001'),'Bacteroides-coprophilus':('UG27_EXP2','UG27_EXP2__C001'),'Bacteroides-sp-9-1-42FAA':('UG27_EXP3','UG27_EXP3__C001'),'Bacteroides-sp-2-1-16':('UG27_EXP4','UG27_EXP4__C001'),'Blautia-hansenii':('UG27_EXP5','UG27_EXP5__C001')}
    unit_rows=[]; summary=[]
    for sp in spans:
        s=sp['system']; assembly,scoring=sysmap[s]; uu=[u for u in units if u['system']==s]; rr=bydual.get((assembly,scoring),[]); ar=byall.get((assembly,scoring),[])
        hit_units=0; hit_dual_ids=set(); hit_acnv_only_units=0
        for u in uu:
            u0=int(u['start0']); u1=int(u['end0'])
            dh=[]
            for r in rr:
                r0=int(r['start0']); r1=int(r['end0'])
                if strict_overlap(u0,u1,r0,r1): dh.append((r,overlap_bp(u0,u1,r0,r1)))
            ah=[]
            for r in ar:
                r0=int(r['start0']); r1=int(r['end0'])
                if strict_overlap(u0,u1,r0,r1): ah.append((r,overlap_bp(u0,u1,r0,r1)))
            if dh:
                hit_units+=1
                for r,_ in dh: hit_dual_ids.add(r['region_id'])
            elif any('ACNV_FAMILY' in r['support_families'] for r,_ in ah):
                hit_acnv_only_units+=1
            best=sorted(dh,key=lambda z:int(z[0]['region_rank']))[0] if dh else None
            unit_rows.append({'system':s,'unit_id':u['unit_id'],'unit_start0':u0,'unit_end0':u1,'dual_arm_direct':bool(dh),'dual_arm_overlap_n':len(dh),'best_dual_region_rank':best[0]['region_rank'] if best else '','best_dual_region_id':best[0]['region_id'] if best else '','best_dual_overlap_bp':best[1] if best else 0,'acnv_only_direct_without_dual':(not dh and any('ACNV_FAMILY' in r['support_families'] for r,_ in ah))})
        a0=int(sp['array_start0']); a1=int(sp['array_end0'])
        span_hits=[r for r in rr if strict_overlap(a0,a1,int(r['start0']),int(r['end0']))]
        nearest=min(rr,key=lambda r:gap_bp(a0,a1,int(r['start0']),int(r['end0']))) if rr else None
        ng=gap_bp(a0,a1,int(nearest['start0']),int(nearest['end0'])) if nearest else ''
        summary.append({'system':s,'assembly_id':assembly,'scoring_id':scoring,'exact_unit_n':len(uu),'exact_unit_dual_arm_overlap_n':hit_units,'exact_unit_dual_arm_overlap_fraction':hit_units/len(uu) if uu else 0,'primary_exact_unit_dual_arm_positive':hit_units>0,'exact_array_span_start0':a0,'exact_array_span_end0':a1,'array_span_dual_arm_overlap':bool(span_hits),'array_span_dual_arm_overlap_n':len(span_hits),'direct_dual_region_ids':';'.join(sorted(hit_dual_ids)),'exact_units_acnv_only_direct_n':hit_acnv_only_units,'nearest_dual_arm_gap_bp':ng,'nearest_dual_region_rank':nearest['region_rank'] if nearest else '','nearest_dual_region_id':nearest['region_id'] if nearest else ''})
    write_tsv(out/'UG27_EXACT_UNIT_OVERLAPS.tsv',unit_rows,['system','unit_id','unit_start0','unit_end0','dual_arm_direct','dual_arm_overlap_n','best_dual_region_rank','best_dual_region_id','best_dual_overlap_bp','acnv_only_direct_without_dual'])
    write_tsv(out/'UG27_EXACT_ARRAY_EVALUATION_SUMMARY.tsv',summary,['system','assembly_id','scoring_id','exact_unit_n','exact_unit_dual_arm_overlap_n','exact_unit_dual_arm_overlap_fraction','primary_exact_unit_dual_arm_positive','exact_array_span_start0','exact_array_span_end0','array_span_dual_arm_overlap','array_span_dual_arm_overlap_n','direct_dual_region_ids','exact_units_acnv_only_direct_n','nearest_dual_arm_gap_bp','nearest_dual_region_rank','nearest_dual_region_id'])
    pos=sum(1 for x in summary if x['primary_exact_unit_dual_arm_positive'])
    final={'stage':'Stage04C2-R4B','gate':'PASS_STAGE04C2_R4B_EXACT_NCRNA_ARRAY_POSTHOC_EVALUATION_COMPLETE_NO_RETUNE','gate_pass':True,'n_systems':len(summary),'primary_exact_unit_dual_arm_positive_n':pos,'primary_exact_unit_dual_arm_positive_fraction':pos/len(summary) if summary else None,'primary_endpoint':'strict overlap of frozen clean dual-arm region with at least one exact author-defined ncRNA unit','nearby_threshold_defined':False,'truth_recomputed_from_detector':False,'score_recomputed':False,'region_board_recomputed':False,'retuning_authorized':False,'science_change':False,'summary_sha256':sha256(out/'UG27_EXACT_ARRAY_EVALUATION_SUMMARY.tsv'),'unit_overlaps_sha256':sha256(out/'UG27_EXACT_UNIT_OVERLAPS.tsv')}
    dump_json(out/'STAGE04C2_R4B_GATE.json',final); print(json.dumps(final,indent=2))
if __name__=='__main__': main()
