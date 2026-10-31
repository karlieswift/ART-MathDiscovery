#!/usr/bin/env python3
import argparse,csv,json
from pathlib import Path
import pandas as pd
from common import *
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); a=ap.parse_args(); root=Path(a.root).resolve(); pol=read_json(root/'config/ENGINE_POLICY.json')
    eng=read_json(root/'outputs/STAGE03A_ENGINE_FREEZE_GATE.json'); fg=read_json(root/'outputs/STAGE03A_FASTA_FREEZE_GATE.json')
    with open(root/'config/VALIDATION_COHORT.tsv') as f: cohort=list(csv.DictReader(f,delimiter='\t'))
    rows=[]; ok=bool(eng.get('gate_pass') and fg.get('gate_pass'))
    for c in cohort:
        acc=c['accession']; out=root/'outputs'/acc; pg=read_json(out/'STAGE03A_PRIMARY_GATE.json'); cg=read_json(out/'STAGE03A_CONFIRMATION_GATE.json'); ok &= bool(pg.get('gate_pass') and cg.get('gate_pass'))
        reg=pd.read_csv(out/'STAGE03A_UNION_REGIONS.tsv',sep='\t'); board=pd.read_csv(out/'STAGE03A_CONFIRMATION_BOARD.tsv',sep='\t') if (out/'STAGE03A_CONFIRMATION_BOARD.tsv').stat().st_size>1 else pd.DataFrame()
        top_reg=reg.sort_values('region_rank').iloc[0] if len(reg) else None; top_conf=board.sort_values('confirmation_rank').iloc[0] if len(board) else None
        rows.append({'accession':acc,'name':c['name'],'n_genomic_windows':pg['n_genomic_windows'],'n_union_regions':pg['n_union_regions'],'primary_union_regions_sha256':pg['union_regions_sha256'],'confirmation_board_sha256':cg.get('confirmation_board_sha256',''),'top_primary_region_id':str(top_reg.region_id) if top_reg is not None else '','top_primary_start0':int(top_reg.start0) if top_reg is not None else '','top_primary_end0':int(top_reg.end0) if top_reg is not None else '','top_primary_support_family_n':int(top_reg.support_family_n) if top_reg is not None else '','top_confirmation_region_id':str(top_conf.region_id) if top_conf is not None else '','top_confirmation_rank':int(top_conf.confirmation_rank) if top_conf is not None else '','truth_used':False,'postcutoff_ART_family_used':False})
    idx=root/'outputs/STAGE03A_BLIND_VALIDATION_INDEX.tsv'; pd.DataFrame(rows).to_csv(idx,sep='\t',index=False)
    gate={'gate':'PASS_STAGE03A_PREDECLARED_L3_BLIND_VALIDATION_SCORES_FROZEN' if ok else 'FAIL_STAGE03A_PREDECLARED_L3_BLIND_VALIDATION_SCORES_FROZEN','gate_pass':ok,'scientific_status':pol['scientific_status'],'blindness_level':pol['blindness_level'],'validation_accessions':[x['accession'] for x in cohort],'n_genomes':len(cohort),'truth_used_before_score_freeze':False,'postcutoff_ART_family_membership_used_before_score_freeze':False,'marsHill_used_for_validation_scoring':False,'bpe_used':False,'engine_freeze_gate_sha256':sha256_file(root/'outputs/STAGE03A_ENGINE_FREEZE_GATE.json'),'fasta_manifest_sha256':sha256_file(root/'inputs/FASTA_MANIFEST.tsv'),'blind_validation_index_sha256':sha256_file(idx),'genomes':rows,'authorization':'If PASS, Stage03A blind scores are immutable. Only now may a separate Stage03B unblinding package import post-cutoff ART truth/family annotations for evaluation. No method, threshold, window, recurrence, anchor, background, GLMY, irregularity, or BPE retuning is authorized.'}
    write_json(root/'outputs/STAGE03A_FINAL_GATE.json',gate)
    # exact file manifest for returned evidence
    files=[]
    for p in sorted(root.rglob('*')):
        if not p.is_file(): continue
        rel=p.relative_to(root).as_posix()
        if rel.startswith('work/'): continue
        if rel=='PACKAGE_MANIFEST.tsv': continue
        files.append({'path':rel,'size_bytes':p.stat().st_size,'sha256':sha256_file(p)})
    with open(root/'PACKAGE_MANIFEST.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['path','size_bytes','sha256'],delimiter='\t'); w.writeheader(); w.writerows(files)
    print(json.dumps(gate,indent=2))
    print(f'[OK] manifest files={len(files)}')
    if not ok: raise SystemExit(2)
if __name__=='__main__': main()
