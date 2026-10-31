#!/usr/bin/env python3
import argparse,glob,json,os,shutil
from pathlib import Path

def link_or_copy(src,dst):
    try: os.link(src,dst)
    except OSError: shutil.copy2(src,dst)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--shards',required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
    shards=Path(a.shards);out=Path(a.out);pred=out/'predictions'
    if out.exists(): shutil.rmtree(out)
    pred.mkdir(parents=True)
    cifs=sorted(shards.glob('gpu*/**/*_model_0.cif'))
    if len(cifs)!=40: raise SystemExit(f'expected 40 model_0 CIFs under {shards}, got {len(cifs)}')
    seen=set()
    for cif in cifs:
        q=cif.name[:-len('_model_0.cif')]
        if q in seen: raise SystemExit(f'duplicate query CIF: {q}')
        seen.add(q);srcdir=cif.parent;dstdir=pred/q;dstdir.mkdir()
        for p in sorted(srcdir.iterdir()):
            if p.is_file(): link_or_copy(p,dstdir/p.name)
    ncif=len(list(pred.glob('*/*_model_0.cif')));nconf=len(list(pred.glob('*/confidence_*_model_0.json')))
    if ncif!=40 or nconf!=40: raise SystemExit(f'canonical P2 assembly failed: CIF={ncif} confidence={nconf}')
    s={'stage':'Stage05D4-R5D-prediction','gate':'PASS_BOLTZ2_P2_40_FULL3D_SINGLEGPU_SHARDS','gate_pass':True,'n_P2':40,'n_cif':ncif,'n_confidence_json':nconf,'boltz_model':'boltz2','seed':20261007,'recycling_steps':3,'sampling_steps':200,'diffusion_samples':1,'msa_server':True,'no_kernels':True,'execution':'4 independent single-GPU shards','scientific_settings_changed':False,'candidate_filtering_performed':False,'dna_score_recomputed':False,'dna_rank_recomputed':False}
    Path('outputs').mkdir(exist_ok=True); json.dump(s,open(Path('outputs')/'GATE_STAGE05D4R5D_PREDICTION.json','w'),indent=2)
    print(json.dumps(s,indent=2))
if __name__=='__main__':main()
