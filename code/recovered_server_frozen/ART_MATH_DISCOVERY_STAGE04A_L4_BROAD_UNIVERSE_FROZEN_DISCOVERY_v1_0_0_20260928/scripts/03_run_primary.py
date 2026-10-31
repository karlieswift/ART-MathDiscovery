#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, json, shutil, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from common import write_json

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); ap.add_argument('--genome-workers',type=int,default=2); ap.add_argument('--inner-threads',type=int,default=16); ap.add_argument('--keep-raw',type=int,default=0); a=ap.parse_args()
    root=Path(a.root).resolve(); engkids=[p for p in (root/'runtime_engine').iterdir() if p.is_dir()]
    if len(engkids)!=1: raise SystemExit('run preflight first')
    eng=engkids[0]; manifest=root/'outputs/STAGE04A_UNIVERSE_MANIFEST.tsv'
    with open(manifest) as f: rows=list(csv.DictReader(f,delimiter='\t'))
    gout=root/'outputs/genomes'; gout.mkdir(parents=True,exist_ok=True); logs=root/'logs'; logs.mkdir(exist_ok=True)
    def one(r):
        aid=r['analysis_id']; t=time.time(); cmd=[sys.executable,str(eng/'scripts/02_primary_scan.py'),'--root',str(eng),'--accession',aid,'--threads',str(a.inner_threads),'--keep-raw',str(a.keep_raw)]
        cp=subprocess.run(cmd,text=True,capture_output=True)
        (logs/f'03_primary_{aid}.log').write_text(cp.stdout+'\n'+cp.stderr)
        if cp.returncode!=0: return aid,False,time.time()-t,cp.returncode
        src=eng/'outputs'/aid; dst=gout/aid
        if dst.exists(): shutil.rmtree(dst)
        shutil.copytree(src,dst)
        return aid,True,time.time()-t,0
    results=[]
    with ThreadPoolExecutor(max_workers=max(1,a.genome_workers)) as ex:
        futs={ex.submit(one,r):r for r in rows}
        for fut in as_completed(futs):
            z=fut.result(); results.append(z); print(f'[{z[0]}] pass={z[1]} elapsed={z[2]:.1f}s',flush=True)
    failed=[z for z in results if not z[1]]
    gate={'gate':'PASS_STAGE04A_ALL_PRIMARY_BOARDS_FROZEN' if not failed else 'FAIL_STAGE04A_PRIMARY_SCAN','gate_pass':not failed,'n_sequences':len(rows),'n_success':len(rows)-len(failed),'n_failed':len(failed),'genome_workers':a.genome_workers,'inner_threads':a.inner_threads,'keep_raw':a.keep_raw,'truth_used':False,'retuning_authorized':False,'failed':[z[0] for z in failed]}
    write_json(root/'outputs/STAGE04A_PRIMARY_SCAN_GATE.json',gate); print(json.dumps(gate,indent=2))
    if failed: raise SystemExit(2)
if __name__=='__main__': main()
