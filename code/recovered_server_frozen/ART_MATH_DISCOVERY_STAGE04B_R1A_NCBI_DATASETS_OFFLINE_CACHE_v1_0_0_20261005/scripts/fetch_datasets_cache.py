#!/usr/bin/env python3
from pathlib import Path
import subprocess, json, hashlib, shutil, sys, time, zipfile

ROOT = Path(__file__).resolve().parents[1]
ACC = ROOT/'inputs/STAGE04B_R1A_ACCESSIONS_174.txt'
ZDIR = ROOT/'work/zips'
UDIR = ROOT/'work/unpacked'
OUT = ROOT/'outputs'
ZDIR.mkdir(parents=True, exist_ok=True); UDIR.mkdir(parents=True, exist_ok=True); OUT.mkdir(parents=True, exist_ok=True)

def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20), b''): h.update(b)
    return h.hexdigest()

def valid_zip(acc, z):
    if not z.exists() or z.stat().st_size == 0: return False
    try:
        with zipfile.ZipFile(z) as q:
            names=set(q.namelist())
            need={'ncbi_dataset/data/genomic.fna','ncbi_dataset/data/data_report.jsonl'}
            return need.issubset(names) and q.testzip() is None
    except Exception:
        return False

def unpack_and_check(acc,z):
    d=UDIR/acc
    if d.exists(): shutil.rmtree(d)
    d.mkdir(parents=True)
    with zipfile.ZipFile(z) as q: q.extractall(d)
    gp=d/'ncbi_dataset/data/genomic.fna'
    rp=d/'ncbi_dataset/data/data_report.jsonl'
    if not gp.exists() or not rp.exists(): return None, 'missing_required_files'
    reports=[]
    for line in rp.read_text(errors='replace').splitlines():
        line=line.strip()
        if not line: continue
        try: reports.append(json.loads(line))
        except Exception: pass
    observed=[]
    for r in reports:
        if isinstance(r,dict):
            a=r.get('accession')
            if a: observed.append(str(a))
    if acc not in observed:
        return None, 'accession_not_in_data_report:' + ','.join(observed[:5])
    return {
        'accession':acc,
        'zip_sha256':sha256(z),
        'zip_bytes':z.stat().st_size,
        'genomic_fna_sha256':sha256(gp),
        'genomic_fna_bytes':gp.stat().st_size,
        'data_report_sha256':sha256(rp),
    }, None

accs=[x.strip() for x in ACC.read_text().splitlines() if x.strip()]
if len(accs)!=174 or len(set(accs))!=174:
    raise SystemExit(f'accession registry invariant failed: n={len(accs)} unique={len(set(accs))}')

rows=[]; errors=[]
for i,acc in enumerate(accs,1):
    z=ZDIR/f'{acc}.zip'
    ok=valid_zip(acc,z)
    if not ok:
        if z.exists(): z.unlink()
        success=False
        last=''
        for attempt in range(1,6):
            cmd=['datasets','download','virus','genome','accession',acc,'--include','all','--filename',str(z)]
            p=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            last=p.stdout[-2000:]
            if p.returncode==0 and valid_zip(acc,z):
                success=True; break
            time.sleep(min(30, 5*attempt))
        if not success:
            errors.append(f'{acc}:download_failed:{last[:500]}')
            print(f'[{i}/174] FAIL {acc}', flush=True)
            continue
    rec,err=unpack_and_check(acc,z)
    if err:
        errors.append(f'{acc}:{err}')
        print(f'[{i}/174] FAIL {acc} {err}', flush=True)
        continue
    rows.append(rec)
    print(f'[{i}/174] PASS {acc} zip={rec["zip_bytes"]}', flush=True)

import csv
with open(OUT/'STAGE04B_R1A_DATASETS_CACHE_MANIFEST.tsv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['accession','zip_sha256','zip_bytes','genomic_fna_sha256','genomic_fna_bytes','data_report_sha256'],delimiter='\t')
    w.writeheader(); w.writerows(rows)

gate_pass=len(rows)==174 and not errors
gate={
 'stage':'Stage04B-R1A',
 'gate':'PASS_STAGE04B_R1A_ALL_174_NCBI_DATASETS_PACKAGES_FROZEN' if gate_pass else 'HOLD_STAGE04B_R1A_DATASETS_CACHE_INCOMPLETE',
 'gate_pass':gate_pass,
 'n_expected_sequences':174,
 'n_cache_rows':len(rows),
 'errors':errors,
 'source_registry_sha256':'73007864fdef7dba662c9b90a02e2b915359c1fbd93a6061edeea35afde07349',
 'score_recomputed':False,'region_ranking_recomputed':False,'candidate_filtering_authorized':False,'retuning_authorized':False,
 'biological_interpretation_authorized':False,
}
(OUT/'STAGE04B_R1A_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
print(json.dumps(gate,indent=2))
sys.exit(0 if gate_pass else 4)
