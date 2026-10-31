#!/usr/bin/env python3
import argparse,csv,hashlib,json,os,shutil
from pathlib import Path
TARGET={'A':13_500_000_000,'B':18_000_000_000,'C':13_500_000_000,'D':9_000_000_000}
OUTFIELDS=['stratum','accession','nucleotide_bp','contig_n','compressed_bytes','zst_sha256','s3_uri']
def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def read_tsv(p):
    with open(p,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def link_or_keep(src,dst):
    src=Path(src).resolve(); dst=Path(dst); dst.parent.mkdir(parents=True,exist_ok=True)
    if dst.exists() or dst.is_symlink():
        try:
            if os.path.samefile(src,dst): return
        except Exception: pass
        if dst.is_file() or dst.is_symlink(): dst.unlink()
    try: os.link(src,dst)
    except OSError:
        try: os.symlink(src,dst)
        except OSError: shutil.copy2(src,dst)
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--stage',required=True); ap.add_argument('--source-r1-stage',required=True); a=ap.parse_args()
    S=Path(a.stage).resolve(); R=Path(a.source_r1_stage).resolve(); old={x.strip() for x in open(R/'inputs/OLD_F3_ACCESSIONS.txt') if x.strip() and not x.startswith('#')}
    allrows=[]; by={}; seen=set()
    for s in 'ABCD':
        p=S/'frozen/shards'/f'{s}_SELECTED.tsv'; m=S/'frozen/shards'/f'{s}_FREEZE.json'
        if not p.is_file() or not m.is_file(): raise SystemExit(f'FATAL: missing completed shard {s}; run 02_acquire_stratum.sh {s} first')
        rows=sorted(read_tsv(p),key=lambda r:int(r['deterministic_rank']))
        bp=sum(int(r['nucleotide_bp']) for r in rows); by[s]=bp
        if bp<TARGET[s]: raise SystemExit(f'FATAL: shard {s} below target {bp} < {TARGET[s]}')
        for r in rows:
            acc=r['accession']
            if acc in old: raise SystemExit(f'FATAL old-universe overlap: {acc}')
            if acc in seen: raise SystemExit(f'FATAL duplicate accepted accession across shards: {acc}')
            seen.add(acc); allrows.append(r)
            src=S/'data/shards'/s/f'{acc}.contigs.fa.zst'
            if not src.is_file(): raise SystemExit(f'FATAL missing shard contig: {src}')
            if sha256(src)!=r['zst_sha256']: raise SystemExit(f'FATAL shard SHA mismatch: {acc}')
            link_or_keep(src,S/'data/contigs'/f'{acc}.contigs.fa.zst')
    final=S/'frozen/F3V2_DELTA_SELECTED_ACCESSIONS.tsv'; final.parent.mkdir(parents=True,exist_ok=True)
    with open(final,'w',newline='') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=OUTFIELDS); w.writeheader();
        for r in allrows: w.writerow({k:r[k] for k in OUTFIELDS})
    # concatenate ledgers with source shard identity retained in detail; deterministic A-D order
    ledger_out=S/'frozen/F3V2_ACQUISITION_LEDGER.tsv'
    lf=['event_id','stratum','deterministic_rank','accession','status','compressed_bytes','nucleotide_bp','contig_n','detail']
    eid=0
    with open(ledger_out,'w',newline='') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=lf); w.writeheader()
        for s in 'ABCD':
            for r in read_tsv(S/'work/shards'/s/'ledger.tsv'):
                eid+=1; r['event_id']=eid; w.writerow({k:r.get(k,'') for k in lf})
    total=sum(by.values()); gate=50_000_000_000<=total<=60_000_000_000
    g={'stage':'Stage05F3V2-R2-distributed-source-freeze','gate':'PASS_F3V2_R2_DISTRIBUTED_54G_SOURCE_BLIND_FREEZE' if gate else 'FAIL_F3V2_R2_TOTAL_BP_GATE','gate_pass':gate,'delta_accession_n':len(allrows),'delta_nucleotide_bp':total,'stratum_bp':by,'stratum_targets_bp':TARGET,'old_accession_overlap_n':0,'source_selection_used_rt_or_candidate_information':False,'source_r1_stage':str(R),'selected_manifest_sha256':sha256(final),'acquisition_ledger_sha256':sha256(ledger_out),'science_change':False,'operational_change':'A/B/C/D acquired independently and concurrently; each stratum preserves the frozen R1 deterministic candidate order and stop-at-quota semantics; final merge is A-D then rank order.','claim_boundary':'Distributed acquisition only. It does not change the source universe, frozen quotas, technical filters, six-frame translation, HMM profiles, clustering, repeat rules, structural thresholds, or novelty rules.'}
    (S/'frozen/F3V2_DELTA_SOURCE_FREEZE.json').write_text(json.dumps(g,indent=2)+'\n'); print(json.dumps(g,indent=2));
    if not gate: raise SystemExit(2)
    print('PASS_F3V2_R2_DISTRIBUTED_54G_DOWNLOAD_AND_SOURCE_FREEZE_COMPLETE')
if __name__=='__main__': main()
