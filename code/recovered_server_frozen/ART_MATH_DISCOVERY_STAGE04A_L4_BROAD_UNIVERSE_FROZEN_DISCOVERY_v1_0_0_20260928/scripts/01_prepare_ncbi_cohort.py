#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, hashlib, json, re
from common import safe_id, accession_root

EXCLUDE_ROOTS={'MW248466','MW349129','MW349128','MW218148'}

def iter_fasta(path):
    h=None; seq=[]
    with open(path) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if h is not None: yield h,''.join(seq).upper()
                h=line[1:].strip(); seq=[]
            else: seq.append(line)
    if h is not None: yield h,''.join(seq).upper()

def normkey(k): return re.sub('[^a-z0-9]','',k.lower())
def find_value(obj, keys):
    want={normkey(x) for x in keys}
    if isinstance(obj,dict):
        for k,v in obj.items():
            if normkey(k) in want and isinstance(v,(str,int,float)): return str(v)
        for v in obj.values():
            x=find_value(v,keys)
            if x: return x
    elif isinstance(obj,list):
        for v in obj:
            x=find_value(v,keys)
            if x: return x
    return ''

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--genomic-fna',required=True); ap.add_argument('--data-report-jsonl')
    ap.add_argument('--arm',required=True,choices=['L4_CURATED','L4_PROSPECTIVE','L4_METAVIROME'])
    ap.add_argument('--source',default='NCBI_DATASETS'); ap.add_argument('--outdir',required=True)
    ap.add_argument('--sample-n',type=int,default=0,help='0=all; otherwise deterministic accession-hash sample')
    a=ap.parse_args(); out=Path(a.outdir).resolve(); fastadir=out/'prepared_fastas'; fastadir.mkdir(parents=True,exist_ok=True)
    meta={}
    if a.data_report_jsonl and Path(a.data_report_jsonl).exists():
        for line in Path(a.data_report_jsonl).read_text().splitlines():
            if not line.strip(): continue
            try: o=json.loads(line)
            except Exception: continue
            acc=find_value(o,['accession','nucleotideAccession'])
            if not acc: continue
            meta[acc.split('.')[0]]={
                'accession':acc,
                'release_date':find_value(o,['releaseDate','release_date','release-date']),
                'virus_name':find_value(o,['virusName','organismName','name'])
            }
    records=[]
    for header,seq in iter_fasta(a.genomic_fna):
        acc=header.split()[0]
        root=accession_root(acc)
        if root in EXCLUDE_ROOTS: continue
        if len(seq)<4096 or (set(seq)-set('ACGT')): continue
        m=meta.get(root,{})
        records.append((acc,header,seq,m.get('release_date',''),m.get('virus_name','')))
    if a.sample_n and len(records)>a.sample_n:
        records=sorted(records,key=lambda r:hashlib.sha256(('STAGE04A_SAMPLE_V1|'+r[0]).encode()).hexdigest())[:a.sample_n]
    rows=[]
    for acc,header,seq,rd,vname in sorted(records,key=lambda r:r[0]):
        aid=safe_id(acc); fp=fastadir/f'{aid}.fasta'
        fp.write_text('>'+header+'\n'+'\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')
        rows.append({'sequence_id':acc,'arm':a.arm,'source':a.source,'release_date':rd,'fasta_path':str(fp),'cluster_id':'','notes':vname})
    cohort=out/f'STAGE04A_{a.arm}_COHORT.tsv'
    with open(cohort,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['sequence_id','arm','source','release_date','fasta_path','cluster_id','notes'],delimiter='\t'); w.writeheader(); w.writerows(rows)
    print(json.dumps({'cohort_tsv':str(cohort),'n_sequences':len(rows),'arm':a.arm,'excluded_known_roots':sorted(EXCLUDE_ROOTS)},indent=2))
if __name__=='__main__': main()
