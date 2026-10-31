#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,time,urllib.parse,urllib.request
from pathlib import Path
from common import *

def fetch_fasta(acc):
    urls=[
      'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urllib.parse.urlencode({'db':'nuccore','id':acc,'rettype':'fasta','retmode':'text'}),
      'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?'+urllib.parse.urlencode({'id':acc,'db':'nuccore','report':'fasta','retmode':'text'})]
    errs=[]
    for url in urls:
        try:
            req=urllib.request.Request(url,headers={'User-Agent':'ART-MathDiscovery/Stage03A'})
            with urllib.request.urlopen(req,timeout=45) as r: txt=r.read().decode('utf-8','replace')
            if txt.lstrip().startswith('>') and len(txt)>1000: return txt,url
            errs.append('nonfasta response')
        except Exception as e: errs.append(repr(e))
        time.sleep(1)
    raise RuntimeError(f'Could not fetch {acc}: {errs}')

def parse_text(txt):
    lines=[x.strip() for x in txt.splitlines() if x.strip()]
    if not lines or not lines[0].startswith('>'): raise ValueError('bad FASTA')
    header=lines[0][1:].split()[0]; seq=''.join(lines[1:]).upper()
    return header,seq

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); a=ap.parse_args(); root=Path(a.root).resolve(); rows=[]; ok=True
    with open(root/'config/VALIDATION_COHORT.tsv') as f: cohort=list(csv.DictReader(f,delimiter='\t'))
    for r in cohort:
        acc=r['accession']; p=root/'inputs'/f'{acc}.frozen.fasta'; source='preexisting'; url=''
        if p.exists(): header,seq=parse_fasta(p)
        else:
            txt,url=fetch_fasta(acc); header,seq=parse_text(txt); source='NCBI_nuccore_raw_fasta'; write_fasta(p,header,seq)
        alphabet=''.join(sorted(set(seq))); valid=set(seq)<=set('ACGT') and len(seq)>=4096
        if not valid: ok=False
        rows.append({'requested_accession':acc,'fasta_header':header,'length_bp':len(seq),'alphabet':alphabet,'fasta_sha256':sha256_file(p),'sequence_sha256':__import__('hashlib').sha256(seq.encode()).hexdigest(),'source':source,'fetch_url':url,'precutoff_annotation_used':False,'postcutoff_truth_used':False,'valid_for_frozen_engine':valid})
    op=root/'inputs/FASTA_MANIFEST.tsv'
    with open(op,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
    gate={'gate':'PASS_STAGE03A_VALIDATION_FASTAS_FROZEN' if ok else 'FAIL_STAGE03A_VALIDATION_FASTAS_FROZEN','gate_pass':ok,'n_genomes':len(rows),'manifest_sha256':sha256_file(op),'genomes':rows,'truth_used':False,'annotation_used':False}
    write_json(root/'outputs/STAGE03A_FASTA_FREEZE_GATE.json',gate); print(json.dumps(gate,indent=2))
    if not ok: raise SystemExit('Validation requires one ACGT-only >=4096bp FASTA per accession. If N is present, stop rather than invent a new masking policy.')
if __name__=='__main__': main()
