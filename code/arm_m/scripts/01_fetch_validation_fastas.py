#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,hashlib,json,os,shutil,time,urllib.parse,urllib.request
from pathlib import Path
import sys

# This file replaces only the Stage03A transport layer. It must be installed
# inside the frozen Stage03A package, where common.py remains unchanged.
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from common import *

USER_AGENT='ART-MathDiscovery/Stage03A-R1 raw-fasta-only'

class FetchFailure(RuntimeError):
    def __init__(self,message,attempts):
        super().__init__(message)
        self.attempts=attempts

def candidate_ids(acc:str):
    out=[acc]
    # Stage00 contains two unversioned accessions. Trying .1 is a transport
    # resolution only; it does not use annotation or ART truth.
    if '.' not in acc:
        out.append(acc+'.1')
    return out

def route_specs(acc:str):
    q1=urllib.parse.urlencode({'db':'nuccore','id':acc,'rettype':'fasta','retmode':'text','tool':'ART_MathDiscovery_Stage03A'})
    q2=urllib.parse.urlencode({'id':acc,'db':'nuccore','report':'fasta','retmode':'text'})
    return [
        ('NCBI_EFETCH', 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+q1),
        ('NCBI_SVIEWER','https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?'+q2),
        ('ENA_BROWSER_API','https://www.ebi.ac.uk/ena/browser/api/fasta/'+urllib.parse.quote(acc,safe='')),
    ]

def parse_text(txt:str):
    lines=[x.strip() for x in txt.splitlines() if x.strip()]
    if not lines or not lines[0].startswith('>'):
        raise ValueError('not FASTA')
    if any(x.startswith('>') for x in lines[1:]):
        raise ValueError('multi-record FASTA is not allowed')
    header=lines[0][1:].strip()
    seq=''.join(lines[1:]).upper()
    if not seq:
        raise ValueError('empty FASTA sequence')
    return header,seq

def preview(txt:str,n=96):
    return ' '.join(txt[:n].replace('\r',' ').replace('\n',' ').split())

def http_fetch(url:str,timeout=45):
    req=urllib.request.Request(url,headers={
        'User-Agent':USER_AGENT,
        'Accept':'text/plain,text/x-fasta,*/*;q=0.1',
        'Cache-Control':'no-cache',
    })
    with urllib.request.urlopen(req,timeout=timeout) as r:
        raw=r.read()
        status=getattr(r,'status',None)
        ctype=r.headers.get('Content-Type','') if getattr(r,'headers',None) else ''
    return raw.decode('utf-8','replace'),status,ctype

def fetch_fasta(acc:str):
    attempts=[]
    for resolved in candidate_ids(acc):
        for route,url in route_specs(resolved):
            rec={'requested_accession':acc,'resolved_query':resolved,'route':route,'url':url}
            try:
                txt,status,ctype=http_fetch(url)
                rec.update({'http_status':status,'content_type':ctype,'bytes':len(txt.encode('utf-8','replace')),'preview':preview(txt)})
                try:
                    header,seq=parse_text(txt)
                    rec.update({'result':'ACCEPTED_FASTA','length_bp':len(seq)})
                    attempts.append(rec)
                    return txt,url,route,resolved,attempts
                except Exception as e:
                    rec.update({'result':'REJECTED_NONFASTA','error':repr(e)})
            except Exception as e:
                rec.update({'result':'FETCH_ERROR','error':repr(e)})
            attempts.append(rec)
            time.sleep(0.5)
    raise FetchFailure('Could not fetch raw FASTA for %s via NCBI/ENA routes' % acc,attempts)

def local_candidates(source_dir:Path,acc:str):
    names=[]
    for rid in candidate_ids(acc):
        for suf in ('.frozen.fasta','.fasta','.fna','.fa'):
            names.append(source_dir/(rid+suf))
    return names

def read_local_source(source_dir:Path,acc:str):
    for p in local_candidates(source_dir,acc):
        if p.is_file():
            header,seq=parse_fasta(p)
            return p,header,seq
    return None

def validate_seq(seq:str):
    alphabet=''.join(sorted(set(seq)))
    valid=set(seq)<=set('ACGT') and len(seq)>=4096
    return valid,alphabet

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',required=True)
    ap.add_argument('--source-dir',default='')
    a=ap.parse_args()
    root=Path(a.root).resolve()
    source_dir=Path(a.source_dir).resolve() if a.source_dir else None
    (root/'inputs').mkdir(parents=True,exist_ok=True)
    (root/'outputs').mkdir(parents=True,exist_ok=True)
    (root/'logs').mkdir(parents=True,exist_ok=True)
    with open(root/'config/VALIDATION_COHORT.tsv') as f:
        cohort=list(csv.DictReader(f,delimiter='\t'))
    rows=[]; all_attempts=[]; ok=True
    for r in cohort:
        acc=r['accession']; p=root/'inputs'/f'{acc}.frozen.fasta'; source=''; url=''; route=''; resolved=acc
        if p.exists():
            header,seq=parse_fasta(p)
            source='preexisting_frozen_fasta'
            route='LOCAL_PREEXISTING'
        else:
            loc=read_local_source(source_dir,acc) if source_dir else None
            if loc is not None:
                src,header,seq=loc
                valid,alphabet=validate_seq(seq)
                if not valid:
                    raise SystemExit(f'Local source for {acc} is invalid: {src}; requires one ACGT-only >=4096bp FASTA.')
                write_fasta(p,header,seq)
                source='local_raw_fasta_source_dir'
                route='LOCAL_SOURCE_DIR'
                url=str(src)
            else:
                try:
                    txt,url,route,resolved,attempts=fetch_fasta(acc)
                    all_attempts.extend(attempts)
                except FetchFailure as e:
                    attempts=e.attempts
                    all_attempts.extend(attempts)
                    # Write transport attempts before failing so diagnosis is preserved.
                    write_attempts(root/'logs/STAGE03A_FASTA_TRANSPORT_ATTEMPTS.tsv',all_attempts)
                    msg=(f'Could not obtain raw FASTA for {acc}. No scientific scoring was run for this accession. '
                         f'Place an ACGT-only raw FASTA in inputs/{acc}.frozen.fasta or rerun with a local source directory as third runner argument.')
                    raise SystemExit(msg)
                header,seq=parse_text(txt)
                valid,alphabet=validate_seq(seq)
                if not valid:
                    raise SystemExit(f'Fetched FASTA for {acc} is not ACGT-only >=4096bp; stop rather than invent a masking policy.')
                write_fasta(p,header,seq)
                source='raw_fasta_transport'
        valid,alphabet=validate_seq(seq)
        if not valid: ok=False
        rows.append({
            'requested_accession':acc,'resolved_query':resolved,'fasta_header':header,'length_bp':len(seq),
            'alphabet':alphabet,'fasta_sha256':sha256_file(p),
            'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),
            'source':source,'transport_route':route,'fetch_url_or_local_source':url,
            'precutoff_annotation_used':False,'postcutoff_truth_used':False,'valid_for_frozen_engine':valid})
    write_attempts(root/'logs/STAGE03A_FASTA_TRANSPORT_ATTEMPTS.tsv',all_attempts)
    op=root/'inputs/FASTA_MANIFEST.tsv'
    with open(op,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
    gate={
        'gate':'PASS_STAGE03A_VALIDATION_FASTAS_FROZEN' if ok else 'FAIL_STAGE03A_VALIDATION_FASTAS_FROZEN',
        'gate_pass':ok,'hotfix':'STAGE03A_R1_FETCH_RESILIENCE','n_genomes':len(rows),
        'manifest_sha256':sha256_file(op),'genomes':rows,'truth_used':False,'annotation_used':False,
        'transport_policy':'PREEXISTING_OR_LOCAL_RAW_FASTA_OR_NCBI_EFETCH_OR_NCBI_SVIEWER_OR_ENA_RAW_FASTA__FAIL_CLOSED'}
    write_json(root/'outputs/STAGE03A_FASTA_FREEZE_GATE.json',gate)
    print(json.dumps(gate,indent=2))
    if not ok:
        raise SystemExit('Validation requires one ACGT-only >=4096bp FASTA per accession. If N is present, stop rather than invent a masking policy.')

def write_attempts(path:Path,rows):
    fields=['requested_accession','resolved_query','route','url','http_status','content_type','bytes','result','length_bp','preview','error']
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(rows)

if __name__=='__main__': main()
