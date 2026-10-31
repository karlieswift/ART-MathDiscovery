#!/usr/bin/env python3
"""Fetch a fixed Pfam panel and freeze exact downloaded models before scoring."""
from __future__ import annotations
import argparse,gzip,re,time,urllib.request,urllib.error,shutil,sys
from pathlib import Path
from common import *
ROOT=Path(__file__).resolve().parents[1]

def split_hmms(text):
    blocks=[];buf=[]
    for line in text.splitlines():
        if not buf and not line.strip():continue
        buf.append(line)
        if line.strip()=='//':blocks.append('\n'.join(buf)+'\n');buf=[]
    if any(s.strip() for s in buf):raise ValueError('Truncated HMM data')
    return blocks

def hmm_meta(block):
    lines=block.splitlines()
    if not lines or not lines[0].startswith('HMMER3/'):raise ValueError('Not a HMMER3 model')
    fields={}
    for l in lines:
        p=l.split(maxsplit=1)
        if len(p)==2 and p[0] in ('ACC','NAME','LENG','ALPH','GA','DESC','DATE'):fields[p[0]]=p[1]
    if fields.get('ALPH')!='amino':raise ValueError('Model alphabet is not amino')
    if 'GA' not in fields:raise ValueError('Missing Pfam gathering thresholds')
    ga=[float(x) for x in fields['GA'].replace(';','').split()]
    if len(ga)!=2 or int(fields.get('LENG','0'))<1:raise ValueError('Bad GA/length')
    if not re.fullmatch(r'PF\d{5}\.\d+',fields.get('ACC','')):raise ValueError('Missing accession.version')
    return {'accession':fields['ACC'].split('.')[0],'accession_version':fields['ACC'],'name':fields['NAME'],'length':int(fields['LENG']),'ga_sequence_bits':ga[0],'ga_domain_bits':ga[1],'description':fields.get('DESC',''),'model_date':fields.get('DATE','')}

def unpack_payload(data):
    if data[:2]==b'\x1f\x8b':data=gzip.decompress(data)
    text=data.decode('ascii')
    if '<html' in text[:500].lower() or '<!doctype' in text[:500].lower():raise ValueError('HTML instead of HMM')
    return split_hmms(text)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--pfam-hmm',type=Path,help='Local authoritative Pfam-A.hmm[.gz]; extract only fixed accessions');a=ap.parse_args()
    pol=read_json(ROOT/'config/METHOD_POLICY.json');rf=ROOT/'refs';rf.mkdir(exist_ok=True);lockp=rf/'REFERENCE_LOCK.json';manifest=rf/'REFERENCE_MANIFEST.tsv'
    if lockp.exists():
        lock=read_json(lockp)
        if lock['method_policy_sha256']!=sha_file(ROOT/'config/METHOD_POLICY.json'):raise RuntimeError('Policy changed after reference freeze')
        assert_manifest(ROOT,manifest)
        if sha_file(manifest)!=lock['reference_manifest_sha256']:raise RuntimeError('Reference manifest changed')
        print('[PASS] existing immutable reference snapshot');return
    local={}
    if a.pfam_hmm:
        op=gzip.open if a.pfam_hmm.suffix=='.gz' else open
        with op(a.pfam_hmm,'rt',encoding='ascii') as f:
            buf=[]
            for line in f:
                buf.append(line)
                if line.strip()=='//':
                    b=''.join(buf);buf=[]
                    ma=re.search(r'^ACC\s+(PF\d{5})\.',b,re.M)
                    if ma and ma[1] in {p['accession'] for p in pol['rt_scan']['profiles']}:local[ma[1]]=b
    models=[];combined=[]
    for item in pol['rt_scan']['profiles']:
        acc=item['accession'];p=rf/(acc+'.hmm');url=f'https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc}/?annotation=hmm'
        if p.exists():blocks=split_hmms(p.read_text('ascii'));source='EXISTING_LOCAL_PANEL_FILE'
        elif a.pfam_hmm:
            if acc not in local:raise RuntimeError('Fixed model missing in local Pfam library: '+acc)
            blocks=[local[acc]];source='LOCAL_PFAM_FILE_SHA256:'+sha_file(a.pfam_hmm)
        else:
            last=None
            for attempt in range(3):
                try:
                    req=urllib.request.Request(url,headers={'User-Agent':'ART-MathDiscovery-RT-evidence/1.0','Accept':'*/*'})
                    with urllib.request.urlopen(req,timeout=60) as r:
                        data=r.read(10*1024*1024+1)
                        if len(data)>10*1024*1024:raise ValueError('Unexpectedly large per-model response')
                    blocks=unpack_payload(data);source=url;break
                except Exception as e:
                    last=e;print(f'[REFERENCE] {acc} attempt={attempt+1} {type(e).__name__}: {e}',flush=True)
                    if attempt<2:time.sleep(5*(attempt+1))
            else:raise RuntimeError('Reference fetch failed before scoring: '+str(last))
        if len(blocks)!=1:raise ValueError('Expected exactly one model for '+acc)
        meta=hmm_meta(blocks[0])
        if meta['accession']!=acc or meta['name']!=item['name']:raise ValueError('Model identity changed: '+str(meta))
        if not p.exists():atomic_text(p,blocks[0])
        meta.update({'role':item['role'],'source':source,'sha256':sha_file(p),'path':str(p.relative_to(ROOT))});models.append(meta);combined.append(p.read_text('ascii'))
        print('[REFERENCE FROZEN]',meta['accession_version'],meta['sha256'],flush=True)
    atomic_text(rf/'RT_PANEL.hmm',''.join(combined))
    write_tsv(rf/'MODEL_CATALOG.tsv',models,list(models[0]))
    mh=output_manifest(ROOT,[rf/(p['accession']+'.hmm') for p in pol['rt_scan']['profiles']]+[rf/'RT_PANEL.hmm',rf/'MODEL_CATALOG.tsv'],manifest)
    write_json(lockp,{'created_utc':now(),'method_policy_sha256':sha_file(ROOT/'config/METHOD_POLICY.json'),'reference_manifest_sha256':mh,'models':models,'reference_scope':'three-model panel, not all RT diversity','scores_observed_in_this_stage':False})
    print('[PASS] reference models byte-frozen before scoring')

if __name__=='__main__':
    try:main()
    except Exception as e:print('[HOLD] '+str(e),file=sys.stderr);sys.exit(2)
