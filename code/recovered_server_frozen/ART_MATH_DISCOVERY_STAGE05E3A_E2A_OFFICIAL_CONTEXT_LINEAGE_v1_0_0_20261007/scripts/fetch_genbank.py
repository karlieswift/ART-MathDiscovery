#!/usr/bin/env python3
import argparse, os, time, urllib.parse, urllib.request, hashlib, json
from pathlib import Path

def valid_gb(p):
    try:
        t=p.read_text(errors='ignore')
        return t.startswith('LOCUS') and '\n//' in t and len(t)>500
    except Exception:
        return False

def fetch(acc, out, api_key='', email=''):
    base='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    params={'db':'nuccore','id':acc,'rettype':'gbwithparts','retmode':'text','tool':'ART_MathDiscovery'}
    if api_key: params['api_key']=api_key
    if email: params['email']=email
    url=base+'?'+urllib.parse.urlencode(params)
    req=urllib.request.Request(url,headers={'User-Agent':'ART-MathDiscovery/Stage05E3A'})
    with urllib.request.urlopen(req,timeout=90) as r:
        data=r.read()
    out.write_bytes(data)
    return valid_gb(out)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--ids',required=True)
    ap.add_argument('--outdir',required=True)
    ap.add_argument('--manifest',required=True)
    a=ap.parse_args()
    ids=[x.strip() for x in Path(a.ids).read_text().splitlines() if x.strip()]
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    api=os.environ.get('NCBI_API_KEY','').strip(); email=os.environ.get('NCBI_EMAIL','').strip()
    sleep_s=0.12 if api else 0.38
    rows=[]; failed=[]
    for i,acc in enumerate(ids,1):
        p=out/f'{acc}.gbff'
        if valid_gb(p):
            status='CACHED_VALID'
        else:
            ok=False; err=''
            for attempt in range(1,6):
                try:
                    ok=fetch(acc,p,api,email)
                    if not ok and '.' in acc:
                        ok=fetch(acc.split('.')[0],p,api,email)
                    if ok: break
                    err='invalid GenBank payload'
                except Exception as e:
                    err=repr(e)
                time.sleep(min(10,attempt*2))
            if not ok:
                failed.append((acc,err)); status='FAILED'
            else: status='FETCHED'
        sha=hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else ''
        rows.append((acc,status,p.stat().st_size if p.exists() else 0,sha))
        print(f'[{i}/{len(ids)}] {acc} {status}',flush=True)
        time.sleep(sleep_s)
    with open(a.manifest,'w') as h:
        h.write('sequence_id\tstatus\tbytes\tsha256\n')
        for r in rows: h.write('\t'.join(map(str,r))+'\n')
    summary={'requested':len(ids),'valid':sum(r[1] in ('CACHED_VALID','FETCHED') for r in rows),'failed':len(failed),'failed_ids':[x[0] for x in failed]}
    print(json.dumps(summary,indent=2))
    if failed: raise SystemExit('GenBank retrieval incomplete: '+repr(failed[:5]))
if __name__=='__main__': main()
