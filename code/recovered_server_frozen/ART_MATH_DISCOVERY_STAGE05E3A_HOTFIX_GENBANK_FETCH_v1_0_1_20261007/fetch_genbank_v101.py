#!/usr/bin/env python3
import argparse, hashlib, io, json, os, time, urllib.parse, urllib.request
from pathlib import Path
from Bio import SeqIO

EUTILS='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
SVIEW='https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi'

def normalize_text(data: bytes) -> str:
    return data.decode('utf-8', errors='replace').lstrip('\ufeff\r\n\t ')

def valid_gb_bytes(data: bytes):
    if not data or len(data) < 300:
        return False, 'too_short_or_empty'
    text = normalize_text(data)
    if not text.startswith('LOCUS'):
        return False, 'does_not_start_LOCUS'
    if '\nFEATURES' not in text:
        return False, 'missing_FEATURES'
    if '\n//' not in text:
        return False, 'missing_record_terminator'
    try:
        recs = list(SeqIO.parse(io.StringIO(text), 'genbank'))
    except Exception as e:
        return False, 'biopython_parse_error:'+repr(e)
    if len(recs) < 1:
        return False, 'biopython_zero_records'
    return True, f'parsed_records={len(recs)}'

def valid_gb_file(p: Path):
    try:
        ok, why = valid_gb_bytes(p.read_bytes())
        return ok, why
    except Exception as e:
        return False, repr(e)

def route_urls(acc, api_key='', email=''):
    base_id = acc.split('.')[0]
    ids = [acc] + ([base_id] if base_id != acc else [])
    routes=[]
    # Official EFetch documentation uses rettype=gb for GenBank flat files.
    for ident in ids:
        params={'db':'nuccore','id':ident,'rettype':'gb','retmode':'text','tool':'ART_MathDiscovery'}
        if api_key: params['api_key']=api_key
        if email: params['email']=email
        routes.append((f'EFETCH_GB:{ident}', EUTILS+'?'+urllib.parse.urlencode(params)))
    # Independent NCBI Sequence Viewer text route as transport fallback.
    for ident in ids:
        params={'id':ident,'db':'nuccore','report':'genbank','retmode':'text'}
        routes.append((f'SVIEW_GENBANK:{ident}', SVIEW+'?'+urllib.parse.urlencode(params)))
    return routes

def fetch_url(url, ua, timeout=90):
    req=urllib.request.Request(url, headers={
        'User-Agent': ua,
        'Accept': 'text/plain,*/*;q=0.1',
        'Cache-Control': 'no-cache',
    })
    with urllib.request.urlopen(req, timeout=timeout) as r:
        data=r.read()
        status=getattr(r,'status',None) or r.getcode()
        ctype=r.headers.get('Content-Type','')
    return data, status, ctype

def safe_first_line(data):
    if not data: return ''
    t=data.decode('utf-8',errors='replace').replace('\r','\n')
    line=next((x.strip() for x in t.split('\n') if x.strip()), '')
    return line[:240].replace('\t',' ')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--ids', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--diagnostics', required=True)
    ap.add_argument('--max-route-attempts', type=int, default=2)
    a=ap.parse_args()

    ids=[x.strip() for x in Path(a.ids).read_text().splitlines() if x.strip()]
    out=Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    api=os.environ.get('NCBI_API_KEY','').strip()
    email=os.environ.get('NCBI_EMAIL','').strip()
    ua='ART-MathDiscovery/Stage05E3A-v1.0.1'
    if email: ua += f' ({email})'
    inter_request_sleep=0.11 if api else 0.36

    rows=[]; diag=[]; failed=[]
    for i,acc in enumerate(ids,1):
        p=out/f'{acc}.gbff'
        ok_existing, why_existing = valid_gb_file(p) if p.exists() else (False,'missing')
        selected_route=''
        if ok_existing:
            status='CACHED_VALID'
            selected_route='EXISTING_VALID'
        else:
            ok=False; last_err=why_existing
            for route,url in route_urls(acc,api,email):
                for attempt in range(1,a.max_route_attempts+1):
                    t0=time.time(); data=b''; code=''; ctype=''; err=''
                    try:
                        data,code,ctype=fetch_url(url,ua)
                        ok,why=valid_gb_bytes(data)
                        if not ok: err=why
                    except Exception as e:
                        ok=False; err=repr(e)
                    diag.append({
                        'sequence_id':acc,'route':route,'attempt':attempt,'http_status':code,
                        'content_type':ctype,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest() if data else '',
                        'valid_genbank':int(ok),'validation_or_error':('PASS' if ok else err),
                        'first_nonempty_line':safe_first_line(data),'elapsed_s':round(time.time()-t0,3)
                    })
                    if ok:
                        tmp=p.with_suffix(p.suffix+'.tmp')
                        tmp.write_bytes(data)
                        ok2,why2=valid_gb_file(tmp)
                        if not ok2:
                            tmp.unlink(missing_ok=True)
                            ok=False; last_err='atomic_revalidation_failed:'+why2
                        else:
                            tmp.replace(p); selected_route=route; status='FETCHED'; break
                    else:
                        last_err=f'{route} attempt {attempt}: {err}'
                    time.sleep(inter_request_sleep + min(2.0,0.4*attempt))
                if ok: break
            if not ok:
                status='FAILED'; failed.append((acc,last_err))
        sha=hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else ''
        size=p.stat().st_size if p.exists() else 0
        rows.append((acc,status,size,sha,selected_route))
        print(f'[{i}/{len(ids)}] {acc} {status} route={selected_route or "NONE"}', flush=True)
        time.sleep(inter_request_sleep)

    with open(a.manifest,'w') as h:
        h.write('sequence_id\tstatus\tbytes\tsha256\tselected_route\n')
        for r in rows: h.write('\t'.join(map(str,r))+'\n')
    with open(a.diagnostics,'w') as h:
        cols=['sequence_id','route','attempt','http_status','content_type','bytes','sha256','valid_genbank','validation_or_error','first_nonempty_line','elapsed_s']
        h.write('\t'.join(cols)+'\n')
        for d in diag:
            h.write('\t'.join(str(d[c]).replace('\t',' ') for c in cols)+'\n')

    summary={
        'requested':len(ids),
        'valid':sum(r[1] in ('CACHED_VALID','FETCHED') for r in rows),
        'failed':len(failed),
        'failed_ids':[x[0] for x in failed],
        'transport_hotfix':'v1.0.1: use documented rettype=gb EFetch first, accession-without-version retry, NCBI sviewer GenBank fallback, atomic validation, diagnostic capture',
        'scientific_settings_changed':False,
    }
    print(json.dumps(summary,indent=2))
    if failed:
        raise SystemExit('GenBank retrieval incomplete after all routes: '+repr(failed[:5]))

if __name__=='__main__':
    main()
