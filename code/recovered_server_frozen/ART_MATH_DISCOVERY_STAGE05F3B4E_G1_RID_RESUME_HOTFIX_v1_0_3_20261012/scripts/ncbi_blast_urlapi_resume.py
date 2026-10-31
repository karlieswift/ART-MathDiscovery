#!/usr/bin/env python3
import argparse, os, re, sys, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ENDPOINT='https://blast.ncbi.nlm.nih.gov/Blast.cgi'
USER_AGENT='ART-MathDiscovery-F3B4E/1.0.3'
TOOL='ART_MathDiscovery_F3B4E'

def read_text(p): return Path(p).read_text(encoding='utf-8')
def local(tag): return tag.rsplit('}',1)[-1]
def child_text(node,name,default=''):
    for c in list(node):
        if local(c.tag)==name: return (c.text or '').strip()
    return default
def descendants(node,name): return [x for x in node.iter() if local(x.tag)==name]

def request(params, timeout=180):
    p=dict(params); p.setdefault('TOOL',TOOL)
    email=os.environ.get('NCBI_EMAIL','').strip()
    if email: p.setdefault('EMAIL',email)
    data=urllib.parse.urlencode(p).encode('utf-8')
    req=urllib.request.Request(ENDPOINT,data=data,headers={'User-Agent':USER_AGENT,'Content-Type':'application/x-www-form-urlencoded'})
    with urllib.request.urlopen(req,timeout=timeout) as r:
        return r.read().decode('utf-8','replace')

def submit(program,database,query_text,expect,hitlist,filter_value,megablast=False,no_db_override=False):
    p={'CMD':'Put','PROGRAM':program,'DATABASE':database,'QUERY':query_text,'EXPECT':str(expect),'HITLIST_SIZE':str(hitlist),'FILTER':filter_value}
    if megablast: p['MEGABLAST']='on'
    if no_db_override: p['NO_DATABASE_OVERRIDE']='true'
    body=request(p)
    m=re.search(r'RID\s*=\s*([A-Za-z0-9_.-]+)',body); t=re.search(r'RTOE\s*=\s*(\d+)',body)
    if not m: raise RuntimeError('URL API submission returned no RID: '+body[:1000].replace('\n',' '))
    return m.group(1), int(t.group(1)) if t else 60

def search_info(rid):
    body=request({'CMD':'Get','RID':rid,'FORMAT_OBJECT':'SearchInfo'})
    m=re.search(r'Status\s*=\s*([A-Za-z]+)',body,re.I)
    status=(m.group(1).upper() if m else 'UNKNOWN_PARSE')
    return status,body

def fetch_xml(rid):
    return request({'CMD':'Get','RID':rid,'FORMAT_TYPE':'XML2','DESCRIPTIONS':'200','ALIGNMENTS':'200'})

def wait_existing_or_new(rid,rtoe,max_wait):
    wait=max(60,min(max(rtoe,1),300))
    if wait>0:
        print(f'[URLAPI] RID={rid} initial_wait={wait}s',flush=True); time.sleep(wait)
    start=time.time(); polls=0
    while True:
        status,body=search_info(rid); polls+=1
        elapsed=int(time.time()-start)
        print(f'[URLAPI] RID={rid} status={status} elapsed_resume={elapsed}s poll={polls}',flush=True)
        if status=='READY':
            xml=fetch_xml(rid)
            if '<BlastXML2' in xml or '<BlastOutput2' in xml: return xml
            if 'ThereAreHits=yes' in body or 'ThereAreHits=no' in body:
                # READY but XML formatting may lag briefly.
                time.sleep(60); continue
            raise RuntimeError('RID READY but XML2 retrieval did not contain BLAST XML2')
        if status in ('FAILED','UNKNOWN'):
            raise RuntimeError(f'RID {rid} status={status}: '+body[:1000].replace('\n',' '))
        if elapsed > max_wait:
            raise TimeoutError(f'RID {rid} exceeded resume max wait {max_wait}s; RID preserved for another resume, do not interpret biologically')
        time.sleep(60)

def extract_accession(id_text,accession):
    m=re.search(r'\|([A-Za-z0-9_]+\.\d+)\|',id_text or '')
    return m.group(1) if m else (accession or (id_text or ''))
def intervals_union_len(intervals):
    if not intervals:return 0
    xs=sorted((min(a,b),max(a,b)) for a,b in intervals); total=0; s,e=xs[0]
    for a,b in xs[1:]:
        if a<=e+1:e=max(e,b)
        else:total+=e-s+1;s,e=a,b
    return total+e-s+1

def parse_xml2_to_rows(xml_text):
    root=ET.fromstring(xml_text); rows=[]
    for search in descendants(root,'Search'):
        qid=child_text(search,'query-id'); qtitle=child_text(search,'query-title'); qlen=int(child_text(search,'query-len','0') or 0)
        qseqid=(qtitle.split()[0] if qtitle else qid)
        if qseqid.startswith('Query_') and qid and not qid.startswith('Query_'):qseqid=qid
        hits=[]
        for h in list(search):
            if local(h.tag)=='hits':hits=[x for x in list(h) if local(x.tag)=='Hit']
        for hit in hits:
            slen=int(child_text(hit,'len','0') or 0); ds=descendants(hit,'HitDescr'); d=ds[0] if ds else None
            hid=child_text(d,'id') if d is not None else ''; acc=child_text(d,'accession') if d is not None else ''
            title=child_text(d,'title') if d is not None else ''; sciname=child_text(d,'sciname') if d is not None else ''; saccver=extract_accession(hid,acc)
            hsps=[]
            for x in list(hit):
                if local(x.tag)=='hsps':hsps=[y for y in list(x) if local(y.tag)=='Hsp']
            qints=[]
            for hsp in hsps:
                qf=int(child_text(hsp,'query-from','0') or 0); qt=int(child_text(hsp,'query-to','0') or 0)
                if qf and qt:qints.append((qf,qt))
            qcovs=100.0*intervals_union_len(qints)/qlen if qlen else 0.0
            for hsp in hsps:
                alen=int(child_text(hsp,'align-len','0') or 0); ident=int(child_text(hsp,'identity','0') or 0); pident=100.0*ident/alen if alen else 0.0
                rows.append([qseqid,saccver,f'{pident:.6f}',str(alen),child_text(hsp,'query-from'),child_text(hsp,'query-to'),child_text(hsp,'hit-from'),child_text(hsp,'hit-to'),child_text(hsp,'evalue'),child_text(hsp,'bit-score'),str(qlen),str(slen),f'{qcovs:.6f}',title,sciname])
    return rows

def read_rid_file(path):
    p=Path(path)
    if not p.exists() or p.stat().st_size==0:return None,None
    d={}
    for line in p.read_text().splitlines():
        a=line.split('\t',1)
        if len(a)==2:d[a[0]]=a[1]
    return d.get('RID'),int(d.get('RTOE','60'))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--program',required=True,choices=['blastp','blastn']); ap.add_argument('--database',required=True); ap.add_argument('--query',required=True)
    ap.add_argument('--expect',required=True); ap.add_argument('--hitlist',type=int,default=200); ap.add_argument('--filter',default='T'); ap.add_argument('--megablast',action='store_true'); ap.add_argument('--no-db-override',action='store_true')
    ap.add_argument('--out',required=True); ap.add_argument('--xml-out',required=True); ap.add_argument('--rid-out',required=True); ap.add_argument('--resume-rid-file'); ap.add_argument('--resume-rid')
    ap.add_argument('--max-wait',type=int,default=64800,help='seconds to wait in this invocation; RID is preserved on timeout')
    ap.add_argument('--resubmit-if-expired',action='store_true')
    args=ap.parse_args()
    rid=args.resume_rid; rtoe=60; source='explicit'
    if not rid and args.resume_rid_file:
        rid,rtoe=read_rid_file(args.resume_rid_file); source='rid-file'
    if rid:
        status,body=search_info(rid)
        print(f'[URLAPI] resume_candidate RID={rid} source={source} status_now={status}',flush=True)
        if status in ('UNKNOWN','FAILED'):
            if not args.resubmit_if_expired: raise RuntimeError(f'existing RID {rid} is {status}; rerun with --resubmit-if-expired only if this RID has no retrievable result')
            rid=None
    if not rid:
        q=read_text(args.query); rid,rtoe=submit(args.program,args.database,q,args.expect,args.hitlist,args.filter,args.megablast,args.no_db_override)
        Path(args.rid_out).write_text(f'RID\t{rid}\nRTOE\t{rtoe}\n',encoding='utf-8')
        print(f'[URLAPI] submitted program={args.program} db={args.database} RID={rid} RTOE={rtoe}',flush=True)
    else:
        Path(args.rid_out).write_text(f'RID\t{rid}\nRTOE\t{rtoe}\nRESUMED\ttrue\n',encoding='utf-8')
    xml=wait_existing_or_new(rid,rtoe,args.max_wait)
    Path(args.xml_out).write_text(xml,encoding='utf-8'); rows=parse_xml2_to_rows(xml)
    with open(args.out,'w',encoding='utf-8') as f:
        for r in rows:f.write('\t'.join(r)+'\n')
    print(f'[URLAPI] complete RID={rid} rows={len(rows)} out={args.out}',flush=True)
if __name__=='__main__':main()
