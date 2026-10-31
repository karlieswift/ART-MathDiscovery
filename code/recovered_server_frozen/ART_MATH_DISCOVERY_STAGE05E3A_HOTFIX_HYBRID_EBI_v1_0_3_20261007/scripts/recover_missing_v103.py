#!/usr/bin/env python3
import argparse, csv, hashlib, json, os, re, subprocess, tempfile, time
from pathlib import Path
from urllib.parse import quote
from Bio import SeqIO

VERSION = "v1.0.3"
REFSEQ_PREFIX = re.compile(r'^(?:NC_|NG_|NM_|NR_|NW_|NZ_|NT_|XM_|XR_|AC_)')
ENA_Q = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ena_sequence&id={id}&format=embl&style=raw"
ENA_R = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch/ena_sequence/{id}/embl?style=raw"
REF_Q = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqn&id={id}&format=refseq&style=raw"
REF_R = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch/refseqn/{id}/refseq?style=raw"

def sha256_bytes(b): return hashlib.sha256(b).hexdigest()
def first_nonempty(text):
    for line in text.splitlines():
        if line.strip(): return line.strip()[:300]
    return ""

def curl_fetch(url, out_path, hdr_path, timeout=120):
    cmd=["curl","-fsSL","--retry","5","--retry-all-errors","--retry-delay","3",
         "--connect-timeout","20","--max-time",str(timeout),
         "-A","ART-MathDiscovery-Stage05E3A-EBI-v1.0.3",
         "-D",str(hdr_path),"-o",str(out_path),url]
    p=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    return p.returncode,p.stderr.strip()

def parse_status_and_type(hdr):
    status=ctype=""
    if not hdr.exists(): return status,ctype
    for ln in hdr.read_text(errors="replace").splitlines():
        if ln.startswith("HTTP/"):
            sp=ln.split(); status=sp[1] if len(sp)>=2 else status
        elif ln.lower().startswith("content-type:"):
            ctype=ln.split(":",1)[1].strip()
    return status,ctype

def accession_ok(rec, requested):
    rb=requested.split('.')[0]
    vals=[str(rec.id),str(rec.name)]
    vals += [str(x) for x in rec.annotations.get('accessions',[]) or []]
    vals += [str(rec.annotations.get('sequence_version',''))]
    return any(v.split('.')[0]==rb for v in vals if v)

def validate_refseq(path, requested):
    data=path.read_bytes() if path.exists() else b""
    text=data.decode('utf-8','replace')
    if not text.lstrip().startswith('LOCUS'):
        return None,'does_not_start_LOCUS',data
    try: recs=list(SeqIO.parse(str(path),'genbank'))
    except Exception as e: return None,f'biopython_genbank_parse_error:{type(e).__name__}:{e}',data
    if len(recs)!=1: return None,f'genbank_record_count_{len(recs)}',data
    rec=recs[0]
    if len(rec.seq)<=0: return None,'empty_sequence',data
    if not accession_ok(rec,requested): return None,f'accession_mismatch_requested_{requested}_record_{rec.id}',data
    return rec,'OK',data

def validate_embl(path, requested):
    data=path.read_bytes() if path.exists() else b""
    text=data.decode('utf-8','replace')
    if not text.lstrip().startswith('ID'):
        return None,'does_not_start_EMBL_ID',data
    try: recs=list(SeqIO.parse(str(path),'embl'))
    except Exception as e: return None,f'biopython_embl_parse_error:{type(e).__name__}:{e}',data
    if len(recs)!=1: return None,f'embl_record_count_{len(recs)}',data
    rec=recs[0]
    if len(rec.seq)<=0: return None,'empty_sequence',data
    if not accession_ok(rec,requested): return None,f'accession_mismatch_requested_{requested}_record_{rec.id}',data
    return rec,'OK',data

def write_genbank(rec,out):
    if 'molecule_type' not in rec.annotations: rec.annotations['molecule_type']='DNA'
    tmp=out.with_suffix(out.suffix+'.tmp')
    if tmp.exists(): tmp.unlink()
    SeqIO.write(rec,str(tmp),'genbank')
    rs=list(SeqIO.parse(str(tmp),'genbank'))
    if len(rs)!=1 or len(rs[0].seq)<=0: raise RuntimeError('generated GenBank validation failed')
    os.replace(tmp,out)

def direct_copy_genbank(src,out):
    tmp=out.with_suffix(out.suffix+'.tmp')
    tmp.write_bytes(src.read_bytes())
    rs=list(SeqIO.parse(str(tmp),'genbank'))
    if len(rs)!=1 or len(rs[0].seq)<=0: raise RuntimeError('copied RefSeq validation failed')
    os.replace(tmp,out)

def routes(sid):
    base=sid.split('.')[0]
    if REFSEQ_PREFIX.match(sid):
        ids=[sid] + ([] if base==sid else [base])
        out=[]
        for q in ids:
            out += [
              (f'EBI_REFSEQN_QUERY:{q}','REFSEQ',REF_Q.format(id=quote(q,safe='._-'))),
              (f'EBI_REFSEQN_REST:{q}','REFSEQ',REF_R.format(id=quote(q,safe='._-'))),
            ]
        return out
    ids=[sid] + ([] if base==sid else [base])
    out=[]
    for q in ids:
        out += [
          (f'EBI_ENA_QUERY:{q}','EMBL',ENA_Q.format(id=quote(q,safe='._-'))),
          (f'EBI_ENA_REST:{q}','EMBL',ENA_R.format(id=quote(q,safe='._-'))),
        ]
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--ids',required=True)
    ap.add_argument('--outdir',required=True)
    ap.add_argument('--manifest',required=True)
    ap.add_argument('--diagnostics',required=True)
    ap.add_argument('--sleep',type=float,default=1.0)
    a=ap.parse_args()
    ids=[x.strip() for x in Path(a.ids).read_text().splitlines() if x.strip()]
    outdir=Path(a.outdir); outdir.mkdir(parents=True,exist_ok=True)
    manifest=Path(a.manifest); manifest.parent.mkdir(parents=True,exist_ok=True)
    diagnostics=Path(a.diagnostics); diagnostics.parent.mkdir(parents=True,exist_ok=True)
    dfields=['sequence_id','route','payload_format','http_status','content_type','bytes','sha256','valid_record','validation_or_error','first_nonempty_line','curl_error','elapsed_s']
    mfields=['sequence_id','status','bytes','sha256','selected_route','provider','source_format','output_format']
    drows=[]; mrows=[]; failed=[]
    for i,sid in enumerate(ids,1):
        selected=None; last=''
        for route,fmt,url in routes(sid):
            with tempfile.TemporaryDirectory(prefix='e3a_v103_') as td0:
                td=Path(td0); body=td/'body.txt'; hdr=td/'headers.txt'
                t0=time.time(); rc,cerr=curl_fetch(url,body,hdr); elapsed=time.time()-t0
                status,ctype=parse_status_and_type(hdr); data=body.read_bytes() if body.exists() else b''
                rec=None
                if rc==0:
                    if fmt=='REFSEQ': rec,verr,data=validate_refseq(body,sid)
                    else: rec,verr,data=validate_embl(body,sid)
                else: verr=f'curl_rc_{rc}'
                drows.append({'sequence_id':sid,'route':route,'payload_format':fmt,'http_status':status,'content_type':ctype,'bytes':len(data),
                              'sha256':sha256_bytes(data) if data else '','valid_record':int(rec is not None),'validation_or_error':verr,
                              'first_nonempty_line':first_nonempty(data.decode('utf-8','replace')) if data else '',
                              'curl_error':cerr[:500],'elapsed_s':round(elapsed,3)})
                if rec is not None:
                    out=outdir/f'{sid}.gbff'
                    if fmt=='REFSEQ': direct_copy_genbank(body,out)
                    else: write_genbank(rec,out)
                    ob=out.read_bytes(); selected=(route,fmt,len(ob),sha256_bytes(ob)); break
                last=f'{route}:{verr}:{cerr[:180]}'
            time.sleep(max(0.0,a.sleep))
        if selected:
            route,fmt,n,h=selected
            provider='EMBL-EBI_DBFETCH_REFSEQN' if fmt=='REFSEQ' else 'EMBL-EBI_ENA'
            source='RefSeq_GenBank_flatfile' if fmt=='REFSEQ' else 'EMBL'
            output='GenBank_direct_validated' if fmt=='REFSEQ' else 'GenBank_converted_by_Biopython'
            mrows.append({'sequence_id':sid,'status':'VALID','bytes':n,'sha256':h,'selected_route':route,'provider':provider,'source_format':source,'output_format':output})
            print(f'[{i}/{len(ids)}] {sid} FETCHED route={route}',flush=True)
        else:
            failed.append((sid,last))
            mrows.append({'sequence_id':sid,'status':'FAILED','bytes':0,'sha256':'','selected_route':'NONE','provider':'EMBL-EBI','source_format':'AUTO','output_format':'GenBank'})
            print(f'[{i}/{len(ids)}] {sid} FAILED route=NONE',flush=True)
    with diagnostics.open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=dfields,delimiter='\t'); w.writeheader(); w.writerows(drows)
    with manifest.open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=mfields,delimiter='\t'); w.writeheader(); w.writerows(mrows)
    summary={'requested':len(ids),'valid':len(ids)-len(failed),'failed':len(failed),'failed_ids':[x[0] for x in failed],
             'transport_hotfix':'v1.0.3 hybrid EBI recovery: RefSeq nucleotide dbfetch for RefSeq accessions; ENA dbfetch retry for INSDC accession; preserve existing v1.0.2 validated records',
             'scientific_settings_changed':False}
    print(json.dumps(summary,indent=2))
    if failed: raise SystemExit(f'v1.0.3 recovery incomplete: {failed[:5]}')
if __name__=='__main__': main()
