#!/usr/bin/env python3
"""Retrieve only the exact previously frozen myRT HMM, with SHA256 verification."""
import argparse,hashlib,json,os,sys,urllib.request
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent))
from j2_pipeline import RVT_SHA,sha,write_json,now
URL='https://raw.githubusercontent.com/mgtools/myRT/main/Models/HMM/RVT-All.hmm'
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args();dest=Path(a.output).resolve()
 if dest.exists():
  if sha(dest)!=RVT_SHA:raise SystemExit('HOLD existing file has wrong SHA256; not overwritten')
  print(json.dumps({'status':'EXACT_FROZEN_HMM_ALREADY_PRESENT','path':str(dest),'sha256':RVT_SHA}));return
 dest.parent.mkdir(parents=True,exist_ok=True);tmp=dest.with_name(dest.name+'.download-part')
 if tmp.exists():raise SystemExit('HOLD incomplete download path exists; review it before retry')
 req=urllib.request.Request(URL,headers={'User-Agent':'ART-MathDiscovery-J2-reference-recovery/1.0'})
 try:
  with urllib.request.urlopen(req,timeout=60) as r,tmp.open('xb') as f:
   while True:
    block=r.read(1024*1024)
    if not block:break
    f.write(block)
  observed=sha(tmp)
  if observed!=RVT_SHA:raise ValueError('download does not match frozen HMM: '+observed)
  os.replace(tmp,dest)
 except Exception:
  tmp.unlink(missing_ok=True);raise
 report={'status':'EXACT_FROZEN_HMM_RECOVERED','url':URL,'downloaded_utc':now(),
  'sha256':RVT_SHA,'path':str(dest),'bytes':dest.stat().st_size,
  'license':'upstream myRT terms apply; this downloader does not relicense the model'}
 write_json(dest.with_suffix(dest.suffix+'.source.json'),report);print(json.dumps(report,indent=2))
if __name__=='__main__':main()
