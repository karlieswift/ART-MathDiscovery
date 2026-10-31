#!/usr/bin/env python3
import argparse,hashlib,json,csv
from pathlib import Path

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for c in iter(lambda:f.read(1<<20),b''):h.update(c)
 return h.hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--out',required=True);a=p.parse_args()
 root=Path(a.root); manifest=root/'SOURCE_INPUTS_SHA256.tsv';errors=[];n=0
 with manifest.open() as f:
  for v in csv.DictReader(f,delimiter='\t'):
   path=root/v['path']; n+=1
   if not path.is_file() or sha(path)!=v['sha256']:errors.append(v['path'])
 status={'status':'PASS_RELEASE_SOURCE_MANIFEST' if not errors else 'FAIL_RELEASE_SOURCE_MANIFEST','source_inputs_checked':n,'errors':errors,'public_github_repository_created':False,'zenodo_doi_created':False,'license_approved':False}
 Path(a.out).write_text(json.dumps(status,indent=2)+'\n');print(json.dumps(status,indent=2))
 if errors:raise SystemExit(2)
if __name__=='__main__':main()
