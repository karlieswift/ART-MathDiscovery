#!/usr/bin/env python3
"""Obtain official author document and locate dataset links/identifier mentions.
An accession mention is not a confirmed catalogue member. Missing mentions never
establish novelty. Source-PDF raw bytes are preserved for later manual review.
"""
import argparse,json,re,subprocess,shutil
from pathlib import Path
from urllib.parse import urlparse
from common import *

def extract_pdf(pdf,out,accessions):
 import fitz
 pdf=Path(pdf);out=Path(out);out.mkdir(parents=True,exist_ok=True)
 if not pdf.read_bytes()[:5]==b'%PDF-':raise ValueError('not a PDF file')
 doc=fitz.open(pdf)
 if len(doc)<1 or doc.is_encrypted:raise ValueError('empty or encrypted PDF')
 links=[];mentions=[];availability=[];embedded=[];pages=[]
 for i,page in enumerate(doc):
  text=page.get_text();pages.append(text)
  for link in page.get_links():
   uri=link.get('uri','')
   if uri:links.append({'page_1based':i+1,'url':uri,'host':urlparse(uri).netloc,'source':'PDF_URI_ANNOTATION','author_link_is_not_automatically_author_owned_dataset':True})
  for uri in re.findall(r'https?://[^\s<>]+',text):
   links.append({'page_1based':i+1,'url':uri.rstrip('.,;)'),'host':urlparse(uri).netloc,'source':'PDF_TEXT_URL_REQUIRES_REVIEW','author_link_is_not_automatically_author_owned_dataset':True})
  for n,ln in enumerate(text.splitlines(),1):
   if re.search(r'data\s+(and\s+code\s+)?availab|code\s+availab|supplementary\s+(data|table|file)',ln,re.I):availability.append({'page_1based':i+1,'text_line_1based':n,'category':'AVAILABILITY_OR_SUPPLEMENT_LOCATION_REQUIRES_READ'})
  for acc in accessions:
   for kind,term in [('EXACT_VERSION',acc),('ACCESSION_ROOT_ONLY',acc.split('.')[0])]:
    # Root matches are retained even when an exact match is also present.
    count=len(re.findall(r'(?<![A-Za-z0-9_])'+re.escape(term)+r'(?![A-Za-z0-9_])',text))
    if count:mentions.append({'accession':acc,'page_1based':i+1,'match_kind':kind,'occurrences':count,'catalogue_membership_confirmed':False})
 for name in doc.embfile_names():
  b=doc.embfile_get(name);dest=out/'embedded_files'/(sha_bytes(b)[:16]+'_'+re.sub(r'[^A-Za-z0-9._-]','_',name));dest.parent.mkdir(exist_ok=True);dest.write_bytes(b)
  embedded.append({'author_embedded_name':name,'local_file':str(dest.relative_to(out)),'size_bytes':len(b),'sha256':sha_bytes(b),'full_family_catalogue_status':'NOT_YET_PARSED_OR_VALIDATED'})
 write_tsv(out/'AUTHOR_PDF_LINKS.tsv',links,['page_1based','url','host','source','author_link_is_not_automatically_author_owned_dataset'])
 write_tsv(out/'EXACT_ACCESSION_MENTIONS.tsv',mentions,['accession','page_1based','match_kind','occurrences','catalogue_membership_confirmed'])
 write_tsv(out/'AVAILABILITY_SECTION_LOCATORS.tsv',availability,['page_1based','text_line_1based','category'])
 write_tsv(out/'AUTHOR_EMBEDDED_FILES.tsv',embedded,['author_embedded_name','local_file','size_bytes','sha256','full_family_catalogue_status'])
 status=[{'accession':acc,'exact_version_mentioned':any(r['accession']==acc and r['match_kind']=='EXACT_VERSION' for r in mentions),'root_mentioned':any(r['accession']==acc for r in mentions),'interpretation':'TEXT_MENTION_REQUIRES_CONTEXT_NOT_MEMBERSHIP' if any(r['accession']==acc for r in mentions) else 'NOT_FOUND_IN_RETRIEVED_PDF_TEXT_NOT_NOVELTY','novelty_confirmed':False} for acc in accessions]
 write_tsv(out/'TEN_LOCUS_AUTHOR_DOCUMENT_STATUS.tsv',status)
 # The document text is useful locally; do not present it as author catalogue.
 (out.parent.parent/'work/author').mkdir(parents=True,exist_ok=True)
 (out.parent.parent/'work/author/technical_report_text_local.txt').write_text('\n\f\n'.join(pages))
 return {'page_count':len(doc),'PDF_sha256':sha_file(pdf),'n_links':len(links),'n_embedded_files':len(embedded),'n_exact_accessions_mentioned':sum(r['exact_version_mentioned'] for r in status),'full_family_catalogue_acquired_and_validated':False}

def run(root,local_pdf=None):
 root=Path(root);out=root/'outputs/author_source';out.mkdir(parents=True,exist_ok=True);work=root/'work/author';work.mkdir(parents=True,exist_ok=True)
 pol=read_json(root/'config/POLICY.json')['author_source'];accs=[r['sequence_id'] for r in read_tsv(root/'precomputed/TEN_LOCUS_UNAMBIGUOUS_LEDGER.tsv')]
 target=work/'anthropic_ART_technical_report.pdf';error=''
 try:
  if local_pdf:
   src=Path(local_pdf).resolve()
   if target.exists() and sha_file(src)!=sha_file(target):raise ValueError('existing author snapshot differs; preserve both, do not overwrite')
   if not target.exists():shutil.copyfile(src,target)
  elif not target.exists():
   tmp=target.with_suffix('.pdf.part')
   wget=shutil.which('wget');curl=shutil.which('curl')
   if wget:cmd=[wget,'--tries=2','--timeout=45','--waitretry=3','-O',str(tmp),pol['technical_report_url']]
   elif curl:cmd=[curl,'--fail','--location','--retry','1','--connect-timeout','30','--max-time','240','-o',str(tmp),pol['technical_report_url']]
   else:raise RuntimeError('wget/curl missing: supply --pdf local_author_report.pdf')
   cp=subprocess.run(cmd,capture_output=True,text=True,timeout=300);(out/'download.log').write_text(cp.stdout+'\n'+cp.stderr)
   if cp.returncode:raise RuntimeError('official PDF transfer failed: rc='+str(cp.returncode))
   if tmp.stat().st_size<1000 or tmp.read_bytes()[:5]!=b'%PDF-':raise ValueError('download not a valid PDF payload')
   tmp.rename(target)
  info=extract_pdf(target,out,accs)
  old=out/'AUTHOR_SOURCE_GATE.json'
  if old.exists():
   prev=read_json(old)
   if prev.get('PDF_sha256') and prev['PDF_sha256']!=info['PDF_sha256']:raise ValueError('PDF snapshot changed since previous source check')
  gate={'gate':'PASS_AUTHOR_DOCUMENT_RETRIEVED_CATALOGUE_MEMBERSHIP_STILL_UNRESOLVED','document_retrieval_pass':True,'source_url':pol['technical_report_url'],'source_supplied_locally':bool(local_pdf),**info,'catalogue_novelty_gate_pass':False,'accession_absence_is_not_novelty':True,'new_ART_confirmed':0}
 except Exception as e:
  gate={'gate':'HOLD_AUTHOR_DOCUMENT_OR_CATALOGUE_NOT_RETRIEVED','document_retrieval_pass':False,'error':type(e).__name__+': '+str(e),'source_url':pol['technical_report_url'],'catalogue_novelty_gate_pass':False,'new_ART_confirmed':0}
 write_json(out/'AUTHOR_SOURCE_GATE.json',gate);print(json.dumps(gate,indent=2))
 return gate
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--root',required=True);a.add_argument('--pdf');x=a.parse_args();run(x.root,x.pdf)
