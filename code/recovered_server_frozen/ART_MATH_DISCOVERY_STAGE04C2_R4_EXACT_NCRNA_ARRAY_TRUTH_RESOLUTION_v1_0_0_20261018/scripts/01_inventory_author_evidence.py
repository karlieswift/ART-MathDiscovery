#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, json, os, re
from common import sha256, write_tsv, dump_json
try:
    import openpyxl
except Exception as e:
    raise SystemExit('openpyxl required: '+str(e))

TARGETS=['Alistipes-putredinis','Bacteroides-coprophilus','Bacteroides-sp-9-1-42FAA','Bacteroides-sp-2-1-16','Blautia-hansenii']
KEYWORDS=['UG27','ncrna','ncRNA','array','pseudoknot','hairpin','covariance','cmsearch','infernal']
CAND_EXT={'.cm','.sto','.stk','.afa','.aln','.fa','.fasta','.fna','.tsv','.csv','.bed','.gff','.gff3','.json','.txt','.xlsx','.xls','.gb','.gbk','.gbff'}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--ug-root',required=True); ap.add_argument('--out',required=True); ap.add_argument('--extra-root',default='')
    a=ap.parse_args(); ug=Path(a.ug_root); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    xlsx=ug/'inputs/Minerva_SupplementaryTables_preprint.xlsx'
    wb=openpyxl.load_workbook(xlsx,read_only=True,data_only=True)
    hits=[]
    for ws in wb.worksheets:
        for r,row in enumerate(ws.iter_rows(values_only=True),start=1):
            vals=['' if v is None else str(v) for v in row]
            joined=' | '.join(vals)
            low=joined.lower()
            matched=[]
            for t in TARGETS:
                if t.lower() in low: matched.append(t)
            kws=[k for k in KEYWORDS if k.lower() in low]
            if matched or kws:
                hits.append({'sheet':ws.title,'row':r,'targets':';'.join(matched),'keywords':';'.join(kws),'row_text':joined[:20000]})
    write_tsv(out/'WORKBOOK_RELEVANT_CELLS.tsv',hits,['sheet','row','targets','keywords','row_text'])

    roots=[ug/'inputs',ug/'outputs',ug/'work/stage04c2b_r3']
    if a.extra_root: roots.append(Path(a.extra_root))
    inv=[]; seen=set()
    for root in roots:
        if not root.exists(): continue
        for dp,dns,fns in os.walk(root):
            # Avoid deep/huge derived caches if extra-root is broad.
            dns[:] = [d for d in dns if d not in {'__pycache__','.git','blastdb'}]
            for fn in fns:
                p=Path(dp)/fn
                if p in seen: continue
                seen.add(p)
                low=fn.lower(); ext=p.suffix.lower()
                if ext in CAND_EXT and any(k.lower() in low for k in KEYWORDS+TARGETS):
                    try: size=p.stat().st_size; h=sha256(p) if size < 200*1024*1024 else ''
                    except Exception: size=-1; h=''
                    inv.append({'path':str(p),'size':size,'sha256':h,'extension':ext})
    write_tsv(out/'AUTHOR_EVIDENCE_FILE_INVENTORY.tsv',inv,['path','size','sha256','extension'])
    cm=[x for x in inv if x['extension']=='.cm']
    coord=[x for x in inv if x['extension'] in {'.bed','.gff','.gff3','.tsv','.csv'} and any(k in Path(x['path']).name.lower() for k in ['ug27','ncrna','array'])]
    d={'stage':'Stage04C2-R4A-evidence-inventory','workbook_relevant_rows_n':len(hits),'candidate_files_n':len(inv),'candidate_cm_files_n':len(cm),'candidate_coordinate_or_annotation_files_n':len(coord),'note':'Inventory only. Candidate files are not automatically accepted as exact truth. Exact truth resolver applies provenance rules.'}
    dump_json(out/'AUTHOR_EVIDENCE_INVENTORY_SUMMARY.json',d); print(json.dumps(d,indent=2))
if __name__=='__main__': main()
