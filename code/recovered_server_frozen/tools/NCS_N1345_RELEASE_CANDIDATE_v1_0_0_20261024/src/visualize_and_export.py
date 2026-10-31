#!/usr/bin/env python3
"""N3/N4 diagnostic figures and workbook. No frozen manuscript edits."""
import argparse,csv,collections,statistics,json
from pathlib import Path
from openpyxl import Workbook
from openpyxl.styles import PatternFill,Font,Alignment,Border,Side
from openpyxl.utils import get_column_letter

def read(p):
 with open(p,encoding='utf-8') as f:return list(csv.DictReader(f,delimiter='\t'))
def sheets(root,out,inp):
 files=[('N1_Admitted',inp/'N1/N1_ADMISSION.tsv'),('N1_Sources',root/'data/N1_SOURCE_TRIAGE_20261024.tsv'),('N3_Summary',inp/'N3/N3_ROBUSTNESS_SUMMARY.tsv'),('N3_Cases',inp/'N3/N3_ROBUSTNESS_CASES.tsv'),('N3_Top10',inp/'N3/N3_TOP10_MATCHING.tsv'),('N4_Groups',inp/'N4/N4_GROUP_SUMMARY.tsv'),('N4_Queries',inp/'N4/N4_QUERY_NORMALIZATION_AUDIT_69.tsv'),('N4_Strict24',inp/'N4/N4_STRICT_24_PRIORITIZED.tsv')]
 wb=Workbook();ws=wb.active;ws.title='START';ws.append(['Evidence panel','NCS N1 N3 N4 N5; version 20261024']);ws.append(['N1','No newly approved independent ART primary truth; HOLD']);ws.append(['N3','Post-freeze, retrospective q-threshold and rank diagnostics, 6 genomes']);ws.append(['N4','67 available structures x 32 references; 2 historical OOM; no tier changes']);ws.append(['N5','Source manifests passed; external release and licenses require author approval']);ws.append(['Methods','Run scripts are included in the science-only ZIP']);ws.append(['Publication use','Do NOT mislabel sensitivity or structural similarity as functionally validated'])
 for title,path in files:
  ws=wb.create_sheet(title);rows=read(path)
  if rows: headers=list(rows[0]);ws.append(headers)
  for row in rows:
   record=[]
   for col in headers:
    val=row[col]
    if val and col not in {'q_threshold','source_id','accession_version','publication_doi','source_url','genome_sha256','query_id','membership_id','query_length_aa','case_id','archived_tier_unchanged','best_ref_by_tm_min','best_ref_by_tm_query','best_ref_by_tm_reference','best_ref_by_tm_mean','cohort','scope','descriptive_flags','expert_review','biochemical_function','source','reason_not_admitted','evidence_type','data_type','algorithm'}:
     try: val=int(val) if val.isdigit() else float(val)
     except ValueError: pass
    record.append(val)
   ws.append(record)
  ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions
 for ws in wb:
  ws.sheet_view.showGridLines=False
  ws.row_dimensions[1].height=36
  for cell in ws[1]:
   cell.fill=PatternFill('solid',fgColor='143E5A');cell.font=Font(name='Aptos',size=10,color='FFFFFF',bold=True);cell.alignment=Alignment(wrap_text=True,vertical='center')
  for c in range(1,ws.max_column+1):
   values=[str(ws.cell(row=r,column=c).value or '') for r in range(1,min(ws.max_row,80)+1)]
   maxlen=min(57,max([len(v) for v in values]+[12]));ws.column_dimensions[get_column_letter(c)].width=max(16,min(52,int(maxlen*.85)+3))
  for row in ws.iter_rows(min_row=2):
   for c in row:
    c.font=Font(name='Aptos',size=10,color=('087443' if isinstance(c.value,(float,int)) else '384A5A'))
    if isinstance(c.value,float):c.number_format='0.0000'
    if c.value and isinstance(c.value,str) and any(t in c.value for t in ['HOLD','PENDING','not independent','NO_NEW','NOT_TESTED']): c.font=Font(name='Aptos',size=10,color='B66C18')
    if c.row%2==0:c.fill=PatternFill('solid',fgColor='F2F6F8')
 out.parent.mkdir(parents=True,exist_ok=True);wb.save(out)
 return {'sheets':wb.sheetnames,'file':str(out)}

def figures(root,out,inp):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 out.mkdir(parents=True,exist_ok=True)
 rows=read(inp/'N3/N3_ROBUSTNESS_CASES.tsv')
 axx=sorted({float(r['q_threshold']) for r in rows})
 vals={float(x['q_threshold']):x for x in read(inp/'N3/N3_ROBUSTNESS_SUMMARY.tsv')}
 fig,ax=plt.subplots(figsize=(7.5,4.5))
 ax.plot(axx,[float(vals[q]['total_regions']) for q in axx],marker='o',linewidth=2)
 ax.axvline(.005,ls='--',linewidth=1)
 ax.set_xlabel('Empirical q threshold (not genome-wide FDR)');ax.set_ylabel('Candidate regions (6 previously seen genomes)')
 ax.set_title('N3. Candidate burden under post-freeze q perturbation')
 fig.tight_layout()
 for ext in ['png','pdf']:fig.savefig(out/('N3_Total_Regions_20261024.'+ext),dpi=240)
 plt.close(fig)
 fig,ax=plt.subplots(figsize=(8.5,5.4))
 for case in sorted({r['case_id'] for r in rows}):
  sub=sorted([r for r in rows if r['case_id']==case],key=lambda z:float(z['q_threshold']))
  ax.plot([float(z['q_threshold']) for z in sub],[float(z['interval_jaccard']) for z in sub],marker='o',label=case)
 ax.set_ylim(0,1.07);ax.set_xlabel('Empirical q threshold');ax.set_ylabel('Interval Jaccard relative to q=0.005')
 ax.set_title('N3. Footprint stability is threshold dependent');ax.legend(ncol=3,fontsize=8,loc='lower center')
 fig.tight_layout()
 for ext in ['png','pdf']:fig.savefig(out/('N3_Interval_Stability_20261024.'+ext),dpi=240)
 plt.close(fig)
 r4=read(inp/'N4/N4_QUERY_NORMALIZATION_AUDIT_69.tsv')
 fig,ax=plt.subplots(figsize=(8.4,5))
 marker={'historical_6p12Gb':'o','strict_mirror_54Gb':'s','architecture_complete_54Gb':'^'}
 for cohort in marker:
  xx=[z for z in r4 if z['cohort']==cohort and z['pairs_count']=='32']
  ax.scatter([float(z['tm_min_of_best_pair']) for z in xx],[float(z['best_query_fraction']) for z in xx],marker=marker[cohort],label=cohort,s=42,alpha=.8)
 ax.axvline(.3,ls='--',linewidth=1);ax.axhline(.5,ls='--',linewidth=1)
 ax.set_xlabel('Best-reference min-normalized US-align TM-score');ax.set_ylabel('Fraction of query residues aligned')
 ax.set_title('N4. Full-length structural alignment: descriptive QC');ax.legend(fontsize=8)
 fig.tight_layout()
 for ext in ['png','pdf']:fig.savefig(out/('N4_TM_vs_Coverage_20261024.'+ext),dpi=240)
 plt.close(fig)
 return {'figures':6}
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--root',required=True);a.add_argument('--out',required=True);a.add_argument('--results');x=a.parse_args();p=Path(x.out);p.mkdir(parents=True,exist_ok=True);inp=Path(x.results) if x.results else Path(x.root)/'results';res={'excel':sheets(Path(x.root),p/'NCS_N1N3N4_EVIDENCE_20261024.xlsx',inp),'plots':figures(Path(x.root),p,inp)};print(json.dumps(res,indent=2))
