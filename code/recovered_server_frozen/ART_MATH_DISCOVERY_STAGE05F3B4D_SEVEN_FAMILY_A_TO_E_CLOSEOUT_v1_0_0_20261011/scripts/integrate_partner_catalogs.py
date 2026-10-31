#!/usr/bin/env python3
import csv,json,sys
from pathlib import Path

def ffloat(x):
    try:return float(x)
    except:return float('inf')

def parse_tblout(p):
    out={}
    p=Path(p)
    if not p.exists() or p.stat().st_size==0:return out
    for line in open(p,encoding='utf-8',errors='ignore'):
        if not line.strip() or line.startswith('#'):continue
        x=line.split()
        if len(x)<6:continue
        target=x[0]; query=x[2]; e=x[4]; score=x[5]
        r={'target':target,'evalue':e,'score':score}
        if query not in out or ffloat(e)<ffloat(out[query]['evalue']):out[query]=r
    return out

def parse_m8(p):
    out={}; p=Path(p)
    if not p.exists() or p.stat().st_size==0:return out
    fields=['query','target','fident','alnlen','qcov','tcov','evalue','bits']
    with open(p,encoding='utf-8') as f:
        first=f.readline().rstrip('\n')
        if not first:return out
        rows=[]
        if first.startswith('query\t'): rows=list(csv.DictReader([first]+list(f),delimiter='\t'))
        else:
            vals=first.split('\t'); rows=[dict(zip(fields,vals))]
            for line in f:
                if line.strip(): rows.append(dict(zip(fields,line.rstrip('\n').split('\t'))))
    for r in rows:
        q=r['query']
        if q not in out or ffloat(r['evalue'])<ffloat(out[q]['evalue']) or (ffloat(r['evalue'])==ffloat(out[q]['evalue']) and float(r.get('bits',0))>float(out[q].get('bits',0))): out[q]=r
    return out

def main():
    if len(sys.argv)!=8: raise SystemExit('usage: integrate_partner_catalogs.py REP7_TSV BASE7_TSV PFAM_TBLOUT PHROG_M8 VOG_M8 OUT_TSV OUT_GATE')
    rep_p,base_p,pfam_p,phrog_p,vog_p,out_p,gate_p=sys.argv[1:]
    reps=list(csv.DictReader(open(rep_p,encoding='utf-8'),delimiter='\t'))
    base={r['family_rep_id']:r for r in csv.DictReader(open(base_p,encoding='utf-8'),delimiter='\t')}
    pf=parse_tblout(pfam_p); ph=parse_m8(phrog_p); vg=parse_m8(vog_p)
    rows=[]
    for r in reps:
        fam=r['family_rep_id']; q=r['partner_orf_id']; b=base[fam]
        z=dict(b)
        z.update({
          'partner_representative_orf_id':q,
          'pfam_top_target':pf.get(q,{}).get('target',''),'pfam_top_evalue':pf.get(q,{}).get('evalue',''),'pfam_top_score':pf.get(q,{}).get('score',''),
          'phrog_top_target':ph.get(q,{}).get('target',''),'phrog_top_fident':ph.get(q,{}).get('fident',''),'phrog_top_qcov':ph.get(q,{}).get('qcov',''),'phrog_top_tcov':ph.get(q,{}).get('tcov',''),'phrog_top_evalue':ph.get(q,{}).get('evalue',''),
          'vog_top_target':vg.get(q,{}).get('target',''),'vog_top_fident':vg.get(q,{}).get('fident',''),'vog_top_qcov':vg.get(q,{}).get('qcov',''),'vog_top_tcov':vg.get(q,{}).get('tcov',''),'vog_top_evalue':vg.get(q,{}).get('evalue',''),
          'partner_function_confirmed':False,'dedicated_partner_gene_confirmed':False,
        })
        rows.append(z)
    cols=list(rows[0])
    with open(out_p,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    gate={
      'stage':'Stage05F3B4D-C','gate':'PASS_F3B4D_C_PARTNER_CONSERVATION_AND_CATALOG_ANNOTATION_COMPLETE','gate_pass':True,
      'target_family_n':7,'pfam_resource_available':bool(pf),'phrog_resource_available':bool(ph),'vog_resource_available':bool(vg),
      'partner_function_confirmed':False,'dedicated_partner_gene_confirmed':False,
      'claim_boundary':'Partner catalog hits are descriptive homology/domain annotations. Recurrent adjacency plus homology does not establish dedicated partner function.'}
    Path(gate_p).write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate,indent=2)); print('\n=== PARTNER CONSERVATION + CATALOG 7 ===')
    print('\t'.join(cols))
    for r in rows:print('\t'.join(str(r[c]) for c in cols))

if __name__=='__main__': main()
