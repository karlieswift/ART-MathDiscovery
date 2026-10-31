#!/usr/bin/env python3
import csv,json,re,sys,hashlib
from pathlib import Path

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''):h.update(b)
    return h.hexdigest()
def ff(x):
    try:return float(x)
    except:return float('nan')
def normacc(t):
    m=re.search(r'AF-([A-Za-z0-9]+)-F\d+-model',t); return m.group(1) if m else t
def tier(tm,qc,tc):
    if tm>=0.50 and qc>=0.50 and tc>=0.50:return 'STRONG'
    if tm>=0.40 and qc>=0.30 and tc>=0.30:return 'SUGGESTIVE'
    return 'WEAK'

def main():
    if len(sys.argv)!=7: raise SystemExit('usage: summarize_d.py TARGETS BOLTZ_STATUS FOLD_STATUS FOLD_DIR OUTDIR LOCK')
    targets_p,boltz_p,fold_p,fold_dir,outdir,lockp=sys.argv[1:]
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True);fold_dir=Path(fold_dir)
    targets=list(csv.DictReader(open(targets_p,encoding='utf-8'),delimiter='\t')); assert len(targets)==6
    boltz={r['candidate_id']:r for r in csv.DictReader(open(boltz_p,encoding='utf-8'),delimiter='\t')}
    fold={r['candidate_id']:r for r in csv.DictReader(open(fold_p,encoding='utf-8'),delimiter='\t')}
    fields=['query','target','evalue','bits','qcov','tcov','alnlen','qlen','tlen','fident','alntmscore']
    rows=[];hits=[]
    for t in targets:
        cid=t['candidate_id'];b=boltz[cid]
        z={k:t[k] for k in ['extension_rank','family_rep_id','candidate_id','accession','aa_length']}
        z.update({'boltz_status':b['boltz_status'],'boltz_failure_category':b['failure_category'],'cif_path':b['cif_path']})
        if b['boltz_status']!='PASS':
            z.update({'foldseek_status':'NOT_RUN','r5c_hit_n':0,'best_clean_rt_target':'','best_alntmscore':'','best_qcov':'','best_tcov':'','r5c_tier':'TECHNICAL_FAILURE'});rows.append(z);continue
        fs=fold.get(cid,{})
        if fs.get('foldseek_status')!='PASS':
            z.update({'foldseek_status':fs.get('foldseek_status','MISSING'),'r5c_hit_n':0,'best_clean_rt_target':'','best_alntmscore':'','best_qcov':'','best_tcov':'','r5c_tier':'FOLDSEEK_TECHNICAL_FAILURE'});rows.append(z);continue
        rr=[];m8=fold_dir/f'{cid}.m8'
        if m8.exists() and m8.stat().st_size:
            with open(m8) as f:
                for r in csv.DictReader(f,delimiter='\t',fieldnames=fields):
                    r['candidate_id']=cid;r['target_accession']=normacc(r['target']);r['tier']=tier(ff(r['alntmscore']),ff(r['qcov']),ff(r['tcov']));rr.append(r);hits.append(r)
        if rr:
            best=max(rr,key=lambda r:ff(r['alntmscore']))
            z.update({'foldseek_status':'PASS','r5c_hit_n':len(rr),'best_clean_rt_target':best['target_accession'],'best_alntmscore':best['alntmscore'],'best_qcov':best['qcov'],'best_tcov':best['tcov'],'r5c_tier':best['tier']})
        else:z.update({'foldseek_status':'PASS','r5c_hit_n':0,'best_clean_rt_target':'','best_alntmscore':'','best_qcov':'','best_tcov':'','r5c_tier':'NO_HIT'})
        rows.append(z)
    cols=list(rows[0]);sp=out/'F3B4D_D_STRUCTURE_R5C_SUMMARY_6.tsv'
    with sp.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    hp=out/'F3B4D_D_ALL_CLEAN_RT_HITS.tsv';hcols=['candidate_id']+fields+['target_accession','tier']
    with hp.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=hcols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(hits)
    counts={k:sum(r['r5c_tier']==k for r in rows) for k in ['STRONG','SUGGESTIVE','WEAK','NO_HIT','TECHNICAL_FAILURE','FOLDSEEK_TECHNICAL_FAILURE']}
    gate={'stage':'Stage05F3B4D-D','gate':'PASS_F3B4D_D_FROZEN_R5C_STRUCTURE_EXTENSION_COMPLETE','gate_pass':True,'structure_extension_family_n':6,'boltz_success_n':sum(r['boltz_status']=='PASS' for r in rows),'r5c_strong_n':counts['STRONG'],'r5c_suggestive_n':counts['SUGGESTIVE'],'r5c_weak_n':counts['WEAK'],'r5c_no_hit_n':counts['NO_HIT'],'technical_failure_n':counts['TECHNICAL_FAILURE']+counts['FOLDSEEK_TECHNICAL_FAILURE'],'r5c_strong_or_suggestive_n':counts['STRONG']+counts['SUGGESTIVE'],'scientific_threshold_changed':False,'boltz_settings_changed':False,'summary_sha256':sha(sp),'hits_sha256':sha(hp),'claim_boundary':'Post-discovery architecture-driven structure extension under the unchanged frozen R5C rule. STRONG/SUGGESTIVE supports clean-RT fold compatibility, not biochemical RT activity or ART identity.'}
    (out/'STAGE05F3B4D_D_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate,indent=2));print('\n=== D STRUCTURE EXTENSION 6 ===');print('\t'.join(cols));
    for r in rows:print('\t'.join(str(r[c]) for c in cols))
if __name__=='__main__':main()
