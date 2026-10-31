#!/usr/bin/env python3
import argparse,json,csv,hashlib
from pathlib import Path

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def readj(p): return json.loads(Path(p).read_text())
def feat_overlap(f,s,e): return f.get('start1') and f['start1']<=e and f['end1']>=s

def one_locus(fs,tag):
    z=[f for f in fs if f.get('locus_tag')==tag and f['type']=='CDS']
    if len(z)!=1: raise RuntimeError(f'{tag}: expected one CDS, got {len(z)}')
    return z[0]
def rt_like(f):
    t=(f.get('product','')+' '+f.get('gene','')).lower()
    return 'reverse transcriptase' in t or 'rna-directed dna polymerase' in t

def row(control,acc,component,f,source,status='RESOLVED'):
    return dict(control_id=control,accession=acc,truth_status=status,component=component,start1=f['start1'],end1=f['end1'],start0=f['start1']-1,end0=f['end1'],strand=f['strand'],locus_tag=f.get('locus_tag',''),protein_id=f.get('protein_id',''),product=f.get('product',''),coordinate_source=source)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); ap.add_argument('--c0',required=True); ap.add_argument('--r3',required=True); a=ap.parse_args()
    R=Path(a.root); out=[]; audit={}
    # ensure R3 gate is valid and C0 is unchanged registry source
    g=readj(Path(a.r3)/'outputs/STAGE04C0_R3_FINAL_GATE.json')
    if g.get('gate')!='PASS_STAGE04C0_R3_SEQUENCE_LENGTH_CORRECTION_FIXED': raise RuntimeError('R3 gate mismatch')
    c0truth=Path(a.c0)/'truth/CONTROL_LOCUS_TRUTH_SEALED.tsv'
    if not c0truth.exists(): raise RuntimeError('missing original C0 sealed truth')
    # SEN2
    fs=readj(R/'work_SEN2_14028S_features.json')
    rt=one_locus(fs,'STM14_4641'); rca=one_locus(fs,'STM14_4640')
    if not rt_like(rt): raise RuntimeError('SEN2 STM14_4641 not annotated RT-like')
    out += [row('SEN2_14028S','CP001363.1','RT_rrtT_STM14_4641',rt,'DDBJ_CP001363.1_locus_tag+literature'),row('SEN2_14028S','CP001363.1','RcaT_STM14_4640',rca,'DDBJ_CP001363.1_locus_tag+literature')]
    env={'start1':min(rt['start1'],rca['start1']),'end1':max(rt['end1'],rca['end1']),'strand':rt['strand'],'locus_tag':'STM14_4640;STM14_4641','protein_id':'','product':'RcaT+RT coding envelope'}
    out.append(row('SEN2_14028S','CP001363.1','CODING_SYSTEM_ENVELOPE',env,'resolved coding components; ncRNA boundary not invented','PARTIAL_SYSTEM'))
    audit['SEN2_14028S']={'rt':rt,'rcaT':rca}
    # EC78: RT 311 aa and RT-like annotation; neighbor ATPase/HNH evidence within 15 kb
    fs=readj(R/'work_EC78_ECO7_features.json'); cds=[f for f in fs if f['type']=='CDS']
    cands=[f for f in cds if rt_like(f) and (f.get('aa_len') in (310,311,312) or f.get('aa_len') is None)]
    if len(cands)!=1:
        # fallback: all 311 aa CDS then require neighboring ATPase/HNH products
        cands=[f for f in cds if f.get('aa_len')==311]
    good=[]
    for f in cands:
        neigh=[x for x in cds if x['start1']<=f['end1']+15000 and x['end1']>=f['start1']-15000]
        txt=' '.join((x.get('product','')+' '+x.get('gene','')).lower() for x in neigh)
        if ('atpase' in txt or 'atp-binding' in txt or 'abc' in txt) and ('hnh' in txt or 'nuclease' in txt or 'endonuclease' in txt): good.append((f,neigh))
    if len(good)!=1: raise RuntimeError(f'EC78 unique RT+ATPase+HNH neighborhood not resolved: {len(good)}')
    rt,neigh=good[0]; out.append(row('EC78_ECO7','CP026202.1','RT_Ec78',rt,'DDBJ_CP026202.1_annotation+311aa+paper_architecture'))
    # coding system envelope: select same-strand candidate RT and nearest downstream/upstream ATPase+nuclease on same strand
    relevant=[]
    for x in neigh:
        txt=(x.get('product','')+' '+x.get('gene','')).lower()
        if x is rt or 'atpase' in txt or 'atp-binding' in txt or 'hnh' in txt or 'nuclease' in txt or 'endonuclease' in txt: relevant.append(x)
    # only near RT within 8 kb and same strand; take contiguous span conservatively
    relevant=[x for x in relevant if x['strand']==rt['strand'] and x['start1']<=rt['end1']+8000 and x['end1']>=rt['start1']-8000]
    env={'start1':min(x['start1'] for x in relevant),'end1':max(x['end1'] for x in relevant),'strand':rt['strand'],'locus_tag':';'.join(x.get('locus_tag','') for x in relevant),'protein_id':'','product':'Ec78 coding-component envelope'}
    out.append(row('EC78_ECO7','CP026202.1','CODING_SYSTEM_ENVELOPE',env,'DDBJ annotation neighborhood; ncRNA 144nt boundary not invented','PARTIAL_SYSTEM'))
    audit['EC78_ECO7']={'rt':rt,'relevant':relevant}
    # DRT2 frozen exact C0 truth anchor; verify annotation overlap and RT identity
    fs=readj(R/'work_DRT2_KPN_features.json'); s,e=598753,600030
    ov=[f for f in fs if f['type']=='CDS' and feat_overlap(f,s,e)]
    exact=[f for f in ov if f['start1']==s and f['end1']==e]
    if len(exact)!=1: raise RuntimeError(f'DRT2 expected exact CDS {s}..{e}, got {len(exact)}')
    rt=exact[0]
    if not (rt_like(rt) or rt.get('protein_id')=='WP_012737279.1'): raise RuntimeError('DRT2 exact CDS identity mismatch')
    out.append(row('DRT2_KPN','AP006725.1','RT_DRT2',rt,'C0_PADLOC_exact_interval_verified_against_DDBJ_AP006725.1'))
    audit['DRT2_KPN']={'rt':rt}
    # write TSV
    cols=['control_id','accession','truth_status','component','start1','end1','start0','end0','strand','locus_tag','protein_id','product','coordinate_source']
    op=R/'outputs/C1B0_CONTROL_TRUTH_RESOLVED.tsv'; op.parent.mkdir(exist_ok=True)
    with op.open('w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=cols,delimiter='\t'); w.writeheader(); w.writerows(out)
    (R/'outputs/C1B0_RESOLUTION_AUDIT.json').write_text(json.dumps(audit,indent=2))
    gate={'gate':'PASS_STAGE04C1B0_RT_ANCHOR_TRUTH_RESOLUTION','gate_pass':True,'truth_resolution_after_score_freeze':True,'score_files_read':False,'controls_resolved_primary_rt_anchor_n':3,'primary_anchor':'RT coding interval','ncRNA_boundaries_invented':False,'retuning_authorized':False,'truth_tsv_sha256':sha(op)}
    (R/'outputs/C1B0_FINAL_GATE.json').write_text(json.dumps(gate,indent=2))
    print(json.dumps(gate,indent=2))
if __name__=='__main__': main()
