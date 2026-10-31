#!/usr/bin/env python3
import argparse,json,csv,hashlib,re
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

def compact(f):
    return {k:f.get(k) for k in ['start1','end1','strand','locus_tag','gene','product','protein_id','aa_len','has_YADD','has_YxDD']}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); ap.add_argument('--c0',required=True); ap.add_argument('--r3',required=True); a=ap.parse_args()
    R=Path(a.root); out=[]; audit={}
    g=readj(Path(a.r3)/'outputs/STAGE04C0_R3_FINAL_GATE.json')
    if g.get('gate')!='PASS_STAGE04C0_R3_SEQUENCE_LENGTH_CORRECTION_FIXED': raise RuntimeError('R3 gate mismatch')
    c0truth=Path(a.c0)/'truth/CONTROL_LOCUS_TRUTH_SEALED.tsv'
    if not c0truth.exists(): raise RuntimeError('missing original C0 sealed truth')

    # SEN2: literature-linked stable locus tags in CP001363.1.
    fs=readj(R/'work_SEN2_14028S_features.json')
    rt=one_locus(fs,'STM14_4641'); rca=one_locus(fs,'STM14_4640')
    if not rt_like(rt): raise RuntimeError('SEN2 STM14_4641 not annotated RT-like')
    out += [row('SEN2_14028S','CP001363.1','RT_rrtT_STM14_4641',rt,'CP001363.1 locus_tag STM14_4641 + literature'),row('SEN2_14028S','CP001363.1','RcaT_STM14_4640',rca,'CP001363.1 locus_tag STM14_4640 + literature')]
    env={'start1':min(rt['start1'],rca['start1']),'end1':max(rt['end1'],rca['end1']),'strand':rt['strand'],'locus_tag':'STM14_4640;STM14_4641','protein_id':'','product':'RcaT+RT coding envelope'}
    out.append(row('SEN2_14028S','CP001363.1','CODING_SYSTEM_ENVELOPE',env,'resolved coding components; ncRNA boundary not invented','PARTIAL_SYSTEM'))
    audit['SEN2_14028S']={'rt':compact(rt),'rcaT':compact(rca)}

    # EC78/Eco7: GenBank products can be generic/hypothetical. Resolve from the published
    # sequence-level properties of the native ECONIH5 RT: 311 aa and exact YADD catalytic core.
    fs=readj(R/'work_EC78_ECO7_features.json'); cds=[f for f in fs if f['type']=='CDS']
    yadd=[f for f in cds if f.get('has_YADD')]
    near_len=[f for f in yadd if f.get('aa_len') is not None and 310 <= int(f['aa_len']) <= 312]
    exact=[f for f in near_len if int(f.get('aa_len'))==311]
    audit['EC78_ECO7']={
      'all_YADD_CDS_n':len(yadd),
      'YADD_310_312aa_n':len(near_len),
      'YADD_exact_311aa_n':len(exact),
      'YADD_310_312aa_candidates':[compact(x) for x in near_len]
    }
    if len(exact)!=1:
        apath=R/'outputs/C1B0_R1_EC78_CANDIDATE_AUDIT.json'; apath.parent.mkdir(exist_ok=True)
        apath.write_text(json.dumps(audit['EC78_ECO7'],indent=2))
        raise RuntimeError(f'EC78 unique exact 311-aa YADD RT candidate not resolved: {len(exact)}; see {apath}')
    rt=exact[0]
    out.append(row('EC78_ECO7','CP026202.1','RT_Ec78',rt,'CP026202.1 unique 311-aa CDS with exact YADD catalytic motif + published Ec78 RT properties'))
    # Do not invent/derive a PtuA/PtuB coding envelope from generic genome annotations.

    # DRT2: frozen exact C0 PADLOC anchor, verified against source record.
    fs=readj(R/'work_DRT2_KPN_features.json'); s,e=598753,600030
    ov=[f for f in fs if f['type']=='CDS' and feat_overlap(f,s,e)]
    exact2=[f for f in ov if f['start1']==s and f['end1']==e]
    if len(exact2)!=1: raise RuntimeError(f'DRT2 expected exact CDS {s}..{e}, got {len(exact2)}')
    rt2=exact2[0]
    if not (rt_like(rt2) or rt2.get('protein_id')=='WP_012737279.1'): raise RuntimeError('DRT2 exact CDS identity mismatch')
    out.append(row('DRT2_KPN','AP006725.1','RT_DRT2',rt2,'C0 PADLOC exact interval verified against AP006725.1'))
    audit['DRT2_KPN']={'rt':compact(rt2)}

    cols=['control_id','accession','truth_status','component','start1','end1','start0','end0','strand','locus_tag','protein_id','product','coordinate_source']
    op=R/'outputs/C1B0_CONTROL_TRUTH_RESOLVED.tsv'; op.parent.mkdir(exist_ok=True)
    with op.open('w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=cols,delimiter='\t'); w.writeheader(); w.writerows(out)
    (R/'outputs/C1B0_RESOLUTION_AUDIT.json').write_text(json.dumps(audit,indent=2))
    gate={'gate':'PASS_STAGE04C1B0_RT_ANCHOR_TRUTH_RESOLUTION','gate_pass':True,'hotfix_revision':'R1_EC78_SEQUENCE_MOTIF_v1.0.1','truth_resolution_after_score_freeze':True,'score_files_read':False,'controls_resolved_primary_rt_anchor_n':3,'primary_anchor':'RT coding interval','ncRNA_boundaries_invented':False,'ec78_resolution_rule':'unique exact 311-aa CDS with exact YADD motif in CP026202.1','ec78_product_annotation_required':False,'ec78_ptu_annotation_required':False,'retuning_authorized':False,'truth_tsv_sha256':sha(op)}
    (R/'outputs/C1B0_FINAL_GATE.json').write_text(json.dumps(gate,indent=2))
    print(json.dumps(gate,indent=2))
if __name__=='__main__': main()
