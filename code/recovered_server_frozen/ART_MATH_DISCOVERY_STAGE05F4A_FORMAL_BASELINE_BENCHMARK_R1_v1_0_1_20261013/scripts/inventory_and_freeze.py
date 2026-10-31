#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, os, shutil, subprocess
from pathlib import Path
from common import *

EXPECTED_RVT_SHA='74556bd21fd92b947af2759878c3433361b442f2499a62492166a5416971e25e'

START_KEYS=['start0','window_start0','genomic_start0','region_start0','union_start0','merged_start0','proposal_start0','start','window_start','genomic_start','region_start','union_start','merged_start','proposal_start','begin0','begin']
END_KEYS=['end0','window_end0','genomic_end0','region_end0','union_end0','merged_end0','proposal_end0','end','window_end','genomic_end','region_end','union_end','merged_end','proposal_end','stop0','stop']
RANK_KEYS=['rank','region_rank','priority_rank','union_rank','final_rank']

def split_paths(s): return [x for x in str(s).split('|') if x]

def first_existing(root:Path, rels):
    for rel in rels:
        p=root/rel
        if p.is_file(): return p
    return None

def numeric(v):
    try: return float(str(v).strip())
    except Exception: return None

def pick_coord(row, keys):
    low={str(k).strip().lower():v for k,v in row.items()}
    for k in keys:
        if k in low:
            z=numeric(low[k])
            if z is not None: return int(z)
    # cautious fallback: exact suffix only, excluding truth/rt/partner fields
    bad=('truth','rt_','partner','motif','lag')
    for k,v in low.items():
        if any(b in k for b in bad): continue
        if keys is START_KEYS and (k.endswith('_start0') or k.endswith('_start')):
            z=numeric(v)
            if z is not None:return int(z)
        if keys is END_KEYS and (k.endswith('_end0') or k.endswith('_end')):
            z=numeric(v)
            if z is not None:return int(z)
    return None

def family_of(row):
    txt=' '.join([str(k)+'='+str(v) for k,v in row.items()]).lower()
    a='acnv' in txt
    r='recurr' in txt
    if a and r:return 'ACNV+RECURRENCE'
    if r:return 'RECURRENCE'
    if a:return 'ACNV'
    return 'OTHER'

def q_of(row,fam):
    low={str(k).strip().lower():v for k,v in row.items()}
    preferred=[]
    if 'RECURRENCE' in fam:
        preferred += [k for k in low if 'q' in k and 'recurr' in k]
    if 'ACNV' in fam:
        preferred += [k for k in low if 'q' in k and 'acnv' in k]
    preferred += [k for k in ['q','qvalue','q_value','family_q','proposal_q'] if k in low]
    preferred += [k for k in low if ('qvalue' in k or 'q_value' in k or k.endswith('_q'))]
    seen=set()
    for k in preferred:
        if k in seen:continue
        seen.add(k); z=numeric(low.get(k,''))
        if z is not None:return z
    return 1.0

def rank_of(row,idx):
    low={str(k).strip().lower():v for k,v in row.items()}
    for k in RANK_KEYS:
        if k in low:
            z=numeric(low[k])
            if z is not None and z>=1:return int(z)
    return idx

def overlap(a,b,c,d): return a<d and b>c

def canonicalize_proposals(src:Path,dst:Path,cid:str):
    raw=read_tsv(src); out=[]; bad=0
    for i,r in enumerate(raw,1):
        s=pick_coord(r,START_KEYS); e=pick_coord(r,END_KEYS)
        if s is None or e is None or e<=s:
            bad+=1; continue
        fam=family_of(r); q=q_of(r,fam)
        out.append({'case_id':cid,'source_row':i,'start0':s,'end0':e,'family':fam,'qvalue':q,'source_path':str(src),'source_sha256':sha256_file(src)})
    fields=['case_id','source_row','start0','end0','family','qvalue','source_path','source_sha256']
    write_tsv(dst,out,fields)
    return out,{'source_rows':len(raw),'parsed_rows':len(out),'unparsed_rows':bad,'acnv_rows':sum('ACNV' in x['family'] for x in out),'recurrence_rows':sum('RECURRENCE' in x['family'] for x in out)}

def canonicalize_union(src:Path,dst:Path,cid:str,props):
    raw=read_tsv(src); out=[]; bad=0
    ac=[x for x in props if 'ACNV' in x['family']]
    rr=[x for x in props if 'RECURRENCE' in x['family']]
    h=sha256_file(src)
    for i,r in enumerate(raw,1):
        s=pick_coord(r,START_KEYS); e=pick_coord(r,END_KEYS)
        if s is None or e is None or e<=s:
            bad+=1; continue
        ra=rank_of(r,i)
        ha=int(any(overlap(s,e,x['start0'],x['end0']) for x in ac))
        hr=int(any(overlap(s,e,x['start0'],x['end0']) for x in rr))
        out.append({'case_id':cid,'source_row':i,'source_rank':ra,'start0':s,'end0':e,'has_acnv':ha,'has_recurrence':hr,'dual_arm':int(ha and hr),'source_path':str(src),'source_sha256':h})
    out.sort(key=lambda x:(x['source_rank'],x['source_row']))
    fields=['case_id','source_row','source_rank','start0','end0','has_acnv','has_recurrence','dual_arm','source_path','source_sha256']
    write_tsv(dst,out,fields)
    return out,{'source_rows':len(raw),'parsed_rows':len(out),'unparsed_rows':bad,'dual_arm_rows':sum(x['dual_arm'] for x in out)}

def version_line(path:Path,name):
    for args in ([str(path),'--version'],[str(path),'-v'],[str(path),'-h']):
        try:
            q=subprocess.run(args,text=True,capture_output=True,timeout=10)
            s=(q.stdout+'\n'+q.stderr).strip().splitlines()
            if s:return s[0][:300]
        except Exception:pass
    return None

def resolve_tool(name,env_name=None):
    cand=[]
    if env_name and os.environ.get(env_name): cand.append(Path(os.environ[env_name]))
    q=shutil.which(name)
    if q:cand.append(Path(q))
    for pat in [f'/home/hgq/conda-envs/*/bin/{name}',f'/home/hgq/.conda/envs/*/bin/{name}',f'/opt/anaconda3/envs/*/bin/{name}']:
        cand += sorted(Path('/').glob(pat.lstrip('/')))
    seen=set()
    for p in cand:
        try:p=p.resolve()
        except Exception:continue
        if str(p) in seen:continue
        seen.add(str(p))
        if p.is_file() and os.access(p,os.X_OK):
            return {'present':True,'path':str(p),'sha256':sha256_file(p),'version':version_line(p,name)}
    return {'present':False,'path':None,'sha256':None,'version':None}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('project_root'); a=ap.parse_args()
    root=Path(a.project_root).resolve(); here=Path(__file__).resolve().parents[1]
    cases=read_tsv(here/'config/CASES.tsv'); policy={r['case_id']:r for r in read_tsv(here/'config/SOURCE_POLICY.tsv')}
    out=here/'frozen'; canroot=out/'canonical'; canroot.mkdir(parents=True,exist_ok=True)
    seq_rows=[]; board_rows=[]; all_seq=True; all_board=True
    for c in cases:
        cid=c['case_id']; pol=policy[cid]
        sp=first_existing(root,split_paths(pol['sequence_relpaths']))
        status='MISSING'; hdr=''; seq=''; fsha=''; ssha=''; n=''
        if sp:
            try:
                rec=list(fasta_records(sp))
                if len(rec)==1:
                    hdr,seq=rec[0]; n=len(seq); fsha=sha256_file(sp); ssha=seq_sha(seq); status='RESOLVED_EXACT_POLICY'
                else: status='INVALID_MULTIRECORD_FASTA'
            except Exception as e: status='INVALID_FASTA:'+type(e).__name__
        if status!='RESOLVED_EXACT_POLICY':all_seq=False
        seq_rows.append({'case_id':cid,'reference_accession':c['reference_accession'],'sequence_accession':c['sequence_accession'],'status':status,'chosen_fasta':str(sp) if sp else '','fasta_sha256':fsha,'sequence_sha256':ssha,'length_bp':n,'source_policy':'EXPLICIT_PRECEDENCE'})
        pp=first_existing(root,split_paths(pol['proposal_relpaths'])); up=first_existing(root,split_paths(pol['union_relpaths']))
        cdir=canroot/cid; cdir.mkdir(parents=True,exist_ok=True)
        pst={'source_rows':0,'parsed_rows':0,'unparsed_rows':0,'acnv_rows':0,'recurrence_rows':0}; ust={'source_rows':0,'parsed_rows':0,'unparsed_rows':0,'dual_arm_rows':0}
        props=[]; unions=[]; bstat='MISSING'
        try:
            if pp and up:
                props,pst=canonicalize_proposals(pp,cdir/'CANONICAL_PROPOSALS.tsv',cid)
                unions,ust=canonicalize_union(up,cdir/'CANONICAL_UNION_REGIONS.tsv',cid,props)
                bstat='RESOLVED_CANONICAL' if pst['source_rows']>0 and pst['parsed_rows']>0 and ust['source_rows']>0 and ust['parsed_rows']>0 else 'PARSE_HOLD'
            else:bstat='SOURCE_MISSING'
        except Exception as e:bstat='PARSE_ERROR:'+type(e).__name__
        if bstat!='RESOLVED_CANONICAL':all_board=False
        board_rows.append({'case_id':cid,'reference_accession':c['reference_accession'],'sequence_accession':c['sequence_accession'],'proposal_source':str(pp) if pp else '','proposal_sha256':sha256_file(pp) if pp else '','union_source':str(up) if up else '','union_sha256':sha256_file(up) if up else '','proposal_source_rows':pst['source_rows'],'proposal_parsed_rows':pst['parsed_rows'],'acnv_proposal_n':pst['acnv_rows'],'recurrence_proposal_n':pst['recurrence_rows'],'union_region_n':ust['parsed_rows'],'dual_arm_region_n':ust['dual_arm_rows'],'board_status':bstat})
    write_tsv(out/'BASELINE_SEQUENCE_MANIFEST.tsv',seq_rows,list(seq_rows[0]))
    write_tsv(out/'BASELINE_BOARD_INVENTORY.tsv',board_rows,list(board_rows[0]))

    hmm_candidates=[
      root/'ART_MATH_DISCOVERY_STAGE05D1_NOVELTY_FIRST_MINING_v1_0_0_20261006/references/myRT/RVT-All.hmm',
      root/'ART_MATH_DISCOVERY_STAGE05D4R5B_RT_FOCUSED_STRUCTURAL_PANEL_v1_0_0_20261007/references/RVT-All.hmm',
      root/'ART_MATH_DISCOVERY_STAGE05F3B_DE_NOVO_ANNOTATION_INDEPENDENT_RT_DISCOVERY_v1_0_1_20261011/references/RVT-All.hmm',
      root/'ART_MATH_DISCOVERY_STAGE05F3B_DE_NOVO_ANNOTATION_INDEPENDENT_RT_DISCOVERY_v1_0_0_20261010/references/RVT-All.hmm']
    rvt_hits=[]
    for p in hmm_candidates:
        if p.is_file():rvt_hits.append((p,sha256_file(p)))
    exact=[x for x in rvt_hits if x[1]==EXPECTED_RVT_SHA]; rvt=exact[0] if exact else None
    tools={'python':resolve_tool('python'),'prodigal':resolve_tool('prodigal','F4A_PRODIGAL'),'hmmscan':resolve_tool('hmmscan','F4A_HMMSCAN'),'hmmpress':resolve_tool('hmmpress'),'trf':resolve_tool('trf')}
    tool_obj={'rvt_status':'RESOLVED' if rvt else 'MISSING_EXPECTED_SHA','rvt_expected_sha256':EXPECTED_RVT_SHA,'rvt_path':str(rvt[0]) if rvt else None,'rvt_sha256':rvt[1] if rvt else None,'all_rvt_candidates':[{'path':str(p),'sha256':h} for p,h in rvt_hits],'tools':tools}
    safe_json(out/'BASELINE_TOOLCHAIN.json',tool_obj)
    truth_fields=['case_id','reference_accession','sequence_accession','role','group','primary_eval','truth_mode','truth_start0','truth_end0','notes']
    write_tsv(out/'BASELINE_TRUTH_SEAL.tsv',cases,truth_fields)
    gate_pass=all_seq and all_board and bool(rvt) and tools['prodigal']['present'] and tools['hmmscan']['present']
    gate={'stage':'Stage05F4A-B0-R1-input-board-toolchain-freeze','version':'v1.0.1','gate':'PASS_STAGE05F4A_B0_R1_INPUTS_BOARDS_TOOLCHAIN_FROZEN_SCORE_READY' if gate_pass else 'HOLD_STAGE05F4A_B0_R1_INPUTS_BOARDS_OR_TOOLCHAIN_INCOMPLETE','gate_pass':gate_pass,'project_root':str(root),'all_case_sequences_resolved':all_seq,'all_historical_boards_canonicalized':all_board,'rvt_expected_profile_resolved':bool(rvt),'prodigal_present':tools['prodigal']['present'],'prodigal_path':tools['prodigal']['path'],'hmmscan_present':tools['hmmscan']['present'],'hmmscan_path':tools['hmmscan']['path'],'trf_optional_present':tools['trf']['present'],'truth_opened_for_scoring':False,'retuning_authorized':False,'hotfix_scope':'exact-source precedence + corrected Stage04C1 scored-source accession binding + C1_SELECTED_PROPOSALS support + canonical frozen-board extraction + exact runtime binding; no scientific score/threshold changes'}
    safe_json(out/'BASELINE_B0_GATE.json',gate); print(json.dumps(gate,indent=2))
if __name__=='__main__':main()
