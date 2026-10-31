#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, shutil, subprocess
from pathlib import Path
from common import *

FA_EXT={'.fa','.fna','.fasta','.fas'}
EXPECTED_RVT_SHA='74556bd21fd92b947af2759878c3433361b442f2499a62492166a5416971e25e'


def tool_info(name):
    path=shutil.which(name)
    if not path: return {'present':False,'path':None,'version':None}
    ver=None
    for args in ([name,'--version'],[name,'-v'],[name,'-h']):
        try:
            q=subprocess.run(args,text=True,capture_output=True,timeout=10)
            s=(q.stdout+'\n'+q.stderr).strip().splitlines()
            if s: ver=s[0][:300]; break
        except Exception: pass
    return {'present':True,'path':path,'version':ver}


def case_keys(c):
    return {c['case_id'].upper(),c['accession'].upper(),c['accession'].split('.')[0].upper()}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('project_root'); a=ap.parse_args()
    root=Path(a.project_root).resolve(); here=Path(__file__).resolve().parents[1]
    cases=read_tsv(here/'config/CASES.tsv')
    out=here/'frozen'; out.mkdir(exist_ok=True)
    keymap={c['case_id']:case_keys(c) for c in cases}
    fasta_hits={c['case_id']:[] for c in cases}
    prop_hits={c['case_id']:[] for c in cases}
    union_hits={c['case_id']:[] for c in cases}
    rvt_hits=[]

    # Single project-tree walk. Only candidate filenames/paths are opened.
    for p in root.rglob('*'):
        if not p.is_file(): continue
        up=str(p).upper(); suf=p.suffix.lower(); name=p.name.upper()
        if p.name=='RVT-All.hmm':
            try:rvt_hits.append((p,sha256_file(p)))
            except Exception: pass
            continue
        if suf in FA_EXT:
            matched=[]
            for c in cases:
                if any(k in up for k in keymap[c['case_id']]): matched.append(c)
            if not matched: continue
            try: rec=list(fasta_records(p))
            except Exception: continue
            if len(rec)!=1: continue
            hdr,seq=rec[0]; sh=seq_sha(seq); fh=sha256_file(p)
            for c in matched:
                if c['case_id']=='L0050':
                    if len(seq)!=228907: continue
                else:
                    a0=c['accession'].split('.')[0].upper()
                    if a0 not in hdr.upper() and a0 not in p.name.upper(): continue
                fasta_hits[c['case_id']].append((p,hdr,len(seq),sh,fh))
            continue
        is_prop=('SELECTED_PROPOSAL_WINDOWS' in name and suf=='.tsv')
        is_union=('UNION_REGIONS' in name and suf=='.tsv')
        if not (is_prop or is_union): continue
        # first use path match; only then inspect prefix once if needed
        matched_ids=[c['case_id'] for c in cases if any(k in up for k in keymap[c['case_id']])]
        txt=''
        if not matched_ids:
            try: txt=p.read_text(errors='replace')[:200000].upper()
            except Exception: txt=''
            matched_ids=[c['case_id'] for c in cases if any(k in txt for k in keymap[c['case_id']])]
        for cid in matched_ids:
            (prop_hits if is_prop else union_hits)[cid].append(p)

    seq_rows=[]; board_rows=[]; all_seq_ok=True
    for c in cases:
        cand=fasta_hits[c['case_id']]; seqs={x[3] for x in cand}
        status='MISSING' if not cand else ('AMBIGUOUS_SEQUENCE' if len(seqs)>1 else 'RESOLVED')
        ch=sorted(cand,key=lambda x:str(x[0]))[0] if status=='RESOLVED' else None
        if status!='RESOLVED': all_seq_ok=False
        seq_rows.append({'case_id':c['case_id'],'accession':c['accession'],'status':status,
            'chosen_fasta':str(ch[0]) if ch else '','fasta_sha256':ch[4] if ch else '',
            'sequence_sha256':ch[3] if ch else '','length_bp':ch[2] if ch else '',
            'candidate_copy_n':len(cand)})
        ps=sorted(set(prop_hits[c['case_id']])); us=sorted(set(union_hits[c['case_id']]))
        board_rows.append({'case_id':c['case_id'],'accession':c['accession'],
            'proposal_file_n':len(ps),'proposal_files':'|'.join(map(str,ps)),
            'union_file_n':len(us),'union_files':'|'.join(map(str,us)),
            'board_status':'RESOLVED' if ps and us else 'MISSING_OR_PARTIAL'})
    write_tsv(out/'BASELINE_SEQUENCE_MANIFEST.tsv',seq_rows,list(seq_rows[0]))
    write_tsv(out/'BASELINE_BOARD_INVENTORY.tsv',board_rows,list(board_rows[0]))

    exact=[x for x in rvt_hits if x[1]==EXPECTED_RVT_SHA]
    rvt_st='RESOLVED' if exact else ('MISSING_EXPECTED_SHA' if rvt_hits else 'MISSING')
    rvt_ch=sorted(exact,key=lambda x:str(x[0]))[0] if exact else None
    tools={n:tool_info(n) for n in ['python','prodigal','hmmscan','hmmpress','trf']}
    tool_obj={'rvt_status':rvt_st,'rvt_expected_sha256':EXPECTED_RVT_SHA,
      'rvt_path':str(rvt_ch[0]) if rvt_ch else None,'rvt_sha256':rvt_ch[1] if rvt_ch else None,
      'all_rvt_candidates':[{'path':str(p),'sha256':h} for p,h in rvt_hits],'tools':tools}
    safe_json(out/'BASELINE_TOOLCHAIN.json',tool_obj)
    truth_fields=['case_id','accession','role','group','primary_eval','truth_mode','truth_start0','truth_end0','notes']
    write_tsv(out/'BASELINE_TRUTH_SEAL.tsv',cases,truth_fields)
    gate_pass=all_seq_ok and rvt_st=='RESOLVED' and tools['prodigal']['present'] and tools['hmmscan']['present']
    gate={'stage':'Stage05F4A-B0-input-and-toolchain-freeze',
      'gate':'PASS_STAGE05F4A_B0_INPUTS_FROZEN_SCORE_READY' if gate_pass else 'HOLD_STAGE05F4A_B0_INPUTS_OR_TOOLCHAIN_INCOMPLETE',
      'gate_pass':gate_pass,'project_root':str(root),'single_tree_walk':True,
      'all_case_sequences_resolved':all_seq_ok,'rvt_expected_profile_resolved':rvt_st=='RESOLVED',
      'prodigal_present':tools['prodigal']['present'],'hmmscan_present':tools['hmmscan']['present'],
      'trf_optional_present':tools['trf']['present'],'truth_opened_for_scoring':False,
      'retuning_authorized':False,
      'note':'Board inventory may be partial; score-freeze can run only after sequence/HMM/toolchain PASS. Frozen math-board comparison is finalized post-score-freeze.'}
    safe_json(out/'BASELINE_B0_GATE.json',gate); print(json.dumps(gate,indent=2))
if __name__=='__main__': main()
