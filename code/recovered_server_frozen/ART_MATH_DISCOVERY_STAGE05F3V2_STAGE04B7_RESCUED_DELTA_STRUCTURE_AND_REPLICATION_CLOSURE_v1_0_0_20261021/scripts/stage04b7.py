#!/usr/bin/env python3
import argparse,csv,json,hashlib,os,re,shutil,subprocess,sys
from pathlib import Path
from collections import defaultdict,Counter

EXPECTED={'target18':18,'primary12':12,'edge1':1,'absent5':5,'mixed23':23}
HIST5_DEFAULT={'F3BI37262050','F3BI37356863','F3BI27373723','F3BI55368798','F3BI55915933'}

def tsv_read(p):
    with open(p,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def tsv_write(p,rows,cols):
    p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def fasta_dict(p):
    d={}; name=None; seq=[]
    with open(p) as f:
        for line in f:
            line=line.strip()
            if not line:continue
            if line.startswith('>'):
                if name is not None:d[name]=''.join(seq)
                name=line[1:].split()[0];seq=[]
            else:seq.append(line)
    if name is not None:d[name]=''.join(seq)
    return d

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def find_b6(P):
    hs=sorted(P.glob('ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B6_ARCH41_CROSSWALK_AND_FRAGMENT_RESCUE_v1_0_1_20261020'))
    if not hs: hs=sorted(P.glob('ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B6_ARCH41_CROSSWALK_AND_FRAGMENT_RESCUE_*'))
    good=[]
    for h in hs:
        g=h/'outputs/STAGE04B6_CLOSEOUT_GATE.json'
        if g.exists():
            try:
                z=json.load(open(g))
                if z.get('gate_pass'):good.append(h)
            except:pass
    if not good: raise SystemExit('HOLD Stage04B6 completed root not found')
    return good[-1].resolve()

def load_old_map(path):
    rows=tsv_read(path); hdr=list(rows[0]) if rows else []
    cands=[x for x in hdr if any(k in x.lower() for k in ('candidate_id','instance_id','query_id','member_id'))]
    fams=[x for x in hdr if 'family' in x.lower()]
    best={}; bestcols=None
    for cc in cands:
      for fc in fams:
        d={r.get(cc,''):r.get(fc,'') for r in rows if r.get(cc,'') and r.get(fc,'').startswith('F3BI')}
        if len(d)>len(best):best=d;bestcols=(cc,fc)
    if not best:raise SystemExit('HOLD old family map columns unresolved')
    return best,bestcols

def run_mmseqs(qfaa,tfaa,out,tmp,threads):
    mm=shutil.which('mmseqs')
    if not mm:raise SystemExit('HOLD mmseqs not found')
    if tmp.exists():shutil.rmtree(tmp)
    cmd=[mm,'easy-search',str(qfaa),str(tfaa),str(out),str(tmp),'--min-seq-id','0.3','-c','0.7','--cov-mode','0','-s','7.5','--max-seqs','1000','--threads',str(threads),'--format-output','query,target,fident,qcov,tcov,evalue,bits']
    subprocess.run(cmd,check=True)

def prepare(a):
    P=Path(a.project_root).resolve(); H=Path(a.home).resolve(); O=H/'outputs'; W=H/'work'; I=H/'inputs'
    O.mkdir(exist_ok=True);W.mkdir(exist_ok=True);I.mkdir(exist_ok=True)
    b6=find_b6(P); b6o=b6/'outputs'; b6f=b6/'frozen'
    rg=json.load(open(b6o/'STAGE04B6_RESCUE_GATE.json'))
    if (rg.get('target_n'),rg.get('stop_bounded_both_n'),rg.get('edge_censored_unique_n'),rg.get('absent_n'))!=(18,12,1,5):
        raise SystemExit(f'HOLD unexpected B6 rescue counts {rg}')
    rr=tsv_read(b6o/'DELTA18_UNITIG_RESCUE.tsv')
    if len(rr)!=18:raise SystemExit('HOLD B6 rescue row count')
    # parse rescued FASTA where header is <family>|source_candidate=<cid>
    seqraw=fasta_dict(b6o/'DELTA18_RESCUED_MAXIMAL_STOP_TO_STOP.faa')
    famseq={}
    for hdr,s in seqraw.items():
        fam=hdr.split('|source_candidate=',1)[0]
        famseq[fam]=s
    role=[]
    qmap=[]
    qn=0
    for r in rr:
        st=r['rescue_status']; fam=r['delta_operational_family_id']
        if st=='UNIQUE_STOP_BOUNDED_BOTH': role_name='PRIMARY_COMPLETE_STOP_TO_STOP'
        elif st=='UNIQUE_EDGE_CENSORED': role_name='SENSITIVITY_EDGE_CENSORED'
        else: role_name='UNRESOLVED_NO_STRUCTURE'
        row=dict(r);row['b7_role']=role_name;role.append(row)
        if st.startswith('UNIQUE_'):
            if fam not in famseq:raise SystemExit(f'HOLD rescued FASTA missing {fam}')
            qn+=1;qid=f'B7Q{qn:04d}';seq=famseq[fam]
            qmap.append(dict(query_id=qid,delta_operational_family_id=fam,candidate_id=r['candidate_id'],accession=r['accession'],rescue_status=st,b7_role=role_name,aa_length=len(seq),recovered_motif_c_like=r.get('recovered_motif_c_like',''),sequence=seq))
    if Counter(x['b7_role'] for x in qmap)!=Counter({'PRIMARY_COMPLETE_STOP_TO_STOP':12,'SENSITIVITY_EDGE_CENSORED':1}):raise SystemExit('HOLD B7 role counts')
    tsv_write(O/'DELTA18_STRUCTURE_ROLE_MANIFEST.tsv',role,list(role[0].keys()))
    tsv_write(O/'DELTA13_STRUCTURAL_INPUT_MANIFEST.tsv',[{k:v for k,v in r.items() if k!='sequence'} for r in qmap],[k for k in qmap[0] if k!='sequence'])
    # YAMLs and 4 deterministic shards
    yd=I/'boltz_yaml';shutil.rmtree(yd,ignore_errors=True);yd.mkdir(parents=True)
    for si in range(4):(I/f'shard{si}').mkdir(parents=True,exist_ok=True)
    for idx,r in enumerate(qmap):
        y=yd/f"{r['query_id']}.yaml"
        y.write_text('version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: '+r['sequence']+'\n')
        sy=I/f'shard{idx%4}'/y.name
        if sy.exists() or sy.is_symlink():sy.unlink()
        sy.symlink_to(y.resolve())
    # Mixed23 family-wide replication crosswalk
    pg=json.load(open(b6o/'STAGE04B6_PREPARE_GATE.json'))
    deltafaa=Path(pg['delta_formal_faa']); oldfaa=Path(pg['old_exact_aa_unique_faa']); oldmapfile=Path(pg['old_family_map_source'])
    daa=fasta_dict(deltafaa); oldmap,_=load_old_map(oldmapfile)
    mem=tsv_read(b6f/'B5_SELECTED41_MEMBERS.tsv')
    mix=[r for r in mem if r['arm']=='MIXED_REPLICATION_ARM']
    mixfams=sorted(set(r['delta_operational_family_id'] for r in mix))
    if len(mixfams)!=23:raise SystemExit(f'HOLD mixed family count {len(mixfams)}')
    qfaa=W/'mixed23_all_members.faa'; qlookup={}
    with open(qfaa,'w') as f:
        for i,r in enumerate(mix,1):
            cid=r['candidate_id']
            if cid not in daa:raise SystemExit(f'HOLD delta formal missing {cid}')
            q=f'MX{i:06d}';qlookup[q]=(r['delta_operational_family_id'],cid)
            f.write(f'>{q}\n{daa[cid]}\n')
    m8=W/'mixed23_all_members_vs_old10329.m8';run_mmseqs(qfaa,oldfaa,m8,W/'mmseqs_mixed23_all_tmp',a.threads)
    byfam=defaultdict(list)
    with open(m8) as f:
        for line in f:
            q,t,ident,qcov,tcov,e,bits=line.rstrip('\n').split('\t')[:7]
            if q not in qlookup:continue
            fam,cid=qlookup[q];of=oldmap.get(t,'')
            byfam[fam].append(dict(candidate_id=cid,target=t,old_family=of,fident=float(ident),qcov=float(qcov),tcov=float(tcov),evalue=float(e),bits=float(bits)))
    hist5=set(x.strip() for x in open(b6f/'HISTORICAL5_STRUCTURE_POSITIVE_SYSTEMS.txt') if x.strip()) if (b6f/'HISTORICAL5_STRUCTURE_POSITIVE_SYSTEMS.txt').exists() else HIST5_DEFAULT
    hist7=set(x.strip() for x in open(b6f/'HISTORICAL7_STRICT_RECURRENT_PARTNER_FAMILIES.txt') if x.strip())
    cr=[]
    for fam in mixfams:
        hs=byfam.get(fam,[]); ofs=sorted({h['old_family'] for h in hs if h['old_family']})
        h5=sorted(set(ofs)&hist5);h7=sorted(set(ofs)&hist7)
        cr.append(dict(delta_operational_family_id=fam,delta_member_n=sum(r['delta_operational_family_id']==fam for r in mix),direct_pair_hit_n=len(hs),mapped_old_family_n=len(ofs),historical7_overlap=';'.join(h7),historical5_overlap=';'.join(h5),all_mapped_old_families=';'.join(ofs),familywide_status=('HISTORICAL5_REPLICATION' if h5 else 'HISTORICAL7_COMPARATOR_REPLICATION' if h7 else 'OLD_FAMILY_REPLICATION' if ofs else 'NO_DIRECT_OLD_HIT')))
    tsv_write(O/'MIXED23_FAMILYWIDE_HISTORICAL_CROSSWALK.tsv',cr,list(cr[0].keys()))
    mat=[]
    for h in sorted(hist5):
        fs=[r['delta_operational_family_id'] for r in cr if h in (r['historical5_overlap'].split(';') if r['historical5_overlap'] else [])]
        mat.append(dict(historical_family=h,replicated_by_mixed_family_n=len(fs),replicated=bool(fs),mixed_delta_families=';'.join(fs)))
    tsv_write(O/'HISTORICAL5_REPLICATION_MATRIX.tsv',mat,list(mat[0].keys()))
    gate={'stage':'Stage05F3V2-Stage04B7','gate':'PASS_STAGE04B7_PREPARE_FROZEN','gate_pass':True,'b6_root':str(b6),'primary_complete_n':12,'sensitivity_edge_censored_n':1,'unresolved_absent_n':5,'mixed23_familywide_crosswalk_complete':True,'historical5_replicated_n':sum(x['replicated'] for x in mat),'historical7_represented_n':len({x for r in cr for x in (r['historical7_overlap'].split(';') if r['historical7_overlap'] else [])}),'structure_target_n':13,'rules_changed':False,'retuning_performed':False}
    json.dump(gate,open(O/'STAGE04B7_PREPARE_GATE.json','w'),indent=2);print(json.dumps(gate,indent=2));print('PASS_STAGE04B7_PREPARE_FROZEN')

def hit_tier(tm,qc,tc):
    if tm>=0.50 and qc>=0.50 and tc>=0.50:return 'STRONG'
    if tm>=0.40 and qc>=0.30 and tc>=0.30:return 'SUGGESTIVE'
    return 'WEAK'
ORDER={'NO_HIT':0,'WEAK':1,'SUGGESTIVE':2,'STRONG':3}

def classify(m8):
    hs=[]
    if not Path(m8).exists() or Path(m8).stat().st_size==0:return ('NO_HIT',None,[])
    with open(m8) as f:
        for ln in f:
            z=ln.rstrip('\n').split('\t')
            if len(z)<11:continue
            q,t,e,bits,qc,tc,aln,ql,tl,fi,tm=z[:11]
            h=dict(query=q,target=t,evalue=float(e),bits=float(bits),qcov=float(qc),tcov=float(tc),alnlen=int(float(aln)),qlen=int(float(ql)),tlen=int(float(tl)),fident=float(fi),tm=float(tm))
            h['tier']=hit_tier(h['tm'],h['qcov'],h['tcov']);hs.append(h)
    if not hs:return ('NO_HIT',None,[])
    tier=max((h['tier'] for h in hs),key=lambda x:ORDER[x]);same=[h for h in hs if h['tier']==tier]
    rep=sorted(same,key=lambda h:(-h['tm'],-h['qcov'],-h['tcov'],h['evalue'],h['target']))[0]
    return tier,rep,hs

def finalize(a):
    H=Path(a.home).resolve();O=H/'outputs';M=tsv_read(O/'DELTA13_STRUCTURAL_INPUT_MANIFEST.tsv')
    by={r['query_id']:r for r in M}; rows=[]
    for q,r in by.items():
        m8=H/'work/foldseek'/f'{q}.m8'
        tier,rep,hs=classify(m8)
        role=r['b7_role']
        if role=='PRIMARY_COMPLETE_STOP_TO_STOP':
            auth='PRIMARY_POSITIVE_CLEAN_RT_COMPATIBILITY' if tier in ('STRONG','SUGGESTIVE') else 'PRIMARY_NO_CLEAN_RT_COMPATIBILITY_UNDER_FROZEN_THRESHOLDS'
        else:
            auth='SENSITIVITY_POSITIVE_SUPPORT' if tier in ('STRONG','SUGGESTIVE') else 'SENSITIVITY_NEGATIVE_NOT_AUTHORIZED'
        rows.append(dict(query_id=q,delta_operational_family_id=r['delta_operational_family_id'],candidate_id=r['candidate_id'],rescue_status=r['rescue_status'],b7_role=role,aa_length=r['aa_length'],tier=tier,representative_target='' if rep is None else rep['target'],alntmscore='' if rep is None else f"{rep['tm']:.6g}",qcov='' if rep is None else f"{rep['qcov']:.6g}",tcov='' if rep is None else f"{rep['tcov']:.6g}",evalue='' if rep is None else f"{rep['evalue']:.6g}",returned_hit_n=len(hs),authorized_interpretation=auth))
    rows.sort(key=lambda r:r['query_id']);tsv_write(O/'DELTA13_CORRECTED_STRUCTURAL_COMPATIBILITY.tsv',rows,list(rows[0].keys()))
    prim=[r for r in rows if r['b7_role']=='PRIMARY_COMPLETE_STOP_TO_STOP'];sens=[r for r in rows if r['b7_role']=='SENSITIVITY_EDGE_CENSORED']
    c=Counter(r['tier'] for r in prim);sc=Counter(r['tier'] for r in sens)
    pg=json.load(open(O/'STAGE04B7_PREPARE_GATE.json'))
    gate={'stage':'Stage05F3V2-Stage04B7','gate':'PASS_STAGE04B7_RESCUED_DELTA_STRUCTURE_AND_REPLICATION_CLOSURE','gate_pass':True,'mixed23_familywide_historical5_replicated_n':pg['historical5_replicated_n'],'mixed23_familywide_historical7_represented_n':pg['historical7_represented_n'],'primary_complete_structure_n':len(prim),'primary_strong_n':c['STRONG'],'primary_suggestive_n':c['SUGGESTIVE'],'primary_weak_n':c['WEAK'],'primary_no_hit_n':c['NO_HIT'],'primary_positive_n':c['STRONG']+c['SUGGESTIVE'],'sensitivity_edge_censored_n':len(sens),'sensitivity_positive_n':sc['STRONG']+sc['SUGGESTIVE'],'unresolved_absent_n':5,'classification_rule':'per-hit thresholds then query=max tier over all returned hits','motif_used_for_selection':False,'rules_changed':False,'retuning_performed':False}
    json.dump(gate,open(O/'STAGE04B7_CLOSEOUT_GATE.json','w'),indent=2)
    with open(O/'SHA256SUMS.txt','w') as f:
        for p in sorted(x for x in O.iterdir() if x.is_file() and x.name!='SHA256SUMS.txt'):f.write(f'{sha(p)}  {p.name}\n')
    print(json.dumps(gate,indent=2));print('PASS_STAGE04B7_RESCUED_DELTA_STRUCTURE_AND_REPLICATION_CLOSURE')

def selftest():
    hs=[dict(tm=.55,qcov=.8,tcov=.2,evalue=1e-8,target='A'),dict(tm=.46,qcov=.5,tcov=.4,evalue=1e-5,target='B')]
    tiers=[hit_tier(x['tm'],x['qcov'],x['tcov']) for x in hs]
    assert tiers==['WEAK','SUGGESTIVE']
    assert max(tiers,key=lambda x:ORDER[x])=='SUGGESTIVE'
    assert hit_tier(.51,.51,.51)=='STRONG'
    print('PASS_STAGE04B7_SELFTEST')

def main():
    ap=argparse.ArgumentParser();sp=ap.add_subparsers(dest='cmd',required=True)
    p=sp.add_parser('prepare');p.add_argument('--project-root',required=True);p.add_argument('--home',required=True);p.add_argument('--threads',type=int,default=32)
    f=sp.add_parser('finalize');f.add_argument('--home',required=True)
    sp.add_parser('selftest')
    a=ap.parse_args()
    if a.cmd=='prepare':prepare(a)
    elif a.cmd=='finalize':finalize(a)
    else:selftest()
if __name__=='__main__':main()
