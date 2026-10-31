#!/usr/bin/env python3
import argparse, csv, hashlib, json, os, re, subprocess, sys, urllib.request
from pathlib import Path
from collections import Counter, defaultdict

STRICT_MIN_UPSTREAM=264
PARTNER_MIN_DOWNSTREAM=93
K=31

def read_tsv(p):
    with open(p, newline='') as f: return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(p, rows, fields):
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=fields,extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def rc(s):
    tab=str.maketrans('ACGTNacgtn','TGCANtgcan')
    return s.translate(tab)[::-1]

def parse_cid(cid):
    z=cid.split('|')
    if len(z)<7: raise ValueError(f'candidate_id format: {cid}')
    acc,contig,strand,frame,coords,aa,ordv=z[:7]
    a,b=map(int,coords.split('-'))
    return dict(candidate_id=cid,accession=acc,contig_id=contig,strand=strand,start0=a,end0=b,aa_length=int(aa[2:]))

def open_text(p):
    p=str(p)
    if p.endswith('.zst'):
        return subprocess.Popen(['zstdcat',p],stdout=subprocess.PIPE,text=True,bufsize=1)
    if p.endswith('.gz'):
        import gzip
        return gzip.open(p,'rt')
    return open(p)

def iter_fasta(p):
    proc=None; fh=None
    try:
        if str(p).endswith('.zst'):
            proc=subprocess.Popen(['zstdcat',str(p)],stdout=subprocess.PIPE,text=True,bufsize=1)
            fh=proc.stdout
        elif str(p).endswith('.gz'):
            import gzip; fh=gzip.open(p,'rt')
        else: fh=open(p)
        name=None; seq=[]; header=''
        for line in fh:
            line=line.rstrip('\n')
            if line.startswith('>'):
                if name is not None: yield name,''.join(seq).upper(),header
                header=line[1:]
                name=header.split()[0]; seq=[]
            else: seq.append(line.strip())
        if name is not None: yield name,''.join(seq).upper(),header
    finally:
        if fh is not None:
            try: fh.close()
            except Exception: pass
        if proc is not None:
            try: proc.wait(timeout=10)
            except Exception:
                try: proc.terminate(); proc.wait(timeout=3)
                except Exception:
                    try: proc.kill()
                    except Exception: pass

def find_b4b4_results(root):
    cands=sorted(root.glob('ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B4_NATIVE_ARCHITECTURE_SATURATION_CLOSEOUT_v1_0_1_20261020/outputs'))
    if not cands:
        cands=sorted(root.glob('ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B4_NATIVE_ARCHITECTURE_SATURATION_CLOSEOUT*/outputs'))
    for p in reversed(cands):
        if (p/'DELTA12_NATIVE_CONTEXTS.tsv').exists() and (p/'DELTA12_STRUCTURAL_POSITIVES.tsv').exists(): return p
    raise FileNotFoundError('Stage04B4 outputs not found')

def find_old_manifest(root):
    cands=[
      root/'external_db/F3_LARGE_UNIVERSE_20261009/F3A_EXTERNAL_SOURCE_MANIFEST.tsv',
      root/'ART_MATH_DISCOVERY_STAGE05F3A_LARGE_UNIVERSE_FROZEN_QUERY_EXPANSION_v1_0_0_20261009/outputs/STAGE05F3A0_EXTERNAL_SOURCE_MANIFEST.tsv'
    ]
    for p in cands:
        if p.exists(): return p
    xs=list(root.glob('**/F3A_EXTERNAL_SOURCE_MANIFEST.tsv'))+list(root.glob('**/STAGE05F3A0_EXTERNAL_SOURCE_MANIFEST.tsv'))
    if xs: return sorted(xs)[0]
    return None

def find_new_ledger(root):
    cands=list(root.glob('ART_MATH_DISCOVERY_STAGE05F3V2_60G_DELTA_DISTRIBUTED_ACQUISITION_R2_v1_0_2_20261013/frozen/F3V2_DELTA_SELECTED_ACCESSIONS.tsv'))
    if not cands: cands=list(root.glob('**/frozen/F3V2_DELTA_SELECTED_ACCESSIONS.tsv'))
    return sorted(cands)[-1] if cands else None

def numeric(v):
    try:return int(float(str(v).replace(',','')))
    except:return None

def infer_bp(rows):
    if not rows:return None,None
    keys=list(rows[0])
    prefs=['realized_bp','nucleotide_bp','accepted_bp','contigs_sumlen','seqstats_contigs_sumlen','bp','total_bp']
    for k in prefs:
        if k in keys:
            vals=[numeric(r.get(k,'')) for r in rows]; vals=[x for x in vals if x is not None]
            if vals:return sum(vals),k
    for k in keys:
        lk=k.lower()
        if ('bp' in lk or 'sumlen' in lk) and not any(x in lk for x in ['target','cumulative','rank']):
            vals=[numeric(r.get(k,'')) for r in rows]; vals=[x for x in vals if x is not None]
            if vals:return sum(vals),k
    return None,None

def get_acc_key(rows):
    if not rows:return None
    for k in ['accession','run_accession','sra_accession','source_accession']:
        if k in rows[0]: return k
    return None

def summarize_manifest(p, fallback_bp=None):
    if p is None:return {'path':'','n':None,'bp':fallback_bp,'bp_field':'','columns':[]}
    rows=read_tsv(p); ak=get_acc_key(rows)
    n=len({r[ak] for r in rows if ak and r.get(ak)}) if ak else len(rows)
    bp,bpf=infer_bp(rows)
    if bp is None:bp=fallback_bp
    out={'path':str(p),'n':n,'bp':bp,'bp_field':bpf or 'fallback','columns':list(rows[0]) if rows else []}
    for k in ['stratum','source_stratum','bioproject','BioProject','biosample','BioSample','assay_type','organism','taxonomy']:
        if rows and k in rows[0]:
            c=Counter(r.get(k,'') for r in rows if r.get(k,'')); out[k+'_top']=c.most_common(10)
    return out

def find_old_panel(root):
    xs=list(root.glob('ART_MATH_DISCOVERY_STAGE05F3B3C*/frozen/F3B3C_QUERY_MANIFEST.tsv'))
    if not xs: xs=list(root.glob('**/F3B3C_QUERY_MANIFEST.tsv'))
    return sorted(xs)[-1] if xs else None

def find_delta_status(root):
    xs=list(root.glob('ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B3_DELTA32*/technical_repair_r2/F3V2_STAGE04B3_BOLTZ_STATUS_CLOSURE32.tsv'))
    if not xs: xs=list(root.glob('ART_MATH_DISCOVERY_STAGE05F3V2_STAGE04B3_DELTA32*/outputs/F3V2_STAGE04B3_BOLTZ_STATUS.tsv'))
    if not xs: xs=list(root.glob('**/F3V2_STAGE04B3_BOLTZ_STATUS_CLOSURE32.tsv'))
    return sorted(xs)[-1] if xs else None

def extract_panel(p, label):
    if p is None:return []
    rows=read_tsv(p); out=[]
    for i,r in enumerate(rows,1):
        cid=r.get('candidate_id') or r.get('query_id') or r.get('id') or ''
        if not cid or '|' not in cid: continue
        pc=parse_cid(cid)
        aa=numeric(r.get('aa_length','')) or pc['aa_length']
        out.append({'panel':label,'rank':numeric(r.get('shortlist_rank','')) or i,'family_rep_id':r.get('family_rep_id',''),'candidate_id':cid,'accession':pc['accession'],'aa_length':aa})
    return out

def bins(rows):
    c=Counter()
    for r in rows:
        L=int(r['aa_length'])
        c['<=128']+=L<=128; c['129-899']+=129<=L<=899; c['>=900']+=L>=900
    return dict(c)

def head_size(url):
    try:
        req=urllib.request.Request(url,method='HEAD',headers={'User-Agent':'Mozilla/5.0'})
        with urllib.request.urlopen(req,timeout=20) as resp:
            return int(resp.headers.get('Content-Length') or 0), int(getattr(resp,'status',200)), ''
    except Exception as e:return 0,0,str(e)

def search_anchors(unitig_file, anchors):
    # anchors dict side->sequence in candidate-oriented direction. Search exact seq or rc.
    hits={s:[] for s in anchors}
    for uid,seq,header in iter_fasta(unitig_file):
        for side,a in anchors.items():
            if not a: continue
            pos=seq.find(a)
            if pos>=0:
                ext=pos if side=='UPSTREAM' else len(seq)-(pos+len(a))
                hits[side].append((uid,'FWD',pos,ext,len(seq),header))
            ar=rc(a); pos2=seq.find(ar)
            if pos2>=0:
                # orient reverse-complement unitig so anchor equals candidate-oriented anchor
                olen=len(seq); opos=olen-(pos2+len(a))
                ext=opos if side=='UPSTREAM' else olen-(opos+len(a))
                hits[side].append((uid,'REV',opos,ext,olen,header))
    return hits

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--project-root',required=True)
    ap.add_argument('--outdir',required=True)
    ap.add_argument('--unitig-dir',default='')
    ap.add_argument('--network-head',action='store_true')
    args=ap.parse_args()
    root=Path(args.project_root).resolve(); out=Path(args.outdir).resolve(); out.mkdir(parents=True,exist_ok=True)
    b4=find_b4b4_results(root)
    pos=read_tsv(b4/'DELTA12_STRUCTURAL_POSITIVES.tsv')
    mem=read_tsv(b4/'DELTA12_FAMILY_MEMBERSHIP.tsv')
    ctx=read_tsv(b4/'DELTA12_NATIVE_CONTEXTS.tsv')
    cmap={r['candidate_id']:r['family_rep_id'] for r in mem}
    ctxmap={r['candidate_id']:r for r in ctx}
    reps=[]
    for r in pos:
        cid=r['candidate_id']; pc=parse_cid(cid); c=ctxmap[cid]
        clen=int(c['context_length_bp']); rt0=int(c['rt_start0_in_oriented_context']); rt1=int(c['rt_end0_in_oriented_context'])
        up=rt0; down=clen-rt1
        reps.append({
          'family_rep_id':cmap.get(cid,''),'candidate_id':cid,'accession':pc['accession'],'tier':r['tier'],'aa_length':pc['aa_length'],
          'current_context_bp':clen,'current_upstream_bp':up,'current_downstream_bp':down,
          'strict_array_evaluable_current':str(up>=STRICT_MIN_UPSTREAM),
          'partner_evaluable_current':str(down>=PARTNER_MIN_DOWNSTREAM),
          'both_flanks_500bp':str(up>=500 and down>=500),'both_flanks_2000bp':str(up>=2000 and down>=2000),
          'both_flanks_5000bp':str(up>=5000 and down>=5000),'both_flanks_10000bp':str(up>=10000 and down>=10000)
        })
    repfields=list(reps[0]); write_tsv(out/'DELTA12_REP_CONTEXT_COMPLETENESS.tsv',reps,repfields)

    # design audit
    oldm=find_old_manifest(root); newl=find_new_ledger(root)
    oldsum=summarize_manifest(oldm,6121620524); newsum=summarize_manifest(newl,54377819157)
    oldp=extract_panel(find_old_panel(root),'INITIAL_6P12GB_PRIORITY32')
    newp=extract_panel(find_delta_status(root),'DELTA_54G_PRIORITY32')
    write_tsv(out/'SIX_VS_54G_PRIORITY_PANEL_LENGTHS.tsv',oldp+newp,['panel','rank','family_rep_id','candidate_id','accession','aa_length'])
    # lightweight audit of whether the final delta acquisition implementation explicitly references contiguity metrics
    acq_dirs=list(root.glob('ART_MATH_DISCOVERY_STAGE05F3V2_60G_DELTA_DISTRIBUTED_ACQUISITION_R2_v1_0_2_20261013'))
    contig_tokens=[]
    for d in acq_dirs:
        for q in list(d.rglob('*.py'))+list(d.rglob('*.sh'))+list(d.rglob('*.json')):
            try: txt=q.read_text(errors='ignore').lower()
            except Exception: continue
            hits=[tok for tok in ['seqstats_contigs_n50','contigs_n50','seqstats_contigs_maxlen','contigs_maxlen','n50','maxlen'] if tok in txt]
            if hits: contig_tokens.append({'file':str(q),'tokens':sorted(set(hits))})
    audit={
      'stage':'Stage05F3V2-Stage04B4R1',
      'old_source_manifest':oldsum,'delta_source_ledger':newsum,
      'old_frozen_contig_n_from_project_closeout':8920553,
      'old_mean_contig_bp_from_frozen_totals':6121620524/8920553,
      'delta_acquisition_contiguity_tokens_detected':contig_tokens,
      'delta_selected_ledger_has_contiguity_columns': any(('n50' in c.lower() or 'maxlen' in c.lower()) for c in newsum.get('columns',[])),
      'old_bp_per_accession': (oldsum['bp']/oldsum['n'] if oldsum.get('bp') and oldsum.get('n') else None),
      'delta_bp_per_accession': (newsum['bp']/newsum['n'] if newsum.get('bp') and newsum.get('n') else None),
      'delta_vs_old_bp_per_accession_ratio': ((newsum['bp']/newsum['n'])/(oldsum['bp']/oldsum['n']) if all([oldsum.get('bp'),oldsum.get('n'),newsum.get('bp'),newsum.get('n')]) else None),
      'old_priority_length_bins':bins(oldp),'delta_priority_length_bins':bins(newp),
      'delta_representative_context':{
        'n':len(reps),
        'strict_array_evaluable_current_n':sum(r['strict_array_evaluable_current']=='True' for r in reps),
        'partner_evaluable_current_n':sum(r['partner_evaluable_current']=='True' for r in reps),
        'both_flanks_500bp_n':sum(r['both_flanks_500bp']=='True' for r in reps),
        'both_flanks_2000bp_n':sum(r['both_flanks_2000bp']=='True' for r in reps)
      },
      'method_non_equivalences':[
        'The delta experiment excludes old-delta mixed families before prospective priority selection; the initial discovery universe did not have this delta-only restriction.',
        'The 6.12-Gb final ART-like systems were obtained after exhaustive repeat/partner architecture census over the formal RT-locus universe, whereas the 54-Gb native-architecture closeout was applied only after a 32-family protein-priority and structural-compatibility funnel.',
        'The same four-axis 32-family priority rule is therefore a valid frozen protein-space transfer test but is not an architecture-prevalence replication design.',
        'The delta representative panel has a markedly different length distribution if the audited bins differ, which can alter both structural coverage and susceptibility to contig-edge truncation.'
      ],
      'old_structural16_context_completeness_from_frozen_project_closeout':{
        'family_with_at_least_one_both_flanks_500bp':'11/16',
        'family_with_at_least_one_both_flanks_2kb':'8/16',
        'family_with_at_least_one_both_flanks_5kb':'6/16',
        'family_with_at_least_one_both_flanks_10kb':'3/16',
        'note':'Recorded from the frozen 6.12-Gb context-completeness closeout; not recomputed by this lightweight stage.'
      },
      'provisional_interpretation':{
        '54G_protein_space_value':'VALID_AND_USEFUL',
        '54G_family_saturation_value':'VALID_AND_USEFUL',
        '54G_direct_architecture_replication':'NOT_ESTABLISHED_CONTEXT_CENSORED_AND_PIPELINE_NON_EQUIVALENT',
        'source_selection_problem_status':'POSSIBLE_FOR_ARCHITECTURE_NOT_YET_PROVEN_GLOBALLY',
        'recommended_next_step':'12-representative context/unitig preflight first; if context is recoverable, replay frozen architecture rules. For a true architecture replication, consider an exhaustive frozen repeat/partner census across the 54-Gb formal loci rather than only the 12 structural positives, and only then decide whether a contiguity-matched new source panel is required.'
      }
    }
    (out/'SIX_VS_54G_DESIGN_AUDIT.json').write_text(json.dumps(audit,indent=2)+'\n')

    # unitig preflight plan and optional local anchor search
    unitdir=Path(args.unitig_dir).resolve() if args.unitig_dir else None
    plans=[]; anchor_rows=[]
    # read oriented contexts from FNA
    fna={n.split('|',1)[1] if '|' in n else n:s for n,s,h in iter_fasta(b4/'DELTA12_NATIVE_CONTEXTS.fna')}
    # Actually headers are family|candidate_id; split only first pipe makes candidate id unavailable. rebuild map from full header suffix.
    fna={}
    for n,s,h in iter_fasta(b4/'DELTA12_NATIVE_CONTEXTS.fna'):
        # header starts family + | + full candidate id
        parts=n.split('|',1)
        if len(parts)==2:fna[parts[1]]=s
    for r in reps:
        acc=r['accession']; cid=r['candidate_id']; seq=fna.get(cid,'')
        up_anchor=seq[:K] if len(seq)>=K else ''
        dn_anchor=seq[-K:] if len(seq)>=K else ''
        url=f'https://s3.amazonaws.com/logan-pub/u/{acc}/{acc}.unitigs.fa.zst'
        size=status=0; err=''
        if args.network_head:size,status,err=head_size(url)
        local=''
        if unitdir:
            p=unitdir/f'{acc}.unitigs.fa.zst'
            if p.exists():local=str(p)
        plans.append({'family_rep_id':r['family_rep_id'],'candidate_id':cid,'accession':acc,'unitig_url':url,'remote_content_length_bytes':size,'head_status':status,'head_error':err,'local_unitig_file':local})
        if local and up_anchor and dn_anchor:
            hits=search_anchors(Path(local),{'UPSTREAM':up_anchor,'DOWNSTREAM':dn_anchor})
            for side in ['UPSTREAM','DOWNSTREAM']:
                xs=hits[side]
                classification='ABSENT'
                if len(xs)==1:
                    need=max(0,(STRICT_MIN_UPSTREAM-int(r['current_upstream_bp'])) if side=='UPSTREAM' else (PARTNER_MIN_DOWNSTREAM-int(r['current_downstream_bp'])))
                    ext=xs[0][3]
                    classification='UNIQUE_DIRECT_SUFFICIENT' if ext>=need else ('UNIQUE_DIRECT_PARTIAL' if ext>0 else 'UNIQUE_GRAPH_OR_NO_DIRECT_EXTENSION')
                elif len(xs)>1:classification='AMBIGUOUS'
                if not xs:
                    anchor_rows.append({'family_rep_id':r['family_rep_id'],'candidate_id':cid,'accession':acc,'side':side,'anchor':up_anchor if side=='UPSTREAM' else dn_anchor,'match_n':0,'classification':classification,'unitig_id':'','unitig_orientation':'','direct_extension_bp':'','unitig_length_bp':'','unitig_header':''})
                else:
                    for x in xs:
                        anchor_rows.append({'family_rep_id':r['family_rep_id'],'candidate_id':cid,'accession':acc,'side':side,'anchor':up_anchor if side=='UPSTREAM' else dn_anchor,'match_n':len(xs),'classification':classification,'unitig_id':x[0],'unitig_orientation':x[1],'direct_extension_bp':x[3],'unitig_length_bp':x[4],'unitig_header':x[5]})
    write_tsv(out/'DELTA12_UNITIG_DOWNLOAD_PREFLIGHT.tsv',plans,['family_rep_id','candidate_id','accession','unitig_url','remote_content_length_bytes','head_status','head_error','local_unitig_file'])
    if anchor_rows:write_tsv(out/'DELTA12_UNITIG_ANCHOR_MATCHES.tsv',anchor_rows,['family_rep_id','candidate_id','accession','side','anchor','match_n','classification','unitig_id','unitig_orientation','direct_extension_bp','unitig_length_bp','unitig_header'])
    # shell download plan
    with open(out/'DOWNLOAD_DELTA12_UNITIGS.sh','w') as f:
        f.write('#!/usr/bin/env bash\nset -euo pipefail\nOUT=${1:-unitigs}\nmkdir -p "$OUT"\n')
        for acc in sorted({r['accession'] for r in reps}):
            f.write(f'wget -c "https://s3.amazonaws.com/logan-pub/u/{acc}/{acc}.unitigs.fa.zst" -O "$OUT/{acc}.unitigs.fa.zst"\n')
    os.chmod(out/'DOWNLOAD_DELTA12_UNITIGS.sh',0o755)

    gate={
      'gate':'PASS_STAGE04B4R1_CONTEXT_AND_DESIGN_PREFLIGHT',
      'gate_pass':True,
      'representative_n':len(reps),
      'strict_array_evaluable_on_current_contig_n':sum(r['strict_array_evaluable_current']=='True' for r in reps),
      'partner_evaluable_on_current_contig_n':sum(r['partner_evaluable_current']=='True' for r in reps),
      'unitig_anchor_search_performed':bool(anchor_rows),
      'network_head_performed':bool(args.network_head),
      'candidate_membership_changed':False,'threshold_changed':False,'retuning_performed':False,'source_panel_reselection_performed':False,
      'claim_boundary':'This PASS means the audit completed. It does not mean the architecture endpoint is biologically negative or that deterministic extension succeeded.'
    }
    (out/'STAGE04B4R1_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate,indent=2))

if __name__=='__main__':main()
