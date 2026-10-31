#!/usr/bin/env python3
"""Post-freeze biology, never a mathematical detector rerun.
Uniform full-genome ORF + RT screen on every source; all candidates retained.
"""
from __future__ import annotations
import argparse, contextlib, csv, fcntl, gzip, hashlib, io, json, math, os
import platform, re, shutil, subprocess, sys, tarfile, time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, wait, FIRST_COMPLETED
from pathlib import Path
from common import *
from repeat_evidence import find_exact_core_groups, entropy
from approximate_repeats import find_approximate_groups
ROOT=Path(__file__).resolve().parents[1]
GENE_FIELDS=['record_key','sequence_id','prodigal_id','start0','end0','strand','translation_table','partial','start_type','nt_sha256','aa_sha256','aa_length','candidate_keys','protein_source']
DOMAIN_FIELDS=['target_id','target_length','profile_name','profile_accession','profile_role','hmm_length','full_evalue','full_bitscore','domain_index','domain_count','independent_domain_evalue','conditional_domain_evalue','domain_bitscore','hmm_from1','hmm_to1','ali_from1','ali_to1','hmm_coverage','source_scan','evidence_tier']
REPEAT_FIELDS=['candidate_key','sequence_id','repeat_group_id','assay','start0','end0','core_nt','core_length','n_copies','copy_intervals_json','copy_nt_json','mismatches_to_template_json','min_identity_to_template','n_seed_representations','template_entropy_bits','high_complexity','spacer_lengths_json','spacer_sequence_sha256_json','n_spacers','n_distinct_spacer_sequences','distinct_spacer_fraction','n_empty_spacers','spacing_cv','overlap_candidate_bp','predicted_CDS_overlap_span_bp','predicted_CDS_overlap_copies_bp','all_copy_bases_outside_predicted_CDS','full_span_outside_predicted_CDS','official_CDS_assessed','ncRNA_proven','new_ART_proven']
LOCUS_EXTRA=['target_id','rt_core_GA','rt_core_sensitivity_only','rt_aux_GA','rt_aux_sensitivity_only','profile_evidence','exact_reference_AA_ids','mathematical_candidate_neighborhood_captured','RT_activity_proven']

def require(v,msg):
    if not v:raise RuntimeError(msg)
def policy():return read_json(ROOT/'config/POLICY.json')
def safe_path(base,rel):
    base=Path(base).resolve();p=(base/rel).resolve()
    require(p.is_relative_to(base),'Path outside source root: '+str(rel));return p
def verify_file(p,h):
    require(Path(p).is_file() and sha_file(p)==h,'Missing/changed frozen file: '+str(p))
def package_check():
    for r in read_tsv(ROOT/'PACKAGE_MANIFEST.tsv'):
        p=safe_path(ROOT,r['path']);verify_file(p,r['sha256']);require(p.stat().st_size==int(r['size_bytes']),'Size changed '+str(p))

def tools():
    result={}
    for name,args,pattern in [('prodigal',['-v'],r'V?2\.6\.3'),('hmmsearch',['-h'],r'HMMER\s+3\.4\b')]:
        p=shutil.which(name);require(p is not None,'Missing '+name+'; activate existing art_rt_evidence environment')
        cp=subprocess.run([p]+args,capture_output=True,text=True);t=cp.stdout+cp.stderr
        require(bool(re.search(pattern,t)),'Unexpected version '+name+': '+t[:200])
        result[name]={'path':str(Path(p).resolve()),'sha256':sha_file(p),'version_text':'\n'.join(t.splitlines()[:5])}
    return result

def input_check(a05,b05,check_archive=True):
    p=policy();a05=Path(a05).resolve();b05=Path(b05).resolve()
    for rel,expected in [('config/FROZEN_NEW_GENOME_MANIFEST.tsv',p['source_manifest_sha256']),('outputs/STAGE05B_DUAL_ARM_CANDIDATES.tsv',p['source_candidate_sha256']),('outputs/STAGE05B_ALL_TWO_ARM_UNION_REGIONS.tsv',p['source_all_regions_sha256']),('outputs/SCORE_FILE_MANIFEST.tsv',p['source_score_manifest_sha256'])]:verify_file(b05/rel,expected)
    verify_file(a05/'outputs/FROZEN_NEW_GENOME_MANIFEST.tsv',p['source_manifest_sha256'])
    verify_file(b05/'outputs/STAGE05B_FINAL_GATE.json',sha_file(ROOT/'config/STAGE05B_FINAL_GATE.json'))
    verify_file(b05/'outputs/RUN_LOCK.json',sha_file(ROOT/'config/SOURCE_RUN_LOCK.json'))
    if check_archive:verify_file(b05/p['source_score_archive_name'],p['source_score_archive_sha256'])
    gate=read_json(b05/'outputs/STAGE05B_FINAL_GATE.json');require(gate['gate_pass'] is True and gate['n_success_genomes']==p['n_genomes'],'Source score board incomplete')
    sources=read_tsv(ROOT/'config/FROZEN_NEW_GENOME_MANIFEST.tsv');cands=read_tsv(ROOT/'config/STAGE05B_DUAL_ARM_CANDIDATES.tsv')
    require(len(sources)==p['n_genomes'] and len({s['accession'] for s in sources})==len(sources),'Source membership mismatch')
    require(len(cands)==p['n_candidates'] and len({(c['sequence_id'],c['region_id']) for c in cands})==len(cands),'Candidate membership mismatch')
    require(len({c['sequence_id'] for c in cands})==p['n_candidate_genomes'],'Candidate-genome count changed')
    sm={s['accession']:s for s in sources};checkrows=[];oldlock=read_json(ROOT/'config/SOURCE_RUN_LOCK.json')
    expected_score_files={r['path']:r for r in read_tsv(ROOT/'config/SCORE_FILE_MANIFEST.tsv')}
    for s in sources:
        acc=s['accession'];f=safe_path(a05,s['fasta_path']);verify_file(f,s['fasta_sha256']);rr=list(fasta(f))
        require(len(rr)==1 and rr[0][0]==acc and len(rr[0][2])==int(s['length_bp']) and not set(rr[0][2])-set('ACGT'),'FASTA identity/alphabet changed '+acc)
        seq=rr[0][2];require(sha_text(seq)==s['sequence_sha256'] and sha_text(rc(seq))==s['reverse_complement_sha256'],'Sequence hash changed '+acc)
        comp=read_json(b05/'outputs/genomes'/acc/'COMPLETE.json');require(comp['signature']['fasta_sha256']==s['fasta_sha256'] and comp['signature']['runtime']==oldlock['runtime'] and comp['signature']['method_sha256']==oldlock['method_sha256'],'Scoring source signature changed '+acc)
        for filename,h in comp['files'].items():
            rel='outputs/genomes/'+acc+'/'+filename
            require(rel in expected_score_files and expected_score_files[rel]['sha256']==h,'Completion/archived score manifest disagreement '+rel)
            verify_file(b05/rel,h)
        if (len(checkrows)+1)%50==0:print('[INPUT RECHECK]',len(checkrows)+1,'/',len(sources),flush=True)
        checkrows.append(dict(accession=acc,length_bp=len(seq),fasta_sha256=s['fasta_sha256'],sequence_sha256=s['sequence_sha256'],score_files_rehashed=True))
    for c in cands:
        n=int(sm[c['sequence_id']]['length_bp']);require(0<=int(c['start0'])<int(c['end0'])<=n and c['dual_arm_candidate']=='1','Invalid candidate interval')
    require(sum(int(s['length_bp']) for s in sources)==p['total_bp'],'Source bp changed')
    return sources,cands,checkrows

def preflight(a05,b05,cpus):
    require(cpus>0,'cpus must be positive');package_check();ts=tools()
    cpu=len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else (os.cpu_count() or 1)
    require(cpus+1<=cpu,'Requested CPU workers exceed available affinity')
    require(shutil.disk_usage(ROOT).free>=5*2**30,'Need at least 5 GiB free for launch (not a guaranteed final disk bound)')
    sources,cands,checks=input_check(a05,b05)
    lock={'stage':'Stage05C1','policy_sha256':sha_file(ROOT/'config/POLICY.json'),'package_manifest_sha256':sha_file(ROOT/'PACKAGE_MANIFEST.tsv'),'source_root':str(Path(a05).resolve()),'score_root':str(Path(b05).resolve()),'tools':ts,'python':platform.python_version(),'python_executable':sys.executable,'platform':platform.platform(),'candidate_sha256':policy()['source_candidate_sha256'],'score_archive_sha256':policy()['source_score_archive_sha256']}
    dest=ROOT/'outputs/RUN_LOCK.json'
    if dest.exists():require(read_json(dest)==lock,'Execution lock changed. Use same sources/tools/interpreter; do not mix results.')
    else:write_json(dest,lock)
    write_tsv(ROOT/'outputs/SOURCE919_BIOLOGY_INPUT_RECHECK.tsv',checks,list(checks[0]))
    print('[INPUT PASS]',len(sources),'genomes;',len(cands),'candidates; score archive and raw windows rehashed',flush=True)
    return sources,cands,sha_file(dest)

def parse_pred_header(desc):
    parts=[x.strip() for x in desc.split('#')];require(len(parts)>=5,'Unsupported Prodigal FASTA header')
    ident=parts[0].split()[0];a,b,sign=int(parts[1])-1,int(parts[2]),int(parts[3]);require(0<=a<b and sign in (-1,1),'Invalid predicted coordinates')
    attrs=dict(x.split('=',1) for x in parts[4].split(';') if '=' in x)
    return ident,a,b,'+' if sign==1 else '-',attrs

def parse_prediction(folder,acc,seq,candidates):
    folder=Path(folder);aa=fasta_dict(folder/'proteins.faa');dna=fasta_dict(folder/'cds.fna');gff=(folder/'genes.gff').read_text()
    tt=set(re.findall(r'transl_table=(\d+)',gff));require(len(tt)==1,'Missing/ambiguous translation table '+acc)
    require(set(aa)==set(dna),'Prodigal amino-acid/CDS identifiers mismatch '+acc)
    gl={}
    for line in gff.splitlines():
        if not line or line.startswith('#'):continue
        f=line.split('\t');require(len(f)==9 and f[0]==acc and f[2]=='CDS','Unexpected GFF entry '+acc)
        at=dict(x.split('=',1) for x in f[8].split(';') if '=' in x);gid=at['ID'];require(gid not in gl,'Duplicate GFF feature ID '+acc)
        gl[gid]=(int(f[3])-1,int(f[4]),f[6])
    rows=[];proteins={};r=policy()['neighborhood_bp']
    for i,(pid,(desc,protein)) in enumerate(aa.items(),1):
        _,a,b,sign,at=parse_pred_header(desc);_,na,nb,ns,ntat=parse_pred_header(dna[pid][0]);gid=at.get('ID','')
        require((a,b,sign)==(na,nb,ns)==gl.get(gid) and gid==ntat.get('ID') and b<=len(seq),'Predicted coordinate/ID disagreement '+acc+'/'+pid)
        nt=seq[a:b] if sign=='+' else rc(seq[a:b]);require(nt==dna[pid][1],'CDS fails exact directional genomic reconstruction '+acc+'/'+pid)
        protein=protein[:-1] if protein.endswith('*') else protein
        require(bool(protein) and not set(protein)-set('ACDEFGHIKLMNPQRSTVWY'),'Empty/noncanonical protein; no stop-to-X replacement permitted '+acc+'/'+pid)
        key=acc+f'__PRODIGAL_{i:05d}';cks=[acc+'::'+c['region_id'] for c in candidates if overlaps((a,b),(max(0,int(c['start0'])-r),min(len(seq),int(c['end0'])+r)))]
        rows.append(dict(record_key=key,sequence_id=acc,prodigal_id=pid,start0=a,end0=b,strand=sign,translation_table=int(next(iter(tt))),partial=at.get('partial',''),start_type=at.get('start_type',''),nt_sha256=sha_text(nt),aa_sha256=sha_text(protein),aa_length=len(protein),candidate_keys=';'.join(cks),protein_source='PRODIGAL_META_NEW919_NOT_OFFICIAL'))
        proteins[key]=protein
    require(len(rows)==len(gl),'Missing GFF/protein counterpart '+acc)
    return rows,proteins

def checkpoint(folder,fingerprint,files):
    marker=Path(folder)/'COMPLETE.json'
    if not marker.exists():return False
    j=read_json(marker);require(j['fingerprint']==fingerprint,'Checkpoint signature changed: '+str(folder))
    for f in files:verify_file(Path(folder)/f,j['files'][f])
    return True

def predict_one(s,cands,a05,lock):
    acc=s['accession'];f=safe_path(a05,s['fasta_path']);verify_file(f,s['fasta_sha256']);seq=next(fasta(f))[2]
    d=ROOT/'work/prodigal'/acc;d.mkdir(parents=True,exist_ok=True);raw=['proteins.faa','cds.fna','genes.gff'];fp={'run_lock':lock,'fasta_sha256':s['fasta_sha256'],'mode':'meta'}
    partial=d/'JOB_SIGNATURE.json'
    if partial.exists():require(read_json(partial)==fp,'Different partial ORF job '+acc)
    else:write_json(partial,fp)
    if not checkpoint(d,fp,raw):
        run_cmd(['prodigal','-i',f,'-p','meta','-f','gff','-o',d/'genes.gff','-a',d/'proteins.faa','-d',d/'cds.fna','-q'],ROOT/'logs'/('prodigal_'+acc+'.log'))
    rr,aa=parse_prediction(d,acc,seq,cands)
    verify_file(f,s['fasta_sha256']);write_json(d/'COMPLETE.json',{'fingerprint':fp,'files':{n:sha_file(d/n) for n in raw},'n_genes':len(rr)})
    print('[ORF FROZEN]',acc,'genes=',len(rr),flush=True)
    return acc,rr,aa

def limited_jobs(items,workers,fn):
    itr=iter(items);futures={};results=[]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for _ in range(workers):
            x=next(itr,None)
            if x is not None:futures[pool.submit(fn,x)]=x
        while futures:
            done,_=wait(list(futures),return_when=FIRST_COMPLETED)
            for f in done:
                futures.pop(f);results.append(f.result())
                x=next(itr,None)
                if x is not None:futures[pool.submit(fn,x)]=x
    return results

def predict_all(sources,cands,a05,cpus,lock):
    by=defaultdict(list)
    for c in cands:by[c['sequence_id']].append(c)
    def one(s):return predict_one(s,by[s['accession']],a05,lock)
    res=limited_jobs(sources,cpus,one);res.sort(key=lambda x:x[0]);rr=[];aas={};cov=[]
    for acc,rows,seqs in res:
        rr.extend(rows);aas.update(seqs);cov.append({'accession':acc,'predicted_genes':len(rows),'candidate_neighborhood_genes':sum(bool(r['candidate_keys']) for r in rows),'n_candidates':len(by[acc])})
    raw_files=[]
    for acc,_,_ in res:
        for name in ['proteins.faa','cds.fna','genes.gff','COMPLETE.json']:
            f=ROOT/'work/prodigal'/acc/name;raw_files.append(dict(sequence_id=acc,path=str(f.relative_to(ROOT)),size_bytes=f.stat().st_size,sha256=sha_file(f)))
    write_tsv(ROOT/'outputs/ALL919_PRODIGAL_FILE_MANIFEST.tsv',raw_files,['sequence_id','path','size_bytes','sha256'])
    write_tsv(ROOT/'outputs/FULL919_GENE_REGISTRY.tsv',rr,GENE_FIELDS);write_tsv(ROOT/'outputs/ORF_GENOME_COVERAGE.tsv',cov,list(cov[0]))
    write_fasta(ROOT/'outputs/FULL919_PREDICTED_PROTEINS.faa',[(r['record_key'],aas[r['record_key']]) for r in rr])
    return rr,aas

def build_targets(genes,aas,lock):
    targets={};fan=[]
    for r in genes:
        s=aas[r['record_key']];require(sha_text(s)==r['aa_sha256'],'Protein hash mismatch');tid='AA_'+r['aa_sha256'];targets[tid]=s
        fan.append(dict(target_id=tid,record_key=r['record_key'],sequence_id=r['sequence_id'],candidate_keys=r['candidate_keys'],source_branch='NEW919_PREDICTED'))
    refs=fasta_dict(ROOT/'inputs/controls/RT_REFERENCE_PANEL.faa')
    for rid,(_,s) in refs.items():
        tid='AA_'+sha_text(s);targets[tid]=s;fan.append(dict(target_id=tid,record_key=rid,sequence_id='',candidate_keys='',source_branch='REGISTERED_REFERENCE_NOT_NEW_DISCOVERY'))
    stage=ROOT/'work/targets';stage.mkdir(parents=True,exist_ok=True);fa=stage/'RT_SCAN_TARGETS.faa';ft=stage/'RT_TARGET_FANOUT.tsv'
    write_fasta(fa,sorted(targets.items()));write_tsv(ft,fan,list(fan[0]))
    targetlock={'run_lock':lock,'targets_sha256':sha_file(fa),'fanout_sha256':sha_file(ft),'n_unique_scan_targets_including_references':len(targets),'n_full919_predicted_records':len(genes),'n_references':len(refs),'n_full919_unique_AA':len({r['aa_sha256'] for r in genes}),'n_neighborhood_records':sum(bool(r['candidate_keys']) for r in genes),'n_neighborhood_unique_AA':len({r['aa_sha256'] for r in genes if r['candidate_keys']})}
    lp=ROOT/'outputs/TARGET_UNIVERSE_LOCK.json'
    if lp.exists():
        require(read_json(lp)==targetlock,'Targets changed after first scan freeze')
        verify_file(ROOT/'outputs'/fa.name,targetlock['targets_sha256']);verify_file(ROOT/'outputs'/ft.name,targetlock['fanout_sha256'])
    else:
        shutil.copyfile(fa,ROOT/'outputs'/fa.name);shutil.copyfile(ft,ROOT/'outputs'/ft.name);write_json(lp,targetlock)
    return targets,fan,refs,targetlock

def parse_domains(path,targets,catalog,mode):
    lookup={m['name']:m for m in catalog};s=policy()['rt_scan']['sensitivity'];rows=[]
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith('#'):continue
        x=line.split(maxsplit=22);require(len(x)>=22,'Truncated HMMER domain table')
        require(x[0] in targets and x[3] in lookup,'Unknown scan target/model')
        m=lookup[x[3]];L,HL=int(x[2]),int(x[5]);hf,ht,af,at=map(int,[x[15],x[16],x[17],x[18]])
        require(L==len(targets[x[0]]) and HL==int(m['length']) and x[4]==m['accession_version'],'Domain source mismatch')
        require(1<=hf<=ht<=HL and 1<=af<=at<=L,'Domain interval invalid')
        ev,bs,ie,ce,ds=map(float,[x[6],x[7],x[12],x[11],x[13]]);cov=(ht-hf+1)/HL
        require(all(math.isfinite(z) for z in [ev,bs,ie,ce,ds,cov]),'Nonfinite domain values')
        tier='REPORTED_ONLY';prefix='RT_CORE' if m['role']=='RT_CORE' else 'RT_LIKE_AUXILIARY'
        if mode=='GA':
            require(bs>=float(m['ga_sequence_bits'])-1e-8 and ds>=float(m['ga_domain_bits'])-1e-8,'GA output fails declared model bitscore thresholds')
            tier=prefix+'_GA'
        elif ev<=s['full_evalue_max'] and ie<=s['independent_domain_evalue_max'] and cov>=s['hmm_coverage_min']:tier=prefix+'_SENSITIVITY'
        rows.append(dict(target_id=x[0],target_length=L,profile_name=x[3],profile_accession=x[4],profile_role=m['role'],hmm_length=HL,full_evalue=ev,full_bitscore=bs,domain_index=int(x[9]),domain_count=int(x[10]),independent_domain_evalue=ie,conditional_domain_evalue=ce,domain_bitscore=ds,hmm_from1=hf,hmm_to1=ht,ali_from1=af,ali_to1=at,hmm_coverage=cov,source_scan=mode,evidence_tier=tier))
    return rows

def hmm_scan(targets,cpus,lock):
    catalog=read_tsv(ROOT/'refs/MODEL_CATALOG.tsv');m=ROOT/'refs/RT_PANEL.hmm';t=ROOT/'outputs/RT_SCAN_TARGETS.faa';d=ROOT/'outputs/hmmer';d.mkdir(parents=True,exist_ok=True);allrows=[]
    for mode in ['GA','SENSITIVITY']:
        job=d/mode;job.mkdir(exist_ok=True);files=['report.txt','hits.domtbl','hits.tbl'];fp={'run_lock':lock,'target_sha256':sha_file(t),'model_sha256':sha_file(m),'mode':mode,'Z':policy()['rt_scan']['database_size_Z']}
        if not checkpoint(job,fp,files):
            cmd=['hmmsearch','--cpu',str(cpus),'--seed','42','-Z',str(fp['Z']),'--noali','-o',job/'report.txt','--domtblout',job/'hits.domtbl','--tblout',job/'hits.tbl']
            cmd+=['--cut_ga'] if mode=='GA' else ['-E','10','--domE','10'];cmd+=[m,t]
            print('[HMMER]',mode,'unique_targets=',len(targets),'Z=',fp['Z'],flush=True);run_cmd(cmd,ROOT/'logs'/('hmmer_'+mode+'.log'))
        txt=(job/'report.txt').read_text();require('[ok]' in txt and sorted(re.findall(r'^Query:\s+(\S+)\s+\[M=\d+\]',txt,re.M))==sorted(c['name'] for c in catalog),'Incomplete HMMER execution')
        rows=parse_domains(job/'hits.domtbl',targets,catalog,mode);allrows.extend(rows);write_json(job/'COMPLETE.json',{'fingerprint':fp,'files':{f:sha_file(job/f) for f in files},'domain_rows':len(rows)})
    allrows.sort(key=lambda r:(r['target_id'],r['source_scan'],r['profile_accession'],r['ali_from1']))
    write_tsv(ROOT/'outputs/RT_ALL_REPORTED_DOMAINS.tsv',allrows,DOMAIN_FIELDS);return allrows

def statuses(domains):
    d=defaultdict(set)
    for r in domains:
        if r['evidence_tier']!='REPORTED_ONLY':d[r['target_id']].add(r['evidence_tier'])
    return d

def rt_loci(genes,refs,domains):
    st=statuses(domains);refsha=defaultdict(list)
    for k,(_,s) in refs.items():refsha[sha_text(s)].append(k)
    rr=[]
    for g in genes:
        tid='AA_'+g['aa_sha256'];flags=st[tid];core=int('RT_CORE_GA' in flags);sens=int('RT_CORE_SENSITIVITY' in flags and not core);aux=int('RT_LIKE_AUXILIARY_GA' in flags);auxs=int('RT_LIKE_AUXILIARY_SENSITIVITY' in flags and not aux)
        if not any([core,sens,aux,auxs]):continue
        rr.append(dict(g,target_id=tid,rt_core_GA=core,rt_core_sensitivity_only=sens,rt_aux_GA=aux,rt_aux_sensitivity_only=auxs,profile_evidence=';'.join(sorted(flags)),exact_reference_AA_ids=';'.join(refsha[g['aa_sha256']]),mathematical_candidate_neighborhood_captured=int(bool(g['candidate_keys'])),RT_activity_proven=False))
    write_tsv(ROOT/'outputs/ALL919_RT_PROFILE_LOCI.tsv',rr,GENE_FIELDS+LOCUS_EXTRA)
    readrefs=read_tsv(ROOT/'inputs/controls/RT_REFERENCE_PANEL.tsv');rrrefs=[]
    for r in readrefs:
        k=r['reference_id'];flag=st['AA_'+sha_text(refs[k][1])]
        rrrefs.append(dict(reference_id=k,reference_class=r.get('reference_class',''),rt_core_GA=int('RT_CORE_GA' in flag),profile_status=';'.join(sorted(flag)),specificity_test=False,novel_discovery=False))
    write_tsv(ROOT/'outputs/REFERENCE_RT_CALIBRATION.tsv',rrrefs,list(rrrefs[0]));return rr,rrrefs

WITNESS_FIELDS=['candidate_key','sequence_id','window_start0','window_end0','window_length','k','lag_bp','n_compared_kmer_pairs','n_matching_kmer_pairs','n_overlapping_pair_blocks','blocks_json','window_predicted_CDS_overlap_bp','matching_kmer_count_is_repeat_copy_count','detector_score_recomputed']
CONTEXT_FIELDS=['candidate_key','sequence_id','candidate_start0','candidate_end0','neighborhood_start0','neighborhood_end0','neighborhood_length','left_clipped','right_clipped','sequence_sha256','n_predicted_neighborhood_genes','official_annotation_used','circular_wrap_used']

def recurrence_witness(seq,c,mask):
    acc=c['sequence_id'];a,b=int(c['recurrence_best_window_start0']),int(c['recurrence_best_window_end0']);k=int(c['recurrence_best_k']);lag=int(c['recurrence_best_lag_bp'])
    require(0<=a<b<=len(seq) and k>0 and lag>=k and b-a-k+1-lag>=32,'Frozen recurrence window invalid')
    z=seq[a:b];npairs=len(z)-k+1-lag;positions=[i for i in range(npairs) if z[i:i+k]==z[i+lag:i+lag+k]]
    groups=[]
    for i in positions:
        if groups and i<groups[-1][-1]+k:groups[-1].append(i)
        else:groups.append([i])
    blocks=[]
    for ps in groups:
        st,en=ps[0],ps[-1]+k;first=[a+st,a+en];second=[a+st+lag,a+en+lag]
        require(seq[first[0]:first[1]]==seq[second[0]:second[1]],'Paired block sequence reconstruction failed')
        blocks.append({'first':first,'second':second,'sequence':seq[first[0]:first[1]],'block_length':en-st,'n_overlapping_kmer_matches':len(ps),'two_blocks_overlap':overlaps(first,second),'first_predicted_CDS_overlap_bp':covered_bp(first,mask),'second_predicted_CDS_overlap_bp':covered_bp(second,mask)})
    return dict(candidate_key=acc+'::'+c['region_id'],sequence_id=acc,window_start0=a,window_end0=b,window_length=b-a,k=k,lag_bp=lag,n_compared_kmer_pairs=npairs,n_matching_kmer_pairs=len(positions),n_overlapping_pair_blocks=len(blocks),blocks_json=canonical_json(blocks),window_predicted_CDS_overlap_bp=covered_bp((a,b),mask),matching_kmer_count_is_repeat_copy_count=False,detector_score_recomputed=False)

def annotate_repeat(seq,c,n0,group,assay,ordinal,mask):
    L=group['repeat_len'];parts=[[n0+s,n0+s+L] for s in group['starts']];span=(parts[0][0],parts[-1][1]);core=group['core'];copyseq=[seq[a:b] for a,b in parts]
    require(len(parts)>=3 and all(parts[i][1]<=parts[i+1][0] for i in range(len(parts)-1)),'Overlapping/insufficient copies')
    mms=[sum(x!=y for x,y in zip(core,s)) for s in copyseq]
    require(all(len(s)==L for s in copyseq),'Repeat crosses source boundary')
    if assay=='EXACT_CORE_PRIMARY_DESCRIPTION':require(max(mms)==0,'Nonexact core emitted')
    else:require(max(mms)<=math.floor(L*policy()['approximate_repeat_sensitivity']['max_substitution_fraction']+1e-12),'Substitution limit exceeded')
    steps=[b[0]-a[0] for a,b in zip(parts,parts[1:])];mean=sum(steps)/len(steps);cv=math.sqrt(sum((v-mean)**2 for v in steps)/len(steps))/mean
    spacers=[seq[a[1]:b[0]] for a,b in zip(parts,parts[1:])]
    pc=covered_bp(span,mask);cc=sum(covered_bp(p,mask) for p in parts);ent=entropy(core);key=c['sequence_id']+'::'+c['region_id'];short='EX' if assay.startswith('EXACT') else 'AP'
    return dict(candidate_key=key,sequence_id=c['sequence_id'],repeat_group_id=key+f'::{short}{ordinal:05d}',assay=assay,start0=span[0],end0=span[1],core_nt=core,core_length=L,n_copies=len(parts),copy_intervals_json=canonical_json(parts),copy_nt_json=canonical_json(copyseq),mismatches_to_template_json=canonical_json(mms),min_identity_to_template=1-max(mms)/L,n_seed_representations=group['seed_count'],template_entropy_bits=ent,high_complexity=int(ent>=1.5),spacer_lengths_json=canonical_json([len(z) for z in spacers]),spacer_sequence_sha256_json=canonical_json([sha_text(z) for z in spacers]),n_spacers=len(spacers),n_distinct_spacer_sequences=len(set(spacers)),distinct_spacer_fraction=len(set(spacers))/len(spacers),n_empty_spacers=sum(not z for z in spacers),spacing_cv=cv,overlap_candidate_bp=covered_bp(span,[(int(c['start0']),int(c['end0']))]),predicted_CDS_overlap_span_bp=pc,predicted_CDS_overlap_copies_bp=cc,all_copy_bases_outside_predicted_CDS=int(cc==0),full_span_outside_predicted_CDS=int(pc==0),official_CDS_assessed=False,ncRNA_proven=False,new_ART_proven=False)

def repeat_genome_worker(item):
    s,cands,genes,a05,lock=item;acc=s['accession'];path=safe_path(a05,s['fasta_path']);verify_file(path,s['fasta_sha256']);seq=next(fasta(path))[2];mask=merge_intervals([(int(g['start0']),int(g['end0'])) for g in genes]);rr=[];ww=[];contexts=[]
    gh=sha_text(canonical_json(genes));r=policy()['neighborhood_bp']
    for c in cands:
        key=acc+'::'+c['region_id'];a,b=int(c['start0']),int(c['end0']);n0,n1=max(0,a-r),min(len(seq),b+r);d=ROOT/'work/repeats'/acc/c['region_id'];d.mkdir(parents=True,exist_ok=True)
        require(shutil.disk_usage(ROOT).free>=2*2**30,'Disk headroom below 2 GiB; stop without dropping candidates')
        fp={'run_lock':lock,'fasta_sha256':s['fasta_sha256'],'candidate':c,'gene_registry_sha256':gh};files=['groups.tsv','witness.json','context.json']
        if checkpoint(d,fp,files):
            cached=read_tsv(d/'groups.tsv')
            for x in cached:
                for k in ['start0','end0','core_length','n_copies','n_seed_representations','n_spacers','n_distinct_spacer_sequences','n_empty_spacers','high_complexity','overlap_candidate_bp','predicted_CDS_overlap_span_bp','predicted_CDS_overlap_copies_bp','all_copy_bases_outside_predicted_CDS','full_span_outside_predicted_CDS']:x[k]=int(x[k])
                for k in ['min_identity_to_template','template_entropy_bits','spacing_cv','distinct_spacer_fraction']:x[k]=float(x[k])
                for k in ['official_CDS_assessed','ncRNA_proven','new_ART_proven']:x[k]=x[k]=='True'
            rr.extend(cached);ww.append(read_json(d/'witness.json'));contexts.append(read_json(d/'context.json'));continue
        prior=d/'JOB_SIGNATURE.json'
        if prior.exists():require(read_json(prior)==fp,'Changed partial repeat job '+key)
        else:write_json(prior,fp)
        local=[]
        for assay,fn,pol in [('EXACT_CORE_PRIMARY_DESCRIPTION',find_exact_core_groups,policy()['repeat_survey']),('UNGAPPED_SUBSTITUTION_SENSITIVITY',find_approximate_groups,policy()['approximate_repeat_sensitivity'])]:
            found=fn(seq[n0:n1],pol)
            for i,g in enumerate(found,1):local.append(annotate_repeat(seq,c,n0,g,assay,i,mask))
        w=recurrence_witness(seq,c,mask);ct=dict(candidate_key=key,sequence_id=acc,candidate_start0=a,candidate_end0=b,neighborhood_start0=n0,neighborhood_end0=n1,neighborhood_length=n1-n0,left_clipped=int(a-r<0),right_clipped=int(b+r>len(seq)),sequence_sha256=sha_text(seq[n0:n1]),n_predicted_neighborhood_genes=sum(key in g['candidate_keys'].split(';') for g in genes),official_annotation_used=False,circular_wrap_used=False)
        write_tsv(d/'groups.tsv',local,REPEAT_FIELDS);write_json(d/'witness.json',w);write_json(d/'context.json',ct);write_json(d/'COMPLETE.json',{'fingerprint':fp,'files':{f:sha_file(d/f) for f in files},'n_repeat_groups':len(local)})
        rr.extend(local);ww.append(w);contexts.append(ct)
    verify_file(path,s['fasta_sha256']);return acc,rr,ww,contexts

def repeat_all(sources,cands,genes,a05,cpus,lock):
    cg=defaultdict(list);gg=defaultdict(list)
    for c in cands:cg[c['sequence_id']].append(c)
    for g in genes:gg[g['sequence_id']].append(g)
    jobs=[(s,cg[s['accession']],gg[s['accession']],str(a05),lock) for s in sources if cg[s['accession']]];allr=[];allw=[];allc=[];donecount=0
    with ProcessPoolExecutor(max_workers=cpus) as pool:
        it=iter(jobs);pending={}
        for _ in range(cpus):
            x=next(it,None)
            if x:pending[pool.submit(repeat_genome_worker,x)]=x[0]['accession']
        while pending:
            done,_=wait(list(pending),return_when=FIRST_COMPLETED)
            for f in done:
                pending.pop(f);acc,r,w,c=f.result();allr.extend(r);allw.extend(w);allc.extend(c);donecount+=len(c)
                print('[REPEAT FROZEN]',donecount,'/',len(cands),acc,'groups=',len(r),flush=True)
                x=next(it,None)
                if x:pending[pool.submit(repeat_genome_worker,x)]=x[0]['accession']
    order={c['sequence_id']+'::'+c['region_id']:i for i,c in enumerate(cands)}
    allr.sort(key=lambda r:(order[r['candidate_key']],r['repeat_group_id']));allw.sort(key=lambda r:order[r['candidate_key']]);allc.sort(key=lambda r:order[r['candidate_key']])
    require(len(allw)==len(allc)==len(cands) and {w['candidate_key'] for w in allw}==set(order),'Repeat evidence coverage incomplete')
    write_tsv(ROOT/'outputs/REPEAT_GROUPS_ALL_CANDIDATES.tsv',allr,REPEAT_FIELDS);write_tsv(ROOT/'outputs/RECURRENCE_MATCH_BLOCK_WITNESSES.tsv',allw,WITNESS_FIELDS);write_tsv(ROOT/'outputs/CANDIDATE_NEIGHBORHOOD_REGISTRY.tsv',allc,CONTEXT_FIELDS)
    return allr,allw,allc

def reference_repeat_calibration():
    rows=[];group_rows=[]
    for r in read_tsv(ROOT/'inputs/controls/REPEAT_REFERENCE_CONTEXT.tsv'):
        f=ROOT/r['path'];seq=next(fasta(f))[2];require(sha_text(seq)==r['context_sha256'],'Reference nucleotide-content hash mismatch');offset=int(r['context_start0']);core=(int(r['registered_core_start0'])-offset,int(r['registered_core_end0'])-offset)
        require(len(seq)==int(r['context_end0'])-offset and 0<=core[0]<core[1]<=len(seq),'Reference context coordinates inconsistent')
        result={'reference_id':r['reference_id'],'role':'REGISTERED_KNOWN_POSITIVE_DESCRIPTIVE_CAPABILITY_NOT_INDEPENDENT_VALIDATION'}
        for label,fn,pol in [('exact',find_exact_core_groups,policy()['repeat_survey']),('substitution',find_approximate_groups,policy()['approximate_repeat_sensitivity'])]:
            found=fn(seq,pol);hits=0
            for i,g in enumerate(found,1):
                copies=[[s,s+g['repeat_len']] for s in g['starts']];nc=sum(overlaps(x,core) for x in copies);hc=entropy(g['core'])>=1.5;hits+=int(nc>=3 and hc)
                group_rows.append(dict(reference_id=r['reference_id'],assay=label,group_index=i,core_nt=g['core'],core_length=g['repeat_len'],n_copies=len(copies),genomic_copy_intervals_json=canonical_json([[a+offset,b+offset] for a,b in copies]),n_copies_overlapping_registered_core=nc,high_complexity=int(hc)))
            result[label+'_context_groups']=len(found);result[label+'_high_complexity_3copy_core_groups']=hits;result[label+'_registered_core_recovered']=int(hits>0)
        result['parameters_selected_on_new919']=False;result['new_ART_proven']=False;rows.append(result)
    write_tsv(ROOT/'outputs/REFERENCE_REPEAT_CALIBRATION.tsv',rows,list(rows[0]));write_tsv(ROOT/'outputs/REFERENCE_REPEAT_GROUPS.tsv',group_rows,list(group_rows[0]) if group_rows else ['reference_id','assay','group_index'])
    return rows

PAIR_FIELDS=['candidate_key','sequence_id','rt_record_key','rt_start0','rt_end0','rt_strand','rt_evidence','repeat_group_id','repeat_assay','repeat_start0','repeat_end0','repeat_copies','repeat_min_identity_to_template','repeat_high_complexity','repeat_full_span_outside_predicted_CDS','repeat_overlaps_math_candidate','distance_RT_to_repeat_bp','repeat_relative_to_RT_coding_direction','exact_reference_AA_ids','review_evidence_only_not_system_confirmation']

def pair_context(c,locus,rep):
    a,b=int(locus['start0']),int(locus['end0']);x,y=int(rep['start0']),int(rep['end0']);side='OVERLAPS_RT' if overlaps((a,b),(x,y)) else ('LEFT_OF_RT' if y<=a else 'RIGHT_OF_RT')
    relative='OVERLAPS_RT' if side=='OVERLAPS_RT' else ('UPSTREAM_RELATIVE_TO_PREDICTED_RT' if (side=='LEFT_OF_RT')==(locus['strand']=='+') else 'DOWNSTREAM_RELATIVE_TO_PREDICTED_RT')
    return dict(candidate_key=c['sequence_id']+'::'+c['region_id'],sequence_id=c['sequence_id'],rt_record_key=locus['record_key'],rt_start0=a,rt_end0=b,rt_strand=locus['strand'],rt_evidence=locus['profile_evidence'],repeat_group_id=rep['repeat_group_id'],repeat_assay=rep['assay'],repeat_start0=x,repeat_end0=y,repeat_copies=rep['n_copies'],repeat_min_identity_to_template=rep['min_identity_to_template'],repeat_high_complexity=rep['high_complexity'],repeat_full_span_outside_predicted_CDS=rep['full_span_outside_predicted_CDS'],repeat_overlaps_math_candidate=int(rep['overlap_candidate_bp']>0),distance_RT_to_repeat_bp=distance((a,b),(x,y)),repeat_relative_to_RT_coding_direction=relative,exact_reference_AA_ids=locus['exact_reference_AA_ids'],review_evidence_only_not_system_confirmation=True)

def integrate(sources,cands,genes,loci,reps,witnesses,contexts):
    lg=defaultdict(list);rg=defaultdict(list);gg=defaultdict(list);gacc=defaultdict(list);wm={w['candidate_key']:w for w in witnesses};cm={r['candidate_key']:r for r in contexts}
    for l in loci:
        for k in l['candidate_keys'].split(';'):
            if k:lg[k].append(l)
    for r in reps:rg[r['candidate_key']].append(r)
    for g in genes:
        gacc[g['sequence_id']].append(g)
        for k in g['candidate_keys'].split(';'):
            if k:gg[k].append(g)
    summary=[];pairs=[];cards=[];links=[]
    for c in cands:
        k=c['sequence_id']+'::'+c['region_id'];ls=lg[k];rs=rg[k];gs=gg[k];ct=cm[k];w=wm[k];core=[l for l in ls if l['rt_core_GA']];sensitive=[l for l in ls if l['rt_core_sensitivity_only']];ex=[r for r in rs if r['assay'].startswith('EXACT')];ap=[r for r in rs if not r['assay'].startswith('EXACT')]
        ex_out=[r for r in ex if r['high_complexity'] and r['full_span_outside_predicted_CDS']];ap_out=[r for r in ap if r['high_complexity'] and r['full_span_outside_predicted_CDS']]
        if core and ex_out:label='RT_CORE_GA_PLUS_EXACT_REPEAT_ARCHITECTURE_FOR_REVIEW'
        elif core and ap_out:label='RT_CORE_GA_PLUS_APPROXIMATE_REPEAT_SENSITIVITY_FOR_REVIEW'
        elif core:label='RT_CORE_GA_WITHOUT_THIS_PREDICTION_RELATIVE_NONCODING_REPEAT_EVIDENCE'
        elif sensitive and (ex_out or ap_out):label='RT_CORE_SENSITIVITY_PLUS_REPEAT_CONTEXT_FOR_REVIEW'
        elif sensitive:label='RT_CORE_SENSITIVITY_ONLY'
        elif ls:label='RT_LIKE_AUXILIARY_PROFILE_CONTEXT'
        elif ex_out or ap_out:label='REPEAT_CONTEXT_RT_NOT_DETECTED_BY_THIS_PANEL'
        else:label='NO_JOINT_EVIDENCE_IN_THESE_LIMITED_ASSAYS'
        for l in ls:
            for r in rs:pairs.append(pair_context(c,l,r))
        for g in gs:links.append(dict(candidate_key=k,record_key=g['record_key'],sequence_id=g['sequence_id'],start0=g['start0'],end0=g['end0'],predicted_strand=g['strand'],distance_to_math_candidate_bp=distance((g['start0'],g['end0']),(int(c['start0']),int(c['end0'])))))
        added=dict(candidate_key=k,n_predicted_neighborhood_ORFs=len(gs),n_RT_core_GA_loci=len(core),n_RT_core_GA_unique_AA=len({x['aa_sha256'] for x in core}),n_RT_core_sensitivity_only_loci=len(sensitive),n_auxiliary_profile_loci=sum(bool(x['rt_aux_GA'] or x['rt_aux_sensitivity_only']) for x in ls),n_exact_repeat_group_descriptions=len(ex),n_approximate_repeat_group_descriptions=len(ap),n_high_complexity_exact_fullspan_outside_predicted_CDS=len(ex_out),n_high_complexity_approximate_fullspan_outside_predicted_CDS=len(ap_out),n_exact_outside_groups_overlapping_math_candidate=sum(r['overlap_candidate_bp']>0 for r in ex_out),n_approximate_outside_groups_overlapping_math_candidate=sum(r['overlap_candidate_bp']>0 for r in ap_out),recurrence_window_predicted_CDS_overlap_bp=w['window_predicted_CDS_overlap_bp'],recurrence_kmer_pair_matches=w['n_matching_kmer_pairs'],recurrence_paired_blocks=w['n_overlapping_pair_blocks'],left_neighborhood_clipped=ct['left_clipped'],right_neighborhood_clipped=ct['right_clipped'],evidence_label=label,protein_source='UNIFORM_PRODIGAL_NOT_OFFICIAL',known_family_novelty='NOT_ASSESSED',RNA_structure='NOT_COMPUTED',ncRNA_expression='NOT_ASSESSED',partner_function='NOT_PROVEN',ART_identity='NOT_CONFIRMED',negative_call_authorized=False)
        summary.append(dict(c,**added))
        cards.append({'candidate_key':k,'original_mathematical_record':c,'context':ct,'predicted_neighborhood_gene_keys':[g['record_key'] for g in gs],'RT_profile_loci':ls,'exact_repeat_group_ids':[r['repeat_group_id'] for r in ex],'substitution_repeat_group_ids':[r['repeat_group_id'] for r in ap],'recurrence_match_witness':w,'evidence_summary':added,'next_required':['Validate full RT/domain architecture and competing system identities','Verify complete noncoding array architecture including indel/structure alternatives','Compare known ART/RT/partner families and close homologues before novelty claims','Only then prioritize system-level wet-lab candidates']})
    require(len(summary)==len(cands) and all(all(out[key]==c[key] for key in c) for c,out in zip(cands,summary)),'Original candidate fields/order changed')
    write_tsv(ROOT/'outputs/ALL_3674_SYSTEM_EVIDENCE.tsv',summary,list(summary[0]));write_tsv(ROOT/'outputs/RT_REPEAT_SYSTEM_CONTEXT.tsv',pairs,PAIR_FIELDS);write_tsv(ROOT/'outputs/CANDIDATE_PREDICTED_ORF_LINKS.tsv',links,['candidate_key','record_key','sequence_id','start0','end0','predicted_strand','distance_to_math_candidate_bp'])
    write_tsv(ROOT/'outputs/RT_CORE_GA_CANDIDATE_VIEW.tsv',[r for r in summary if r['n_RT_core_GA_loci']>0],list(summary[0]))
    atomic_text(ROOT/'outputs/ALL3674_EVIDENCE_CARDS.jsonl',''.join(canonical_json(c)+'\n' for c in cards))
    neighbors=[]
    for l in loci:
        for g in gacc[l['sequence_id']]:
            if l['record_key']==g['record_key']:continue
            d=distance((l['start0'],l['end0']),(g['start0'],g['end0']))
            if d<=policy()['neighborhood_bp']:neighbors.append(dict(rt_record_key=l['record_key'],neighbor_record_key=g['record_key'],sequence_id=l['sequence_id'],start0=g['start0'],end0=g['end0'],predicted_strand=g['strand'],distance_bp=d,neighbor_is_proven_partner=False,aa_sha256=g['aa_sha256']))
    write_tsv(ROOT/'outputs/RT_NEIGHBOR_ORF_INVENTORY.tsv',neighbors,['rt_record_key','neighbor_record_key','sequence_id','start0','end0','predicted_strand','distance_bp','neighbor_is_proven_partner','aa_sha256'])
    return summary,pairs

def preserve_rt_sources(a05,sources,loci):
    found={l['sequence_id'] for l in loci};rr=[]
    for s in sources:
        if s['accession'] not in found:continue
        acc=s['accession'];dest=ROOT/'outputs/rt_source_evidence'/acc;dest.mkdir(parents=True,exist_ok=True)
        src=safe_path(a05,s['fasta_path']);dst=dest/'genome.fna';shutil.copyfile(src,dst);verify_file(dst,s['fasta_sha256'])
        for name in ['genes.gff','cds.fna','proteins.faa','COMPLETE.json']:shutil.copyfile(ROOT/'work/prodigal'/acc/name,dest/name)
        for f in sorted(dest.iterdir()):rr.append(dict(sequence_id=acc,path=str(f.relative_to(ROOT)),size_bytes=f.stat().st_size,sha256=sha_file(f)))
    write_tsv(ROOT/'outputs/RT_SOURCE_EVIDENCE_MANIFEST.tsv',rr,['sequence_id','path','size_bytes','sha256'])

def finish(sources,cands,genes,targets,loci,reps,summary,refs,refrepeat,lock):
    p=policy();core=[r for r in loci if r['rt_core_GA']];artifact_paths=[f for f in (ROOT/'outputs').rglob('*') if f.is_file() and f.name not in {'GATE_STAGE05C1.json','SUMMARY_STAGE05C1.json','OUTPUT_MANIFEST.tsv'} and not f.name.endswith('.partial')]
    om=output_manifest(ROOT,artifact_paths,ROOT/'outputs/OUTPUT_MANIFEST.tsv')
    gate=dict(stage='Stage05C1',gate='PASS_ALL3674_POSTFREEZE_RT_ORF_REPEAT_EVIDENCE_WITH_SCOPE_LIMITS',gate_pass=True,n_genomes=len(sources),n_frozen_candidates_retained=len(cands),n_full_genome_predicted_ORFs=len(genes),n_candidate_neighborhood_predicted_ORFs=sum(bool(g['candidate_keys']) for g in genes),n_unique_scan_targets_including_6_references=len(targets),n_whole_genome_RT_core_GA_loci=len(core),n_RT_core_GA_loci_captured_by_math_neighborhoods=sum(bool(r['candidate_keys']) for r in core),n_RT_core_GA_loci_outside_math_neighborhoods=sum(not r['candidate_keys'] for r in core),n_candidates_with_RT_core_GA=sum(r['n_RT_core_GA_loci']>0 for r in summary),n_candidates_with_exact_repeat_groups=sum(r['n_exact_repeat_group_descriptions']>0 for r in summary),n_candidates_with_approximate_repeat_groups=sum(r['n_approximate_repeat_group_descriptions']>0 for r in summary),evidence_label_counts=dict(sorted(Counter(r['evidence_label'] for r in summary).items())),n_reference_RT_core_GA=sum(r['rt_core_GA'] for r in refs),n_registered_repeat_cores_exact_recovered=sum(r['exact_registered_core_recovered'] for r in refrepeat),n_registered_repeat_cores_substitution_recovered=sum(r['substitution_registered_core_recovered'] for r in refrepeat),reference_calibration_is_independent_accuracy_estimate=False,mathematical_scores_or_ranks_recomputed=False,candidate_filtering_performed=False,repeat_sensitivity_is_new_exploratory_annotation=True,protein_source='UNIFORM_PRODIGAL_NOT_OFFICIAL',RNA_structure_computed=False,imperfect_repeat_search_exhaustive=False,near_homology_novelty_audit_complete=False,new_ART_systems_confirmed=0,wet_lab_ready_candidates_declared=0,negative_calls_authorized=False,old240_or383_modified=False,UG27_touched=False,sources_unchanged=True,run_lock_sha256=lock,output_manifest_sha256=om,source_score_archive_sha256=p['source_score_archive_sha256'])
    write_json(ROOT/'outputs/SUMMARY_STAGE05C1.json',gate);write_json(ROOT/'outputs/GATE_STAGE05C1.json',gate);print(json.dumps(gate,indent=2),flush=True)

def run(a05,b05,cpus):
    sources,cands,lock=preflight(a05,b05,cpus)
    # These are known reference capability checks, not tuning or new discovery.
    refrepeat=reference_repeat_calibration()
    genes,aas=predict_all(sources,cands,a05,cpus,lock);targets,fan,refs,targetlock=build_targets(genes,aas,lock);domains=hmm_scan(targets,cpus,lock);loci,refresult=rt_loci(genes,refs,domains)
    reps,w,contexts=repeat_all(sources,cands,genes,a05,cpus,lock);summary,pairs=integrate(sources,cands,genes,loci,reps,w,contexts);preserve_rt_sources(a05,sources,loci)
    # Repeat the complete read-only input/score integrity check before closeout.
    input_check(a05,b05);package_check();require(tools()==read_json(ROOT/'outputs/RUN_LOCK.json')['tools'],'Tool binary/version changed during execution');finish(sources,cands,genes,targets,loci,reps,summary,refresult,refrepeat,lock)

def collect():
    package_check();ts=time.strftime('%Y%m%d_%H%M%S');dest=ROOT/f'ART_STAGE05C1_NEW919_BIOLOGY_RESULTS_{ts}.tar.gz';require(not dest.exists(),'Result archive already exists')
    include=['scripts','config','refs','inputs','docs','tests','evidence','outputs','logs','README.md','PACKAGE_MANIFEST.tsv','PACKAGE_MANIFEST.sha256','run_selftest.sh','run_stage05c1.sh','collect_results.sh']
    scope={'all919_raw_genomes_included':False,'raw_full919_mathematical_window_boards_included':False,'all_scan_amino_acid_targets_included':True,'all_predicted_gene_registry_included':True,'all3674_candidate_evidence_included_if_run_completed':True,'RT_supported_source_genomes_and_raw_predictions_included':True,'all_genome_ORF_raw_files_remain_server':'work/prodigal','all_candidate_repeat_checkpoints_remain_server':'work/repeats','no_scoring_or_threshold_change_performed_by_collection':True}
    write_json(ROOT/'evidence/RESULT_COLLECTION_SCOPE.json',scope);paths=[]
    for item in include:
        q=ROOT/item
        if q.is_dir():paths.extend(f for f in q.rglob('*') if f.is_file() and '__pycache__' not in f.parts and not f.name.endswith('.partial'))
        elif q.is_file():paths.append(q)
    manifest=ROOT/'evidence/RESULT_FILES.tsv';paths=[p for p in paths if p!=manifest];output_manifest(ROOT,paths,manifest);paths.append(manifest)
    with tarfile.open(dest.with_suffix(dest.suffix+'.partial'),'w:gz') as tf:
        for q in sorted(set(paths)):require(not q.is_symlink(),'Symlink refused in result package');tf.add(q,arcname=str(q.relative_to(ROOT)),recursive=False)
    os.replace(dest.with_suffix(dest.suffix+'.partial'),dest);h=sha_file(dest);atomic_text(str(dest)+'.sha256',h+'  '+dest.name+'\n');print(dest);print(h)

def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='command',required=True);r=sub.add_parser('run');r.add_argument('--source-root',required=True);r.add_argument('--score-root',required=True);r.add_argument('--cpus',type=int,default=4);sub.add_parser('collect');sub.add_parser('calibrate-repeats');a=ap.parse_args()
    (ROOT/'outputs').mkdir(exist_ok=True);(ROOT/'logs').mkdir(exist_ok=True);(ROOT/'evidence').mkdir(exist_ok=True)
    with (ROOT/'.write.lock').open('a') as f:
        try:fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:print('[HOLD] Another job is using this package root',file=sys.stderr);return 2
        try:
            if a.command=='collect':collect()
            elif a.command=='calibrate-repeats':package_check();print(json.dumps(reference_repeat_calibration(),indent=2))
            else:run(Path(a.source_root),Path(a.score_root),a.cpus)
        except Exception as exc:
            if a.command=='run':write_json(ROOT/'outputs/GATE_STAGE05C1.json',{'stage':'Stage05C1','gate':'HOLD_STAGE05C1_RUNTIME_OR_INTEGRITY_CHECK','gate_pass':False,'error':type(exc).__name__+': '+str(exc),'new_ART_systems_confirmed':0,'negative_calls_authorized':False,'old240_or383_modified':False,'UG27_touched':False})
            import traceback;traceback.print_exc();return 2
    return 0
if __name__=='__main__':raise SystemExit(main())
