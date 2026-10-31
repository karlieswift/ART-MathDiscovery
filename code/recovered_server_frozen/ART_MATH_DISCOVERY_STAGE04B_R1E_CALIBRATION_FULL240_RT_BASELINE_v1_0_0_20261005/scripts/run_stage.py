#!/usr/bin/env python3
"""Calibration and uniform all240 RT-only baseline. No mathematical rescoring."""
from __future__ import annotations
import argparse, fcntl, json, os, re, shutil, subprocess, sys, random, math
from pathlib import Path
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from common import *
from repeat_evidence import find_exact_core_groups, entropy

ROOT=Path(__file__).resolve().parents[1]
AA=set('ACDEFGHIKLMNPQRSTVWY')

def require(condition, message):
    if not condition: raise RuntimeError(message)

def verify_package(root=ROOT):
    assert_manifest(root,root/'PACKAGE_MANIFEST.tsv')

def lock_json(path, obj):
    if path.exists(): require(read_json(path)==obj,'Frozen contract changed: '+str(path))
    else: write_json(path,obj)
    return sha_file(path)

def tools():
    d={}
    for name,args,pat in [('hmmsearch',['-h'],r'HMMER\s+3\.4\b'),('prodigal',['-v'],r'Prodigal\s+V2\.6\.3')]:
        p=shutil.which(name);require(p is not None,'Missing '+name+'; activate art_rt_evidence')
        z=subprocess.run([p]+args,capture_output=True,text=True)
        txt=z.stdout+z.stderr;require(re.search(pat,txt),'Wrong '+name+' version: '+txt[:160])
        d[name]={'path':p,'sha256':sha_file(p),'version_text':'\n'.join(txt.splitlines()[:5])}
    return d

def package_data(root=ROOT):
    pol=read_json(root/'config/METHOD_POLICY.json')
    universe=read_tsv(root/'inputs/STAGE04A_UNIVERSE_MANIFEST.tsv')
    cand=read_tsv(root/'inputs/candidates.tsv')
    require(len(universe)==240 and len({r['sequence_id'] for r in universe})==240,'Bad240manifest')
    require(len(cand)==383 and len({r['sequence_id'] for r in cand})==174,'Bad383registry')
    require(sha_file(root/'inputs/candidates.tsv')==pol['candidate_registry_sha256'],'Candidate hash changed')
    require(sha_file(root/'inputs/STAGE04A_UNIVERSE_MANIFEST.tsv')==pol['stage04a_manifest_sha256'],'Universe hash changed')
    ll={r['sequence_id']:int(r['length_bp']) for r in universe}
    for c in cand:
        require(0<=int(c['start0'])<int(c['end0'])<=ll[c['sequence_id']],'Candidate interval invalid')
        require(int(c['length_bp'])==ll[c['sequence_id']],'Candidate genome length changed')
    return pol,universe,cand,ll

def masks(cand,lengths,radius=10000):
    by=defaultdict(list)
    for c in cand:
        a=c['sequence_id'];by[a].append((max(0,int(c['start0'])-radius),min(lengths[a],int(c['end0'])+radius)))
    return {a:merge_intervals(by[a]) for a in lengths}

def coverage_audit(root=ROOT):
    pol,u,cand,lens=package_data(root);raw=masks(cand,lens,0);flank=masks(cand,lens)
    rr=[]
    for r in u:
        acc=r['sequence_id'];L=lens[acc]
        rr.append({'accession':acc,'length_bp':L,'n_frozen_candidates':sum(c['sequence_id']==acc for c in cand),
                   'candidate_union_bp':sum(y-x for x,y in raw[acc]),'neighborhood_union_bp':sum(y-x for x,y in flank[acc]),
                   'neighborhood_fraction':sum(y-x for x,y in flank[acc])/L})
    out=root/'outputs';out.mkdir(exist_ok=True)
    write_tsv(out/'FULL240_FROZEN_COVERAGE.tsv',rr,list(rr[0]))
    j={'n_genomes':240,'n_candidate_bearing_genomes':174,'n_genomes_without_candidates':66,'n_regions':383,
       'total_bp':sum(lens.values()),'candidate_union_bp':sum(r['candidate_union_bp'] for r in rr),
       'neighborhood_union_bp':sum(r['neighborhood_union_bp'] for r in rr),
       'interpretation':'Region count compression is not DNA footprint, runtime reduction, or biological sensitivity.'}
    j['candidate_union_fraction']=j['candidate_union_bp']/j['total_bp'];j['neighborhood_union_fraction']=j['neighborhood_union_bp']/j['total_bp']
    write_json(out/'COVERAGE_SUMMARY.json',j)
    return j

def repeat_calibration(root=ROOT):
    pol=read_json(root/'config/METHOD_POLICY.json');rr=[];summary=[]
    for r in read_tsv(root/'inputs/controls/REPEAT_REFERENCE_CONTEXT.tsv'):
        seq=next(fasta(root/r['path']))[2];require(sha_text(seq)==r['context_sha256'],'Control DNA changed')
        offset=int(r['context_start0']);core=(int(r['registered_core_start0']),int(r['registered_core_end0']))
        found=find_exact_core_groups(seq,pol['repeat_survey']);n_overlap=0;n_high=0
        for g in found:
            parts=[[offset+s,offset+s+g['repeat_len']] for s in g['starts']];span=(parts[0][0],parts[-1][1])
            ov=overlaps(core,span);hi=entropy(g['core'])>=pol['repeat_survey']['complexity_high_entropy_bits_min']
            copyov=sum(overlaps(core,p) for p in parts)
            # Span alone could bridge across a core: require at least 3 copies actually overlap it.
            recovered=(copyov>=3)
            n_overlap+=recovered;n_high+=recovered and hi
            rr.append({'reference_id':r['reference_id'],'core_nt':g['core'],'n_copies':len(parts),'copy_intervals_json':canonical_json(parts),
                       'span_start0':span[0],'span_end0':span[1],'span_overlaps_registered_core':int(ov),'n_copies_overlapping_registered_core':copyov,
                       'high_complexity':int(hi),'recovers_3copy_registered_core':int(recovered),'RNA_expression_assayed':False})
        summary.append({'reference_id':r['reference_id'],'n_context_repeat_groups':len(found),'n_groups_3copies_overlapping_registered_core':n_overlap,
                        'n_high_complexity_groups_3copies_overlapping_registered_core':n_high,'registered_repeat_core_recovered':int(n_high>0),
                        'role':'known_positive_assay_coverage_not_novel_discovery_or_independent_specificity'})
    write_tsv(root/'outputs/REFERENCE_REPEAT_GROUPS.tsv',rr,list(rr[0]) if rr else ['reference_id'])
    write_tsv(root/'outputs/REFERENCE_REPEAT_CALIBRATION.tsv',summary,list(summary[0]))
    write_json(root/'outputs/REPEAT_CALIBRATION_EXECUTION.json',{'technical_execution_complete':True,'n_reference_contexts':len(summary),
        'n_reference_cores_recovered':sum(r['registered_repeat_core_recovered'] for r in summary),
        'thresholds_changed':False,'known_truth_used_for_assay_evaluation_only':True,'mathematical_scores_recomputed':False})
    return summary

def preflight_genomes(a04,root=ROOT):
    pol,u,c,lens=package_data(root)
    require(sha_file(a04/'outputs/STAGE04A_UNIVERSE_MANIFEST.tsv')==pol['stage04a_manifest_sha256'],'Actual Stage04A universe not canonical')
    sources=[]
    for r in u:
        rel=Path(r['frozen_fasta']);require(not rel.is_absolute() and '..' not in rel.parts,'Unsafe manifest source path')
        src=(a04/rel).resolve();require(src.is_relative_to(a04),'Source outside Stage04A root')
        require(src.is_file(),'Frozen FASTA missing (no guessing or new download): '+str(src))
        require(sha_file(src)==r['fasta_sha256'],'Frozen FASTA file hash mismatch '+r['sequence_id'])
        rec=list(fasta(src));require(len(rec)==1,'Single-record FASTA required')
        seq=rec[0][2];require(len(seq)==lens[r['sequence_id']] and not(set(seq)-set('ACGT')),'Genome length/alphabet mismatch')
        dest=root/'inputs/genomes'/f"{r['sequence_id']}.fna"
        # Exact nucleotide sequence; deterministic header for Prodigal, with both byte identities recorded.
        if dest.exists():require(list(fasta(dest))==[(r['sequence_id'],r['sequence_id'],seq)],'Previously staged genome changed')
        else:write_fasta(dest,[(r['sequence_id'],seq)])
        sources.append({'accession':r['sequence_id'],'source_fasta':str(src),'source_file_sha256':r['fasta_sha256'],
                        'sequence_sha256':sha_text(seq),'staged_file_sha256':sha_file(dest),'length_bp':len(seq)})
    write_tsv(root/'outputs/FULL240_SOURCE_PROVENANCE.tsv',sources,list(sources[0]));return sources

def parse_prediction(root,acc):
    d=root/'work/prodigal'/acc;seq=next(fasta(root/'inputs/genomes'/f'{acc}.fna'))[2]
    aa=fasta_dict(d/'proteins.faa');nt=fasta_dict(d/'cds.fna')
    require(set(aa)==set(nt) and len(aa)>0,'Predicted sequence inventory mismatch '+acc)
    txt=(d/'genes.gff').read_text();tt=set(re.findall(r'transl_table=(\d+)',txt));require(len(tt)==1,'Missing translation table')
    loci=set()
    for l in txt.splitlines():
        if not l or l.startswith('#'):continue
        z=l.split('\t');require(len(z)==9 and z[0]==acc and z[2]=='CDS','Invalid GFF')
        loci.add((int(z[3])-1,int(z[4]),z[6]))
    rows=[];seqs={}
    def header(desc):
        x=[z.strip() for z in desc.split('#')];require(len(x)>=5,'Invalid Prodigal header')
        attrs=dict(z.split('=',1) for z in x[4].split(';') if '=' in z)
        require(int(x[3]) in (-1,1),'Bad Prodigal strand')
        return int(x[1])-1,int(x[2]),'+' if int(x[3])==1 else '-',attrs
    for i,(gid,(desc,protein)) in enumerate(aa.items(),1):
        a,b,strand,attrs=header(desc);a2,b2,s2,_=header(nt[gid][0])
        require((a,b,strand)==(a2,b2,s2) and (a,b,strand) in loci,'CDS/protein/GFF coordinate mismatch')
        require(0<=a<b<=len(seq),'Predicted locus outside genome')
        dna=seq[a:b] if strand=='+' else rc(seq[a:b]);require(dna==nt[gid][1],'Predicted nucleotide mismatch')
        if protein.endswith('*'):protein=protein[:-1]
        require(protein and not(set(protein)-AA),'Predicted protein has unresolved symbols; no replacement permitted')
        rk=acc+f'__PRODIGAL_{i:05d}';seqs[rk]=protein
        rows.append({'record_key':rk,'sequence_id':acc,'prodigal_id':gid,'start0':a,'end0':b,'strand':strand,
                     'translation_table':next(iter(tt)),'partial':attrs.get('partial',''),'aa_length':len(protein),
                     'aa_sha256':sha_text(protein),'nt_sha256':sha_text(dna)})
    require(len(loci)==len(rows),'Prediction count mismatch')
    return rows,seqs

def predictions(root,prior,universe,cpus,runhash):
    locks=read_tsv(root/'config/PRIOR_PRODIGAL_174_RAW_LOCK.tsv'); byacc=defaultdict(list)
    for r in locks:byacc[Path(r['path']).parts[2]].append(r)
    require(len(byacc)==174,'Prior raw lock incomplete')
    # Check every reused file before starting any new prediction.
    for r in locks:require((prior/r['path']).is_file() and sha_file(prior/r['path'])==r['sha256'],'R1D2 raw prediction changed: '+r['path'])
    def one(acc):
        d=root/'work/prodigal'/acc;d.mkdir(parents=True,exist_ok=True)
        fp={'run_lock':runhash,'genome_sha256':sha_file(root/'inputs/genomes'/f'{acc}.fna'),'mode':'meta','reused_R1D2':acc in byacc}
        done=d/'R1E_COMPLETE.json';files=[d/n for n in ('proteins.faa','cds.fna','genes.gff')]
        if done.exists():
            old=read_json(done);require(old['fingerprint']==fp,'Prediction run contract changed')
            require(all(p.exists() and sha_file(p)==old['hashes'].get(p.name) for p in files),'Completed prediction corrupted; preserve and inspect')
        else:
            if acc in byacc:
                for r in byacc[acc]:shutil.copyfile(prior/r['path'],d/Path(r['path']).name)
            else:
                cmd=['prodigal','-i',root/'inputs/genomes'/f'{acc}.fna','-p','meta','-f','gff','-o',d/'genes.gff','-a',d/'proteins.faa','-d',d/'cds.fna','-q']
                run_cmd(cmd,root/'logs'/f'prodigal_{acc}.log')
        rr,ss=parse_prediction(root,acc)
        if not done.exists():write_json(done,{'fingerprint':fp,'hashes':{p.name:sha_file(p) for p in files},'n_genes':len(rr)})
        return acc,rr,ss
    per={}
    with ThreadPoolExecutor(max_workers=cpus) as ex:
        futures=[ex.submit(one,r['sequence_id']) for r in universe]
        for i,f in enumerate(as_completed(futures),1):
            acc,rr,ss=f.result();per[acc]=(rr,ss);print(f'[ORF {i}/240] {acc} n={len(rr)} reused={acc in byacc}',flush=True)
    rows=[];seqs={}
    for acc in sorted(per):rr,ss=per[acc];rows+=rr;seqs.update(ss)
    return rows,seqs

def scan_targets(root,rows,seqs,cand,lens):
    by=defaultdict(list)
    for c in cand:by[c['sequence_id']].append(c)
    fan=[];target={}
    for r in rows:
        ck=[c['sequence_id']+'::'+c['region_id'] for c in by[r['sequence_id']] if overlaps((r['start0'],r['end0']),(max(0,int(c['start0'])-10000),min(lens[r['sequence_id']],int(c['end0'])+10000)))]
        r['candidate_keys']=';'.join(ck);r['in_frozen_neighborhood']=int(bool(ck));r['target_id']='AA_'+r['aa_sha256']
        target[r['target_id']]=seqs[r['record_key']]
        fan.append({'record_key':r['record_key'],'target_id':r['target_id'],'branch':'PREDICTED_FULL240','sequence_id':r['sequence_id'],'candidate_keys':r['candidate_keys']})
    offs=fasta_dict(root/'inputs/official_proteins.faa')
    for r in read_tsv(root/'inputs/official_protein_registry.tsv'):
        k=r['protein_key'];s=offs[k][1];require(sha_text(s)==r['aa_sha256'],'Official amino acid source mismatch')
        t='AA_'+r['aa_sha256'];target[t]=s;fan.append({'record_key':k,'target_id':t,'branch':'OFFICIAL_NEIGHBORHOOD_AUX','sequence_id':r['sequence_id'],'candidate_keys':r['candidate_keys']})
    refs=fasta_dict(root/'inputs/controls/RT_REFERENCE_PANEL.faa')
    for r in read_tsv(root/'inputs/controls/RT_REFERENCE_PANEL.tsv'):
        s=refs[r['reference_id']][1];require(sha_text(s)==r['aa_sha256'],'Reference protein mismatch');t='AA_'+sha_text(s);target[t]=s
        fan.append({'record_key':r['reference_id'],'target_id':t,'branch':'REFERENCE_ONLY','sequence_id':r['source_accession'],'candidate_keys':''})
    write_tsv(root/'outputs/FULL240_PREDICTED_GENE_REGISTRY.tsv',rows,list(rows[0]))
    staged=root/'work/target_staging';staged.mkdir(parents=True,exist_ok=True)
    tf=staged/'SCAN_TARGETS.faa';ff=staged/'TARGET_FANOUT.tsv'
    write_tsv(ff,fan,list(fan[0]));write_fasta(tf,sorted(target.items()))
    j={'target_sha256':sha_file(tf),'fanout_sha256':sha_file(ff),'n_unique_targets':len(target)}
    lp=root/'outputs/TARGET_LOCK.json'
    if lp.exists():
        require(read_json(lp)==j,'Target universe changed after freeze')
        require(sha_file(root/'outputs/SCAN_TARGETS.faa')==j['target_sha256'] and sha_file(root/'outputs/TARGET_FANOUT.tsv')==j['fanout_sha256'],'Frozen scan targets corrupted; preserve and inspect')
    else:
        tf.replace(root/'outputs/SCAN_TARGETS.faa');ff.replace(root/'outputs/TARGET_FANOUT.tsv');write_json(lp,j)
    return target,fan

DOMAIN_FIELDS=['target_id','profile','role','source_scan','tier','full_evalue','independent_domain_evalue','full_bits','domain_bits','hmm_coverage','ali_from1','ali_to1']
def parse_domains(path,catalog,targets,mode):
    mm={r['name']:r for r in catalog};rr=[]
    for l in path.read_text().splitlines():
        if not l or l.startswith('#'):continue
        v=l.split(maxsplit=22);require(len(v)>=22,'Truncated HMMER row')
        require(v[0] in targets and v[3] in mm,'Unknown HMMER target or model')
        m=mm[v[3]];L=int(v[2]);M=int(v[5]);a,b=int(v[17]),int(v[18]);h0,h1=int(v[15]),int(v[16])
        require(L==len(targets[v[0]]) and M==int(m['length']) and v[4]==m['accession_version'],'HMMER identity mismatch')
        require(1<=h0<=h1<=M and 1<=a<=b<=L,'Invalid HMM coordinates')
        fullE,indE,fullB,domB=float(v[6]),float(v[12]),float(v[7]),float(v[13]);cov=(h1-h0+1)/M
        require(all(math.isfinite(x) for x in (fullE,indE,fullB,domB,cov)),'Nonfinite HMMER output')
        tier='REPORTED_ONLY'
        if mode=='GA':
            require(fullB+0.051>=float(m['ga_sequence_bits']) and domB+0.051>=float(m['ga_domain_bits']),'GA row violates model thresholds')
            tier='CORE_GA' if m['role']=='RT_CORE' else 'AUX_GA'
        elif fullE<=.001 and indE<=.001 and cov>=.35:tier='CORE_SENSITIVITY' if m['role']=='RT_CORE' else 'AUX_SENSITIVITY'
        rr.append(dict(zip(DOMAIN_FIELDS,[v[0],v[4],m['role'],mode,tier,fullE,indE,fullB,domB,cov,a,b])))
    return rr

def hmmer(root,targets,cpus,runhash):
    model=root/'refs/RT_PANEL.hmm';catalog=read_tsv(root/'refs/MODEL_CATALOG.tsv');out=root/'outputs/hmmer';out.mkdir(exist_ok=True)
    rows=[]
    for mode in ('GA','SENSITIVITY'):
        report=out/(mode+'.report.txt');dom=out/(mode+'.domtbl');tbl=out/(mode+'.tbl');done=out/(mode+'.COMPLETE.json')
        fp={'runhash':runhash,'targets_sha256':sha_file(root/'outputs/SCAN_TARGETS.faa'),'models_sha256':sha_file(model),'Z':100000,'mode':mode}
        if done.exists():
            j=read_json(done);require(j['fingerprint']==fp,'HMMER targets or settings changed')
            require(all(p.exists() and sha_file(p)==j['hashes'].get(p.name) for p in (report,dom,tbl)),'HMMER completed result changed')
        else:
            cmd=['hmmsearch','--cpu',str(cpus),'--seed','42','-Z','100000','--noali','-o',report,'--domtblout',dom,'--tblout',tbl]
            cmd+=['--cut_ga'] if mode=='GA' else ['-E','10','--domE','10']
            cmd +=[model,root/'outputs/SCAN_TARGETS.faa'];print('[HMMER]',mode,'N=',len(targets),flush=True);run_cmd(cmd,root/'logs'/f'hmmer_{mode}.log')
        text=report.read_text();queries=re.findall(r'^Query:\s+(\S+)\s+\[M=\d+\]',text,re.M)
        require(sorted(queries)==sorted(r['name'] for r in catalog) and text.rstrip().endswith('[ok]'),'HMMER did not finish all models')
        counts=re.findall(r'Target sequences:\s+(\d+)',text);require(len(counts)==3 and all(int(z)==len(targets) for z in counts),'HMMER target count incomplete')
        rr=parse_domains(dom,catalog,targets,mode);rows+=rr
        if not done.exists():write_json(done,{'fingerprint':fp,'hashes':{p.name:sha_file(p) for p in (report,dom,tbl)}})
    write_tsv(root/'outputs/ALL_REPORTED_DOMAINS.tsv',rows,DOMAIN_FIELDS);return rows

def position_null(cand,lens,rt_loci,replicates=1000,seed=20261005):
    by=defaultdict(list)
    for c in cand:by[c['sequence_id']].append(int(c['end0'])-int(c['start0']))
    rt=defaultdict(list)
    for r in rt_loci:rt[r['sequence_id']].append((int(r['start0']),int(r['end0'])))
    observed=masks(cand,lens)
    obs=sum(any(overlaps(locus,interval) for interval in observed[acc]) for acc,ls in rt.items() for locus in ls)
    rng=random.Random(seed);vals=[]
    for b in range(replicates):
        shifted={}
        for acc,widths in by.items():
            intervals=[]
            for width in widths:
                start=rng.randrange(lens[acc]-width+1);intervals.append((max(0,start-10000),min(lens[acc],start+width+10000)))
            shifted[acc]=merge_intervals(intervals)
        k=sum(any(overlaps(locus,interval) for interval in shifted.get(acc,[])) for acc,ls in rt.items() for locus in ls)
        vals.append({'replicate':b+1,'captured_predicted_core_loci':k})
    arr=[r['captured_predicted_core_loci'] for r in vals];mean=sum(arr)/len(arr)
    summary={'observed_captured_predicted_core_loci':obs,'random_mean':mean,'random_min':min(arr),'random_max':max(arr),
             'enrichment_vs_random_mean':obs/mean if mean else None,'monte_carlo_upper_tail_p':(1+sum(x>=obs for x in arr))/(len(arr)+1) if rt_loci else None,
             'replicates':replicates,'seed':seed,'no_informative_RT_reference_loci':len(rt_loci)==0,
             'interpretation':'exploratory_same-genome_same-width_position_null;not_GC/CDS_matched;not_ART_sensitivity;not_prospective_confirmation'}
    return vals,summary

def summarize(root,rows,fan,domains,cand,lens,repeat):
    tiers=defaultdict(set)
    for r in domains:tiers[r['target_id']].add(r['tier'])
    for r in rows:
        ts=tiers[r['target_id']];r.update({'RT_core_GA':int('CORE_GA' in ts),'RT_core_sensitivity':int('CORE_SENSITIVITY' in ts),
                                         'RT_aux_GA':int('AUX_GA' in ts),'RT_aux_sensitivity':int('AUX_SENSITIVITY' in ts)})
    write_tsv(root/'outputs/FULL240_PREDICTED_RT_ASSESSMENT.tsv',rows,list(rows[0]))
    refmap={r['record_key']:r['target_id'] for r in fan if r['branch']=='REFERENCE_ONLY'}
    refrows=[]
    for r in read_tsv(root/'inputs/controls/RT_REFERENCE_PANEL.tsv'):
        ts=tiers[refmap[r['reference_id']]];out=dict(r);out.update({'RT_core_GA':int('CORE_GA' in ts),
          'RT_core_sensitivity':int('CORE_SENSITIVITY' in ts),'aux_GA':int('AUX_GA' in ts),'aux_sensitivity':int('AUX_SENSITIVITY' in ts),
          'assay_scope':'positive_RT_reference;nonarray_reference_is_not_RT_negative','no_hit_is_model_limit_not_negative_label':True});refrows.append(out)
    write_tsv(root/'outputs/REFERENCE_RT_CAPABILITY.tsv',refrows,list(refrows[0]))
    by=defaultdict(list);off=defaultdict(set)
    for r in rows:
        for ck in filter(None,r['candidate_keys'].split(';')):by[ck].append(r)
    for r in fan:
        if r['branch']=='OFFICIAL_NEIGHBORHOOD_AUX' and 'CORE_GA' in tiers[r['target_id']]:
            for ck in filter(None,r['candidate_keys'].split(';')):off[ck].add(r['target_id'])
    out=[]
    for c in cand:
        ck=c['sequence_id']+'::'+c['region_id'];genes=by[ck]
        z=dict(c);z.update({'candidate_key':ck,'n_uniform_predicted_neighborhood_proteins':len(genes),
          'n_uniform_predicted_RT_core_GA':sum(r['RT_core_GA'] for r in genes),'n_official_auxiliary_branch_RT_core_GA_unique_AA':len(off[ck]),
          'ART_identity':'NOT_ASSESSED','novelty':'NOT_ASSESSED','candidate_dropped':False});out.append(z)
    require(all(all(z[k]==c[k] for k in c) for z,c in zip(out,cand)),'Candidate values changed')
    write_tsv(root/'outputs/ALL383_BASELINE_LINKAGE.tsv',out,list(out[0]))
    core=[r for r in rows if r['RT_core_GA']];inside=[r for r in core if r['in_frozen_neighborhood']]
    write_tsv(root/'outputs/CORE_RT_LOCI_INSIDE_AND_OUTSIDE.tsv',core,list(rows[0]))
    vals,pos=position_null(cand,lens,core);write_tsv(root/'outputs/POSITION_NULL_REPLICATES.tsv',vals,list(vals[0]));write_json(root/'outputs/POSITION_NULL_SUMMARY.json',pos)
    totalAA={r['aa_sha256'] for r in rows};selected=[r for r in rows if r['in_frozen_neighborhood']];selAA={r['aa_sha256'] for r in selected}
    coverage=read_json(root/'outputs/COVERAGE_SUMMARY.json')
    summary={'stage':'Stage04B-R1E','gate':'PASS_REFERENCE_CAPABILITY_AND_FULL240_RT_BASELINE_COMPUTED_WITH_SCOPE_LIMITS','gate_pass':True,
      'n_genomes':240,'n_candidates_retained':383,'n_full240_predicted_protein_records':len(rows),'n_whole_genome_predicted_unique_AA':len(totalAA),
      'n_neighborhood_predicted_records':len(selected),'n_neighborhood_predicted_unique_AA':len(selAA),
      'record_workload_fraction':len(selected)/len(rows),'unique_sequence_workload_fraction':len(selAA)/len(totalAA),
      'n_predicted_core_GA_loci_full240':len(core),'n_core_GA_loci_captured_by_frozen_neighborhoods':len(inside),'n_core_GA_loci_outside':len(core)-len(inside),
      'profile_positive_locus_capture_fraction':len(inside)/len(core) if core else None,
      'n_core_GA_loci_in_zero_candidate_genomes':sum(not any(c['sequence_id']==r['sequence_id'] for c in cand) for r in core),
      'n_RT_reference_controls':len(refrows),'n_ART_associated_reference_core_GA':sum(r['RT_core_GA'] for r in refrows if r['reference_class']=='REGISTERED_ART_ASSOCIATED_RT'),
      'n_nonarray_RT_reference_core_GA':sum(r['RT_core_GA'] for r in refrows if r['reference_class']=='REGISTERED_NONARRAY_RT'),
      'n_registered_ART_repeat_cores_recovered':sum(r['registered_repeat_core_recovered'] for r in repeat),
      'reference_controls_are_independent_specificity_test':False,'ART_sensitivity_estimated':False,'method_superiority_established':False,
      'mathematical_rescoring':False,'retuning':False,'UG27_touched':False,'new_genome_discovery_expansion_performed':False,
      'baseline_uses_existing_240_genome_universe':True,'HMMER_Z':100000,'R1D2_primary_result_unchanged':True,'coverage':coverage}
    write_json(root/'outputs/SUMMARY_STAGE04B_R1E.json',summary);return summary

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage04a');ap.add_argument('--r1d2');ap.add_argument('--cpus',type=int,default=4);ap.add_argument('--precompute-only',action='store_true');a=ap.parse_args()
    require(1<=a.cpus<=32,'Use 1..32 CPUs');verify_package();(ROOT/'outputs').mkdir(exist_ok=True);(ROOT/'logs').mkdir(exist_ok=True)
    with open(ROOT/'work.lock','w') as lf:
        try:fcntl.flock(lf,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError('Another process owns this new-stage directory')
        try:
            coverage_audit();repeat=repeat_calibration()
            if a.precompute_only:
                print(json.dumps({'coverage':read_json(ROOT/'outputs/COVERAGE_SUMMARY.json'),'repeat_calibration':repeat},indent=2));return
            require(a.stage04a and a.r1d2,'Explicit --stage04a and --r1d2 paths required')
            a04=Path(a.stage04a).resolve();prior=Path(a.r1d2).resolve()
            require(a04!=ROOT and prior!=ROOT and not ROOT.is_relative_to(a04) and not ROOT.is_relative_to(prior),'New stage must be a separate directory')
            versions=tools()
            original=read_json(ROOT/'config/ORIGINAL_TOOL_IDENTITIES.json')
            for name in ('prodigal','hmmsearch'):
                require(versions[name]['sha256']==original[name]['sha256'],'Tool binary differs from frozen R1D2: '+name+'; do not silently mix implementations')
            pol,u,c,lens=package_data();src=preflight_genomes(a04)
            for r in read_tsv(ROOT/'config/PRIOR_PRODIGAL_174_RAW_LOCK.tsv'):
                require(sha_file(prior/r['path'])==r['sha256'],'Prior prediction changed')
            runhash=lock_json(ROOT/'outputs/RUN_LOCK.json',{'package_manifest_sha256':sha_file(ROOT/'PACKAGE_MANIFEST.tsv'),
              'policy_sha256':sha_file(ROOT/'config/METHOD_POLICY.json'),'stage04a':str(a04),'prior_R1D2':str(prior),
              'source_rows':src,'tool_versions':versions,'cpus':a.cpus,'new_HMMER_Z':100000,'mathematical_scores_recomputed':False})
            rows,seqs=predictions(ROOT,prior,u,a.cpus,runhash);targets,fan=scan_targets(ROOT,rows,seqs,c,lens)
            domains=hmmer(ROOT,targets,a.cpus,runhash);j=summarize(ROOT,rows,fan,domains,c,lens,repeat)
            for r in src:require(sha_file(r['source_fasta'])==r['source_file_sha256'],'Frozen source changed during run')
            for r in read_tsv(ROOT/'config/PRIOR_PRODIGAL_174_RAW_LOCK.tsv'):require(sha_file(prior/r['path'])==r['sha256'],'Prior source changed during run')
            verify_package();paths=list((ROOT/'outputs').rglob('*'));dest=ROOT/'OUTPUT_MANIFEST.tsv';mh=output_manifest(ROOT,paths,dest)
            j.update({'output_manifest_sha256':mh,'run_lock_sha256':runhash,'sources_unchanged':True})
            write_json(ROOT/'GATE_STAGE04B_R1E.json',j);print(json.dumps(j,indent=2))
        except Exception as e:
            write_json(ROOT/'GATE_STAGE04B_R1E.json',{'gate':'HOLD_STAGE04B_R1E_EXECUTION_OR_PROVENANCE','gate_pass':False,'error':str(e),
                'mathematical_scores_recomputed':False,'UG27_touched':False,'negative_calls_authorized':False});raise
if __name__=='__main__':main()
