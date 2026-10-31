"""Real Prodigal and HMMER execution; strict, explicit source linkage."""
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor,as_completed
from collections import defaultdict
import re,shutil,json,subprocess
from common import *

def tool_versions():
    result={}
    for name,args,pattern in [('prodigal',['-v'],r'V?2\.6\.3'),('hmmsearch',['-h'],r'HMMER\s+3\.4\b'),('hmmemit',['-h'],r'HMMER\s+3\.4\b')]:
        p=shutil.which(name)
        if p is None:raise RuntimeError('Missing required tool '+name+'; use the isolated environment install command')
        cp=subprocess.run([p]+args,capture_output=True,text=True);text=cp.stdout+cp.stderr
        if not re.search(pattern,text):raise RuntimeError('Unexpected version of '+name+': '+text[:180])
        result[name]={'path':p,'sha256':sha_file(p),'version_output':text[:1500] if name=='prodigal' else '\n'.join(text.splitlines()[:4])}
    return result

def parse_prodigal_header(desc):
    a=[x.strip() for x in desc.split('#')]
    if len(a)<5:raise ValueError('Unsupported Prodigal header: '+desc)
    st,en,sign=int(a[1])-1,int(a[2]),int(a[3])
    if sign not in (-1,1) or st<0 or en<=st:raise ValueError('Invalid predicted coordinates')
    attrs={}
    for z in a[4].split(';'):
        if '=' in z:k,v=z.split('=',1);attrs[k]=v
    return a[0].split()[0],st,en,'+' if sign==1 else '-',attrs

ORF_FIELDS=['record_key','sequence_id','prodigal_id','start0','end0','strand','translation_table','partial','start_type','nt_sha256','aa_sha256','aa_length','candidate_keys','source_branch','official_exact_aa_match','official_same_locus_aa_match']

def load_predicted(root,acc):
    d=Path(root)/'work/prodigal'/acc;seq=next(fasta(Path(root)/'inputs/genomes'/f'{acc}.fna'))[2]
    aa=fasta_dict(d/'proteins.faa');nt=fasta_dict(d/'cds.fna');gff=(d/'genes.gff').read_text();tt=set(re.findall(r'transl_table=(\d+)',gff))
    if len(tt)!=1:raise ValueError('Missing/ambiguous Prodigal translation table '+acc)
    tt=int(next(iter(tt)))
    if set(aa)!=set(nt):raise ValueError('Predicted CDS/protein IDs differ '+acc)
    gff_loci=set()
    for l in gff.splitlines():
        if not l or l.startswith('#'):continue
        f=l.split('\t')
        if len(f)!=9 or f[2]!='CDS':raise ValueError('Unexpected Prodigal GFF feature')
        if f[0]!=acc:raise ValueError('GFF accession mismatch')
        gff_loci.add((int(f[3])-1,int(f[4]),f[6]))
    rows=[];aas={}
    for i,(gid,(desc,protein)) in enumerate(aa.items(),1):
        name,st,en,strand,attrs=parse_prodigal_header(desc)
        ngid,nst,nen,nstrand,_=parse_prodigal_header(nt[gid][0])
        if (st,en,strand)!=(nst,nen,nstrand) or (st,en,strand) not in gff_loci:raise ValueError('Predicted GFF/CDS/protein coordinates differ')
        if en>len(seq):raise ValueError('Predicted coordinates outside genome')
        dna=seq[st:en] if strand=='+' else rc(seq[st:en])
        if nt[gid][1]!=dna:raise ValueError('Predicted CDS fails exact coding-strand genomic validation '+acc+'/'+gid)
        if protein.endswith('*'):protein=protein[:-1]
        if not protein or set(protein)-set('ACDEFGHIKLMNPQRSTVWY'):raise ValueError('Noncanonical amino acid in Prodigal output; no replacement allowed')
        rk=acc+f'__PRODIGAL_{i:05d}'
        rows.append({'record_key':rk,'sequence_id':acc,'prodigal_id':gid,'start0':st,'end0':en,'strand':strand,'translation_table':tt,'partial':attrs.get('partial',''),'start_type':attrs.get('start_type',''),'nt_sha256':sha_text(dna),'aa_sha256':sha_text(protein),'aa_length':len(protein)})
        aas[rk]=protein
    if len(gff_loci)!=len(rows):raise ValueError('Prediction count mismatch '+acc)
    return rows,aas

def run_prodigal_all(root,cpus,lock_hash):
    root=Path(root);sources=read_tsv(root/'inputs/genome_sources.tsv');cands=read_tsv(root/'inputs/candidates.tsv');byacc=defaultdict(list)
    for c in cands:byacc[c['sequence_id']].append(c)
    official=read_tsv(root/'inputs/official_protein_registry.tsv');offhash=defaultdict(set);offloc=defaultdict(set)
    for r in official:
        offhash[r['sequence_id']].add(r['aa_sha256']);offloc[(r['sequence_id'],int(r['cds_start0']),int(r['cds_end0']))].add(r['aa_sha256'])
    def one(s):
        acc=s['accession'];d=root/'work/prodigal'/acc;d.mkdir(parents=True,exist_ok=True);done=d/'COMPLETE.json'
        fingerprint={'run_lock':lock_hash,'genome_file_sha256':sha_file(root/'inputs/genomes'/f'{acc}.fna'),'mode':'meta'}
        raw=[d/'proteins.faa',d/'cds.fna',d/'genes.gff']
        valid=False
        if done.exists():
            j=read_json(done);valid=(j.get('fingerprint')==fingerprint and all(p.is_file() and sha_file(p)==j.get('outputs',{}).get(p.name) for p in raw))
        if not valid:
            cmd=['prodigal','-i',root/'inputs/genomes'/f'{acc}.fna','-p','meta','-f','gff','-o',d/'genes.gff','-a',d/'proteins.faa','-d',d/'cds.fna','-q']
            run_cmd(cmd,root/'logs'/f'prodigal_{acc}.log')
        rows,aas=load_predicted(root,acc)
        for r in rows:
            ck=[]
            for c in byacc[acc]:
                interval=(max(0,int(c['start0'])-10000),min(int(c['length_bp']),int(c['end0'])+10000))
                if overlaps((r['start0'],r['end0']),interval):ck.append(acc+'::'+c['region_id'])
            r.update({'candidate_keys':';'.join(ck),'source_branch':'PRODIGAL_RESCUE' if s['annotation_status']=='NO_OFFICIAL_CDS' else 'PRODIGAL_SENSITIVITY','official_exact_aa_match':int(r['aa_sha256'] in offhash[acc]),'official_same_locus_aa_match':int(r['aa_sha256'] in offloc[(acc,r['start0'],r['end0'])])})
        write_json(done,{'fingerprint':fingerprint,'outputs':{p.name:sha_file(p) for p in raw},'n_predicted_genes':len(rows)})
        return acc,rows,aas
    collected={}
    with ThreadPoolExecutor(max_workers=max(1,cpus)) as ex:
        futs=[ex.submit(one,s) for s in sources]
        for i,f in enumerate(as_completed(futs),1):
            acc,rr,ss=f.result();collected[acc]=(rr,ss);print(f'[ORF {i}/{len(sources)}] {acc} genes={len(rr)}',flush=True)
    allrows=[];seqs={}
    for acc in sorted(collected):r,s=collected[acc];allrows.extend(r);seqs.update(s)
    write_tsv(root/'outputs/PRODIGAL_ALL174_GENE_REGISTRY.tsv',allrows,ORF_FIELDS)
    nearby=[r for r in allrows if r['candidate_keys']]
    write_tsv(root/'outputs/PRODIGAL_NEIGHBORHOOD_REGISTRY.tsv',nearby,ORF_FIELDS)
    write_fasta(root/'outputs/PRODIGAL_NEIGHBORHOOD_PROTEINS.faa',[(r['record_key'],seqs[r['record_key']]) for r in nearby])
    return allrows,nearby,seqs

FANOUT_FIELDS=['target_id','record_key','source_branch','sequence_id','official_protein_accession','feature_type','parent_json_path','start0','end0','parts_json','coding_strand','candidate_keys','product','aa_sha256','aa_length']

def build_targets(root,predicted,aa_map):
    root=Path(root);offaa=fasta_dict(root/'inputs/official_proteins.faa');off=read_tsv(root/'inputs/official_protein_registry.tsv');targets={};fan=[]
    for r in off:
        seq=offaa[r['protein_key']][1]
        if sha_text(seq)!=r['aa_sha256']:raise ValueError('Official peptide SHA mismatch')
        tid='AA_'+sha_text(seq);targets[tid]=seq
        fan.append({'target_id':tid,'record_key':r['protein_key'],'source_branch':'OFFICIAL','sequence_id':r['sequence_id'],'official_protein_accession':r['official_protein_accession'],'feature_type':r['feature_type'],'parent_json_path':r['parent_json_path'],'start0':int(r['cds_start0']),'end0':int(r['cds_end0']),'parts_json':canonical_json(parts_from_registry(r)),'coding_strand':r['cds_strand'],'candidate_keys':r['candidate_keys'],'product':r['product'],'aa_sha256':r['aa_sha256'],'aa_length':len(seq)})
    for r in predicted:
        seq=aa_map[r['record_key']];tid='AA_'+sha_text(seq);targets[tid]=seq
        fan.append({'target_id':tid,'record_key':r['record_key'],'source_branch':r['source_branch'],'sequence_id':r['sequence_id'],'official_protein_accession':'','feature_type':'PREDICTED_CDS','parent_json_path':'','start0':r['start0'],'end0':r['end0'],'parts_json':canonical_json([[r['start0'],r['end0']]]),'coding_strand':r['strand'],'candidate_keys':r['candidate_keys'],'product':'Prodigal predicted CDS; function unassigned','aa_sha256':r['aa_sha256'],'aa_length':len(seq)})
    staging=root/'work/target_staging';staging.mkdir(parents=True,exist_ok=True)
    sf=staging/'RT_SCAN_TARGETS.faa';st=staging/'RT_TARGET_FANOUT.tsv'
    write_fasta(sf,sorted(targets.items()));write_tsv(st,fan,FANOUT_FIELDS)
    target_lock={'target_fasta_sha256':sha_file(sf),'fanout_sha256':sha_file(st),'n_unique_targets':len(targets),'n_official_records':len(off),'n_predicted_neighborhood_records':len(predicted),'deduplication_is_computational_only':True,'created_before_first_candidate_scan':True}
    lp=root/'outputs/TARGET_UNIVERSE_LOCK.json'
    if lp.exists():
        if read_json(lp)!=target_lock:raise RuntimeError('Target universe changed after its freeze; do not overwrite or rescore')
        for final,expected in [(root/'outputs/RT_SCAN_TARGETS.faa',target_lock['target_fasta_sha256']),(root/'outputs/RT_TARGET_FANOUT.tsv',target_lock['fanout_sha256'])]:
            if not final.is_file() or sha_file(final)!=expected:raise RuntimeError('Previously frozen scan targets changed or missing: '+str(final))
        sf.unlink();st.unlink()
    else:
        sf.replace(root/'outputs/RT_SCAN_TARGETS.faa');st.replace(root/'outputs/RT_TARGET_FANOUT.tsv');write_json(lp,target_lock)
    return targets,fan

DOMAIN_FIELDS=['target_id','target_length','profile_name','profile_accession_version','profile_role','hmm_length','full_evalue','full_bitscore','domain_index','domain_count','conditional_domain_evalue','independent_domain_evalue','domain_bitscore','hmm_from1','hmm_to1','ali_from1','ali_to1','env_from1','env_to1','hmm_coverage','mean_alignment_posterior','source_scan','evidence_tier']

def parse_domtbl(path,catalog,targets,scan):
    rows=[];names={r['name']:r for r in catalog};policy=read_json(Path(__file__).resolve().parents[1]/'config/METHOD_POLICY.json')['rt_scan']['sensitivity']
    with Path(path).open() as f:
        for l in f:
            if not l.strip() or l.startswith('#'):continue
            v=l.strip().split(maxsplit=22)
            if len(v)<22:raise ValueError('Truncated domtblout')
            if v[0] not in targets or v[3] not in names:raise ValueError('Unexpected target/profile in domtblout')
            m=names[v[3]];L=int(v[2]);hl=int(v[5]);h0,h1=int(v[15]),int(v[16]);a0,a1=int(v[17]),int(v[18])
            if L!=len(targets[v[0]]) or hl!=int(m['length']) or v[4]!=m['accession_version']:raise ValueError('domtbl source identity mismatch')
            if not(1<=h0<=h1<=hl and 1<=a0<=a1<=L):raise ValueError('Invalid domain coordinates')
            cov=(h1-h0+1)/hl;fullE=float(v[6]);inde=float(v[12]);tier='REPORTED_ONLY'
            if scan=='GA':tier='RT_CORE_GA' if m['role']=='RT_CORE' else 'RT_LIKE_AUXILIARY_GA'
            elif fullE<=policy['full_evalue_max'] and inde<=policy['independent_domain_evalue_max'] and cov>=policy['hmm_coverage_min']:tier='RT_CORE_SENSITIVITY' if m['role']=='RT_CORE' else 'RT_LIKE_AUXILIARY_SENSITIVITY'
            rows.append({'target_id':v[0],'target_length':L,'profile_name':v[3],'profile_accession_version':v[4],'profile_role':m['role'],'hmm_length':hl,'full_evalue':fullE,'full_bitscore':float(v[7]),'domain_index':int(v[9]),'domain_count':int(v[10]),'conditional_domain_evalue':float(v[11]),'independent_domain_evalue':inde,'domain_bitscore':float(v[13]),'hmm_from1':h0,'hmm_to1':h1,'ali_from1':a0,'ali_to1':a1,'env_from1':int(v[19]),'env_to1':int(v[20]),'hmm_coverage':cov,'mean_alignment_posterior':float(v[21]),'source_scan':scan,'evidence_tier':tier})
    return rows

def run_hmmer(root,cpus,targets,lock_hash):
    root=Path(root);catalog=read_tsv(root/'refs/MODEL_CATALOG.tsv');model=root/'refs/RT_PANEL.hmm';target=root/'outputs/RT_SCAN_TARGETS.faa';raw=root/'outputs/hmmer';raw.mkdir(exist_ok=True);allrows=[]
    for mode in ('GA','SENSITIVITY'):
        report=raw/(mode+'.report.txt');dom=raw/(mode+'.domtbl');tbl=raw/(mode+'.tbl');done=raw/(mode+'.COMPLETE.json')
        fp={'run_lock':lock_hash,'models_sha256':sha_file(model),'targets_sha256':sha_file(target),'mode':mode}
        valid=False
        if done.exists():
            j=read_json(done);valid=j.get('fingerprint')==fp and all(p.is_file() and j.get('outputs',{}).get(p.name)==sha_file(p) for p in (report,dom,tbl))
        if not valid:
            cmd=['hmmsearch','--cpu',str(cpus),'--seed','42','--noali','-o',report,'--domtblout',dom,'--tblout',tbl]
            cmd+=['--cut_ga'] if mode=='GA' else ['-E','10','--domE','10']
            cmd += [model,target]
            print('[RT HMMER]',mode,'targets=',len(targets),flush=True);run_cmd(cmd,root/'logs'/('hmmer_'+mode+'.log'))
        text=report.read_text();seen=re.findall(r'^Query:\s+(\S+)\s+\[M=\d+\]',text,re.M)
        if sorted(seen)!=sorted(r['name'] for r in catalog) or '[ok]' not in text:raise RuntimeError('HMMER did not complete all fixed profile queries: '+mode)
        rr=parse_domtbl(dom,catalog,targets,mode);allrows.extend(rr)
        write_json(done,{'fingerprint':fp,'outputs':{p.name:sha_file(p) for p in (report,dom,tbl)},'reported_domain_rows':len(rr)})
    # Primary and secondary raw scans remain separate, no post-result threshold selection.
    allrows.sort(key=lambda r:(r['target_id'],r['source_scan'],r['profile_accession_version'],r['ali_from1'],r['ali_to1']))
    write_tsv(root/'outputs/RT_ALL_REPORTED_DOMAINS.tsv',allrows,DOMAIN_FIELDS)
    return allrows
