#!/usr/bin/env python3
"""New-to-project phage source acquisition. No protein, truth, or score inputs."""
from __future__ import annotations
import argparse, csv, datetime as dt, hashlib, json, os, re, shutil, subprocess, sys, time, zipfile
from pathlib import Path
from collections import Counter

ROOT=Path(__file__).resolve().parents[1]

class SourceIneligible(ValueError):
    pass

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def shab(b): return hashlib.sha256(b).hexdigest()
def dump(p,obj):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    q=p.with_suffix(p.suffix+'.part'); q.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n'); q.replace(p)
def rows(p):
    with Path(p).open(newline='',encoding='utf-8-sig') as f: return list(csv.DictReader(f,delimiter='\t'))
def tsv(p,data,fields):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(data)
def field(o,camel,snake=None,default=None):
    a=o.get(camel); b=o.get(snake) if snake else None
    if a is not None and b is not None and a!=b: raise ValueError('conflicting aliases: '+camel)
    return a if a is not None else b if b is not None else default
def accroot(a): return a.split('.')[0]
def reverse_complement(s): return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def parse_fasta(text):
    data=[]; name=None; seq=[]
    for line in text.splitlines():
        line=line.strip()
        if not line: continue
        if line.startswith('>'):
            if name is not None: data.append((name,''.join(seq).upper()))
            name=line[1:].split()[0];seq=[]
        else:
            if name is None: raise ValueError('sequence before header')
            seq.append(line)
    if name is not None: data.append((name,''.join(seq).upper()))
    if not data: raise ValueError('empty FASTA')
    return data

def report_objects(path):
    # JSONL reports and wrapped JSON are both accepted, without silently skipping errors.
    text=Path(path).read_text()
    try: objects=[json.loads(text)]
    except json.JSONDecodeError: objects=[json.loads(x) for x in text.splitlines() if x.strip()]
    for obj in objects:
        if 'reports' in obj:
            for r in obj['reports']:
                yield r.get('report',r)
        elif 'report' in obj: yield obj['report']
        else: yield obj

def normalize(r,policy):
    acc=r.get('accession','')
    if not re.fullmatch(r'[A-Z][A-Z0-9_]*\d+\.\d+',acc): return None,'INVALID_ACCESSION_VERSION'
    if str(r.get('completeness','')).upper()!='COMPLETE': return None,'NOT_COMPLETE'
    n=int(r.get('length',0)); bucket=next((b for b in policy['length_bins'] if b['min']<=n<b['max']),None)
    if not bucket: return None,'OUTSIDE_LENGTH_BINS'
    virus=r.get('virus') or {}; host=r.get('host') or {}
    vname=field(virus,'organismName','organism_name',default='')
    vt=field(virus,'taxId','tax_id',default=0)
    if not vt: return None,'MISSING_VIRAL_TAXID'
    vline=virus.get('lineage') or []
    if not any(str(field(x,'taxId','tax_id'))=='2731619' or x.get('name')=='Caudoviricetes' for x in vline):
        return None,'CAUDOVIRICETES_LINEAGE_UNCONFIRMED'
    hline=host.get('lineage') or []; ht=field(host,'taxId','tax_id',default=0)
    tids={str(field(x,'taxId','tax_id',default=0)) for x in hline}|{str(ht)}
    if '2' in tids: hs='BACTERIAL_HOST_LINEAGE'
    elif ht or hline: return None,'NONBACTERIAL_OR_UNRESOLVED_REPORTED_HOST'
    elif 'phage' in vname.lower(): hs='MISSING_HOST_PHAGE_NAME_ONLY'
    else: return None,'HOST_UNRESOLVED'
    return dict(accession=acc,length_bp=n,bin=bucket['name'],viral_taxid=str(vt),organism_name=vname,
                host_status=hs,release_date=field(r,'releaseDate','release_date',default=''),
                source_database=field(r,'sourceDatabase','source_database',default=''),
                selection_hash=shab((policy['selection_seed']+'|'+acc).encode())),''

SELECT_FIELDS=['accession','length_bp','bin','viral_taxid','organism_name','host_status','release_date','source_database','selection_hash']

def select(records,policy,excluded):
    valid={}; rejected=[]
    for r in records:
        try: z,reason=normalize(r,policy)
        except (TypeError,ValueError) as e: raise ValueError('metadata schema issue: '+str(e))
        if not z: rejected.append({'accession':r.get('accession',''),'reason':reason});continue
        if accroot(z['accession']) in excluded: rejected.append({'accession':z['accession'],'reason':'PRIOR_OR_REGISTERED_ACCESSION'});continue
        if z['accession'] in valid and valid[z['accession']]!=z: raise ValueError('conflicting duplicate metadata '+z['accession'])
        valid[z['accession']]=z
    # Explicit viral taxon ID representative, not a claim of species / homology disjointness.
    order=sorted(valid.values(),key=lambda z:(0 if z['source_database'].lower()=='refseq' else 1,z['selection_hash']))
    seen=set();representatives=[]
    for z in order:
        if z['viral_taxid'] in seen: rejected.append({'accession':z['accession'],'reason':'SAME_VIRAL_TAXID_REPRESENTATIVE_EXISTS'});continue
        seen.add(z['viral_taxid']);representatives.append(z)
    chosen=[]
    for b in policy['length_bins']:
        pool=sorted([z for z in representatives if z['bin']==b['name']],key=lambda z:z['selection_hash'])
        chosen.extend(pool[:b['target']])
        rejected.extend({'accession':z['accession'],'reason':'OUTSIDE_FROZEN_HASH_QUOTA'} for z in pool[b['target']:])
    return chosen,rejected

def audit_package():
    mf=ROOT/'PACKAGE_MANIFEST.tsv'
    if not mf.exists(): raise RuntimeError('PACKAGE_MANIFEST.tsv missing; do not run an unverified package')
    bad=[r['path'] for r in rows(mf) if not (ROOT/r['path']).is_file() or sha(ROOT/r['path'])!=r['sha256']]
    if bad: raise RuntimeError('package altered: '+','.join(bad))

def acquire_metadata(args):
    audit_package();out=ROOT/'outputs'; inp=ROOT/'inputs';inp.mkdir(exist_ok=True);out.mkdir(exist_ok=True)
    path=inp/'NCBI_CAUDOVIRICETES_SNAPSHOT.jsonl'; lock=out/'METADATA_SNAPSHOT.json'
    if lock.exists():
        j=json.loads(lock.read_text())
        if not path.exists() or sha(path)!=j['sha256']:raise RuntimeError('metadata snapshot changed')
        print('[REUSE] verified frozen metadata snapshot'); return
    if path.exists(): raise RuntimeError('unlocked snapshot exists; preserve it and use a clean package directory')
    cmd=['datasets','summary','virus','genome','taxon','Caudoviricetes','--complete-only','--as-json-lines']
    tmp=path.with_suffix('.part')
    if args.metadata: shutil.copyfile(Path(args.metadata),tmp);method='USER_PROVIDED_METADATA_SNAPSHOT';version='not fetched by this run'
    else:
        version=subprocess.check_output(['datasets','version'],text=True).strip()
        (ROOT/'logs').mkdir(exist_ok=True)
        with tmp.open('wb') as f, (ROOT/'logs/metadata.stderr.log').open('wb') as err:
            cp=subprocess.run(cmd,stdout=f,stderr=err)
        if cp.returncode:raise RuntimeError('datasets metadata request failed; see logs; .part retained')
        method='NCBI_DATASETS_CLI'
    count=sum(1 for _ in report_objects(tmp))
    if not count:raise RuntimeError('empty metadata response')
    tmp.replace(path)
    dump(lock,{'method':method,'command':cmd,'datasets_version':version,'utc':dt.datetime.now(dt.timezone.utc).isoformat(),
               'sha256':sha(path),'n_records':count,'protein_annotation_used_for_selection':False})
    print('[PASS] metadata snapshot n='+str(count))

def selection(args):
    acquire_metadata(args);pol=json.loads((ROOT/'config/POLICY.json').read_text());out=ROOT/'outputs'
    lock=out/'SELECTION_FREEZE.json';md=ROOT/'inputs/NCBI_CAUDOVIRICETES_SNAPSHOT.jsonl'; dest=out/'SELECTED_SOURCE_ACCESSIONS.tsv'
    if lock.exists():
        j=json.loads(lock.read_text())
        if sha(md)!=j['metadata_sha256'] or sha(dest)!=j['selection_sha256'] or sha(ROOT/'config/POLICY.json')!=j['policy_sha256']: raise RuntimeError('selection freeze changed')
        print('[REUSE] selection frozen n='+str(j['n_selected']));return
    excluded={accroot(x['accession']) for x in rows(ROOT/'config/EXCLUDE_PRIOR240.tsv')}
    excluded.update(x.strip() for x in (ROOT/'config/EXCLUDE_REGISTERED_REFERENCES.txt').read_text().splitlines() if x.strip())
    chosen,rejected=select(report_objects(md),pol,excluded)
    tsv(dest,chosen,SELECT_FIELDS);tsv(out/'METADATA_EXCLUSION_LEDGER.tsv',rejected,['accession','reason'])
    (out/'SELECTED_ACCESSIONS.txt').write_text(''.join(z['accession']+'\n' for z in chosen))
    counts=Counter(z['bin'] for z in chosen)
    j={'gate':'PASS_METADATA_SELECTION_FROZEN' if chosen else 'HOLD_NO_ELIGIBLE_SOURCES','gate_pass':bool(chosen),
       'n_selected':len(chosen),'target_n':pol['target_n'],'by_bin':dict(counts),
       'quota_shortfall':{b['name']:b['target']-counts[b['name']] for b in pol['length_bins']},
       'metadata_sha256':sha(md),'selection_sha256':sha(dest),'policy_sha256':sha(ROOT/'config/POLICY.json'),
       'one_per_viral_taxid_not_species_independence':True,'new_to_project_not_novel_system':True,'math_scores_observed':False}
    dump(lock,j);print(json.dumps(j,indent=2))
    if not chosen:raise RuntimeError('no eligible sources')

def zip_genome(p,expected):
    with zipfile.ZipFile(p) as z:
        if z.testzip():raise ValueError('ZIP CRC failure')
        names=z.namelist(); fs=[n for n in names if n.endswith('/genomic.fna')]
        if len(fs)!=1:raise ValueError('expected single genomic.fna member')
        # Validate all archive-supplied MD5 lines using exact relative members.
        for mn in [x for x in names if x.endswith('md5sum.txt')]:
            parent=str(Path(mn).parent)
            for line in z.read(mn).decode().splitlines():
                if not line.strip():continue
                hh,fn=line.split(None,1);fn=fn.strip().lstrip('*');fn=fn[2:] if fn.startswith('./') else fn
                key=fn if fn in names else str(Path(parent)/fn)
                if key not in names or hashlib.md5(z.read(key)).hexdigest()!=hh:raise ValueError('member MD5 mismatch '+key)
        rr=parse_fasta(z.read(fs[0]).decode('ascii'))
        if len(rr)!=1:raise SourceIneligible('MULTIRECORD_NOT_ACCEPTED_NO_CONCATENATION')
        name,s=rr[0]
        if expected['accession'] not in name.split('|'):raise ValueError('accession version mismatch')
        if not s or set(s)-set('ACGT'):raise SourceIneligible('NON_ACGT_NO_SILENT_BASE_REPAIR')
        if len(s)!=int(expected['length_bp']):raise ValueError('metadata length mismatch')
        return s

def window_n(n):
    total=0
    for l in [128,256,512,1024,2048,4096]:
        if n<l:continue
        d=l//4;q,r=divmod(n-l,d);total+=2*(q+1+bool(r))
    return total

def fetch_and_freeze(args):
    selection(args);out=ROOT/'outputs';pol=json.loads((ROOT/'config/POLICY.json').read_text())
    chosen=rows(out/'SELECTED_SOURCE_ACCESSIONS.tsv');old=rows(ROOT/'config/EXCLUDE_PRIOR240.tsv')
    oldhash={r[k] for r in old for k in ['sequence_sha256','reverse_complement_sha256']}
    genome_dir=ROOT/'inputs/genomes';genome_dir.mkdir(parents=True,exist_ok=True)
    zdir=ROOT/'work/source_zips';zdir.mkdir(parents=True,exist_ok=True);(ROOT/'logs').mkdir(exist_ok=True)
    frozen=out/'SOURCE_COHORT_FREEZE.json'
    if frozen.exists() and json.loads(frozen.read_text()).get('gate_pass'):
        j=json.loads(frozen.read_text());man=out/'FROZEN_NEW_GENOME_MANIFEST.tsv'
        if sha(man)!=j['manifest_sha256']:raise RuntimeError('frozen manifest altered')
        for r in rows(man):
            if sha(ROOT/r['fasta_path'])!=r['fasta_sha256']:raise RuntimeError('frozen FASTA altered')
        print('[REUSE] frozen source cohort verified');return
    accepted=[];ledger=[];seen={};network_errors=[]
    for i,r in enumerate(chosen,1):
        a=r['accession'];zp=zdir/(a+'.zip');err='';s=None;ineligible='' 
        if zp.exists():
            try:s=zip_genome(zp,r)
            except SourceIneligible as e:ineligible=str(e)
            except Exception as e:err=str(e);raise RuntimeError('existing cache invalid; preserve for inspection '+str(zp)+' '+err)
        else:
            for trial in range(1,int(pol['download_retries'])+1):
                part=zdir/(a+'.part.zip')
                if part.exists():part.unlink()
                cmd=['datasets','download','virus','genome','accession',a,'--include','genome','--filename',str(part)]
                try:
                    with (ROOT/'logs'/('fetch_'+a+'.log')).open('a') as log:
                        cp=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=600)
                    if cp.returncode:raise RuntimeError('datasets rc='+str(cp.returncode))
                    s=zip_genome(part,r);part.replace(zp);break
                except SourceIneligible as e:
                    ineligible=str(e);part.replace(zp);break
                except Exception as e:
                    err=str(e);time.sleep(min(3*trial,10))
            if s is None and not ineligible:
                network_errors.append(a);ledger.append(dict(accession=a,status='UNRESOLVED_DOWNLOAD_OR_SEQUENCE_QC',detail=err));print('[HOLD]',a,err,flush=True);continue
        if ineligible:
            ledger.append(dict(accession=a,status='EXCLUDED_PREDECLARED_SEQUENCE_QC',detail=ineligible));continue
        h=shab(s.encode());rh=shab(reverse_complement(s).encode());canon=min(h,rh)
        if h in oldhash or rh in oldhash:
            ledger.append(dict(accession=a,status='EXCLUDED_EXACT_OR_RC_PRIOR240',detail=h));continue
        if canon in seen:
            ledger.append(dict(accession=a,status='EXCLUDED_EXACT_OR_RC_WITHIN_NEW',detail=seen[canon]));continue
        seen[canon]=a
        fp=genome_dir/(a+'.frozen.fasta');payload='>'+a+'\n'+'\n'.join(s[k:k+80] for k in range(0,len(s),80))+'\n'
        if fp.exists() and fp.read_text()!=payload:raise RuntimeError('existing source FASTA changed '+a)
        fp.write_text(payload)
        accepted.append(dict(r,sequence_sha256=h,reverse_complement_sha256=rh,fasta_path=str(fp.relative_to(ROOT)),fasta_sha256=sha(fp),zip_sha256=sha(zp),oriented_windows=window_n(len(s))))
        ledger.append(dict(accession=a,status='FROZEN_ELIGIBLE',detail=''));print(f'[{i}/{len(chosen)}] FREEZE {a} bp={len(s)}',flush=True)
    mf=out/'FROZEN_NEW_GENOME_MANIFEST.tsv'
    tsv(mf,accepted,SELECT_FIELDS+['sequence_sha256','reverse_complement_sha256','fasta_path','fasta_sha256','zip_sha256','oriented_windows'])
    tsv(out/'SOURCE_ACQUISITION_LEDGER.tsv',ledger,['accession','status','detail'])
    good=bool(accepted) and not network_errors
    j={'stage':'Stage05A','gate':'PASS_NEW_SOURCE_COHORT_FROZEN_WITH_ACTUAL_DENOMINATOR' if good else 'HOLD_SOURCE_ACQUISITION_OR_QC_INCOMPLETE',
       'gate_pass':good,'n_metadata_selected':len(chosen),'n_frozen_new_genomes':len(accepted),'target_not_guaranteed':pol['target_n'],
       'n_unresolved':len(network_errors),'unresolved_accessions':network_errors,'no_replacement_after_sequence_exclusion':True,
       'total_bp':sum(int(r['length_bp']) for r in accepted),'oriented_windows':sum(r['oriented_windows'] for r in accepted),
       'manifest_sha256':sha(mf),'selection_sha256':sha(out/'SELECTED_SOURCE_ACCESSIONS.tsv'),
       'novel_ART_systems_claimed':0,'protein_annotation_used_for_selection':False,'mathematical_scoring_performed':False,
       'UG27_touched':False,'primary_383_unchanged':True,'near_homology_independence_established':False,
       'next':'Resource-informed isolated full-genome scoring; freeze score board before protein/array annotation.'}
    dump(frozen,j)
    plan={'n_genomes':len(accepted),'total_oriented_windows':j['oriented_windows'],
          'estimated_aggregate_feature_bytes':j['oriented_windows']*4*(5468+3570),
          'largest_genome_bp':max([int(r['length_bp']) for r in accepted],default=0),
          'not_peak_RAM':True,'do_not_start_unreviewed_thousand_genome_compute':True,
          'selection_cannot_change_based_on_detector_output':True}
    dump(out/'COMPUTE_RESOURCE_PLAN.json',plan)
    print(json.dumps(j,indent=2));print(json.dumps(plan,indent=2))
    if not good:raise RuntimeError('source cohort not complete; rerun resumes valid cached records')

def collect(args):
    import tarfile
    out=ROOT/'outputs'; items=[]
    for folder in ['outputs','logs','config','docs','evidence','scripts','tests']:
        items.extend(p for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    for p in ROOT.iterdir():
        if p.is_file() and p.suffix in ['.md','.sh','.tsv']:items.append(p)
    # Include metadata snapshot; omit genome and ZIP bulk unless --with-genomes is explicit.
    if (ROOT/'inputs/NCBI_CAUDOVIRICETES_SNAPSHOT.jsonl').exists():items.append(ROOT/'inputs/NCBI_CAUDOVIRICETES_SNAPSHOT.jsonl')
    if args.with_genomes:items.extend(p for p in (ROOT/'inputs/genomes').glob('*.fasta'))
    items=sorted(set(items));stamp=dt.datetime.now().strftime('%Y%m%d_%H%M%S')
    pkg=ROOT/('ART_STAGE05A_SOURCE_RESULTS_'+stamp+'.tar.gz')
    with tarfile.open(pkg,'w:gz') as tf:
        for p in items:tf.add(p,arcname=str(p.relative_to(ROOT)),recursive=False)
    (Path(str(pkg)+'.sha256')).write_text(sha(pkg)+'  '+pkg.name+'\n');print(pkg);print(sha(pkg))

def main():
    p=argparse.ArgumentParser();p.add_argument('step',choices=['metadata','select','acquire','collect']);p.add_argument('--metadata');p.add_argument('--with-genomes',action='store_true');a=p.parse_args()
    lock=ROOT/'.RUNNING'
    if a.step!='collect':
        try:fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
        except FileExistsError:raise SystemExit('A source job is running, or stale .RUNNING exists. Check PID before removing only this lock.')
        os.write(fd,str(os.getpid()).encode());os.close(fd)
    try:{'metadata':acquire_metadata,'select':selection,'acquire':fetch_and_freeze,'collect':collect}[a.step](a)
    finally:
        if a.step!='collect':lock.unlink(missing_ok=True)
if __name__=='__main__':
    try:main()
    except Exception as e:print('[FAIL]',type(e).__name__,str(e),file=sys.stderr);sys.exit(2)
