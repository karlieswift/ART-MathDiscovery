#!/usr/bin/env python3
"""Offline, accession-linked official protein provenance audit.
No inferred translation is used as the primary sequence. All 383 mathematical
candidates survive, including the one genome-only annotation gap.
"""
from pathlib import Path
from collections import defaultdict,Counter
import argparse,csv,hashlib,json,re,sys,zipfile
from io_utils import *
HERE=Path(__file__).resolve().parents[1]
POL=read_json(HERE/'config/POLICY.json')
RT=re.compile(r'reverse transcript|rna[- ]directed dna polymerase|\bretron\b|group ii intron|\bmaturase\b|diversity[- ]generating retro|\bdgr\b',re.I)
TRANS=str.maketrans('ACGTNRYKMSWBDHV','TGCANYRMKSWVHDB')
BASE='ncbi_dataset/data/'

def rc(x):return x.translate(TRANS)[::-1]

def get_alias(x,*keys,default=None):
    for k in keys:
        if k in x:return x[k]
    return default

def parse_report(data,acc):
    rows=[json.loads(line) for line in data.decode('utf-8').splitlines() if line.strip()]
    matches=[r for r in rows if r.get('accession')==acc]
    if len(matches)!=1:raise ValueError(f'{acc}: expected exactly one report, got {len(matches)}')
    return matches[0]

def iter_peptide_nodes(annotation):
    """Traverse only schema-declared CDS and maturePeptide branches.

    A processed peptide is NOT another independent gene. Keep parentage, and
    still demand explicit source identifiers and exact sequence linkage.
    """
    def walk(node,path,kind,gene,parent=''):
        if not isinstance(node,dict):
            raise ValueError('Non-object peptide node: '+path)
        yield node,path,kind,gene,parent
        kids=get_alias(node,'maturePeptide','mature_peptide',default=[])
        if not isinstance(kids,list):
            raise ValueError('Non-list maturePeptide: '+path)
        for i,kid in enumerate(kids):
            yield from walk(kid,path+'.maturePeptide['+str(i)+']',
                            'MATURE_PEPTIDE',gene,path)
    genes=annotation.get('genes',[])
    if not isinstance(genes,list):raise ValueError('Non-list genes')
    for gi,gene in enumerate(genes):
        for ci,node in enumerate(gene.get('cds',[])):
            path='genes['+str(gi)+'].cds['+str(ci)+']'
            yield from walk(node,path,'CDS',gene)

def declared_coding_strand(node,gene):
    """Do not infer protein coding direction from cds.fna payload direction.

    Use CDS/peptide metadata first, then explicit gene metadata. Missing is '?'_
    (without the underscore), not silently '+'.
    """
    def value(info):
        strands=set()
        for r in info.get('range',[]):
            o=r.get('orientation')
            if o in ('plus','PLUS','+',1,'1'):strands.add('+')
            elif o in ('minus','MINUS','-',2,'2'):strands.add('-')
            elif o not in (None,'none','NONE','',0,'0'):
                raise ValueError('Unknown orientation enum: '+str(o))
        if len(strands)>1:return '?'
        return next(iter(strands)) if strands else '?'
    v=value(node.get('nucleotide',{}))
    if v!='?':return v,'PEPTIDE_NUCLEOTIDE_METADATA'
    v=value(gene.get('nucleotide',{}))
    return (v,'GENE_NUCLEOTIDE_METADATA') if v!='?' else ('?','NOT_DECLARED')

def fasta_reference_inventory(obj,path='$'):
    """Diagnostic index, never a fallback that invents a feature link."""
    if isinstance(obj,dict):
        for k in ('seqId','seq_id','accessionVersion','accession_version','accession'):
            if k in obj and isinstance(obj[k],(str,int)):
                yield dict(json_path=path+'.'+k,value=str(obj[k]))
        for k,v in obj.items():yield from fasta_reference_inventory(v,path+'.'+k)
    elif isinstance(obj,list):
        for i,v in enumerate(obj):yield from fasta_reference_inventory(v,path+'['+str(i)+']')

def aa_difference_class(old,official):
    if old==official:return 'IDENTICAL'
    if len(old)==len(official) and len(old)>0 and old[1:]==official[1:]:
        return 'FIRST_RESIDUE_ONLY'
    if len(old)==len(official)+1 and old[:-1]==official:
        return 'ONE_EXTRA_TERMINAL_RESIDUE_ONLY'
    if len(official)==len(old)+1 and official[:-1]==old:
        return 'ONE_MISSING_TERMINAL_RESIDUE_ONLY'
    return 'NOT_FIRST_OR_TERMINAL_ONLY'

def select_record(info,records,kind):
    exact=get_alias(info,'seqId','seq_id',default='')
    if exact:
        if exact not in records:raise ValueError(f'{kind} metadata seqId missing from FASTA: {exact}')
        return exact,records[exact]
    accession=get_alias(info,'accessionVersion','accession_version',default='')
    if not accession:raise ValueError(f'{kind} missing seqId/accessionVersion')
    ids=[k for k in records if k.split(':')[0]==accession]
    if len(ids)!=1:raise ValueError(f'{kind} accession cannot uniquely resolve a FASTA record: {accession}')
    return ids[0],records[ids[0]]

def parse_ranges(nuc,acc,n):
    owner=get_alias(nuc,'accessionVersion','accession_version',default='')
    seqid=get_alias(nuc,'seqId','seq_id',default='')
    if owner and owner!=acc:raise ValueError('nucleotide accession differs from requested frozen genome')
    if not owner and (not seqid or seqid.split(':')[0]!=acc):raise ValueError('no exact genomic accession link')
    ranges=nuc.get('range',[])
    if not isinstance(ranges,list) or not ranges:raise ValueError('missing nucleotide ranges')
    parts=[]
    for i,r in enumerate(ranges):
        a=int(r['begin']);b=int(r['end']);lo=min(a,b)-1;hi=max(a,b)
        if lo<0 or hi>n:raise ValueError('provider CDS range outside exact frozen genome')
        # A slippage field needs a dedicated biological reconstruction, not
        # silent interval concatenation. Preserve it in the issue report.
        if int(get_alias(r,'ribosomalSlippage','ribosomal_slippage',default=0) or 0)!=0:
            raise ValueError('explicit ribosomal slippage: separate reconstruction required')
        ori=r.get('orientation',None)
        if ori in ('plus','PLUS','+',1,'1'):strand='+'
        elif ori in ('minus','MINUS','-',2,'2'):strand='-'
        elif ori in (None,'none','NONE','',0,'0'):strand='?'
        else:raise ValueError('unsupported orientation enum: '+str(ori))
        parts.append(dict(start0=lo,end0=hi,strand=strand,order=r.get('order',None),original_index=i))
    return parts

def exact_cds_mapping(genome,nt,parts):
    """Prove orientation/order by exact nucleotide equality, never header guess.
    Provider order, when fully specified, takes precedence. Missing orientations
    allow two uniform alternatives; ambiguity is reported, not resolved by rank.
    """
    orders={p['order'] for p in parts if p['order'] is not None}
    if len(orders)==len(parts) and all(p['order'] is not None for p in parts):
        variants=[sorted(parts,key=lambda p:int(p['order']))]
    else:
        variants=[parts]
        if len(parts)>1:variants.append(list(reversed(parts)))
    matches={}
    for order in variants:
        for default in ('+','-'):
            oriented=[dict(p,strand=p['strand'] if p['strand']!='?' else default) for p in order]
            rebuilt=''.join(genome[p['start0']:p['end0']] if p['strand']=='+' else rc(genome[p['start0']:p['end0']]) for p in oriented)
            if rebuilt==nt:
                key=tuple((p['start0'],p['end0'],p['strand']) for p in oriented);matches[key]=oriented
    if len(matches)!=1:raise ValueError(f'CDS/genome exact orientation mapping count={len(matches)}')
    result=next(iter(matches.values()))
    strands={p['strand'] for p in result}
    if len(strands)!=1:raise ValueError('mixed-strand compound CDS requires separate audit')
    return result,next(iter(strands))

def near_parts(parts,c0,c1,radius,n):
    lo=max(0,c0-radius);hi=min(n,c1+radius)
    hit=any(p['start0']<hi and p['end0']>lo for p in parts)
    def dist(p):
        a,b=p['start0'],p['end0']
        return c0-b if b<=c0 else a-c1 if a>=c1 else 0
    return hit,min(dist(p) for p in parts),lo,hi

def protein_status(seq):
    if '*' in seq:return 'OFFICIAL_INTERNAL_OR_TERMINAL_STOP_REVIEW'
    if not seq or set(seq)-set('ACDEFGHIKLMNPQRSTVWYBXZJUO'):return 'OFFICIAL_INVALID_AA_REVIEW'
    return 'OFFICIAL_SEQUENCE'

def freeze_run_inputs(project,a04,r1a,r1b,r1c):
    reg=HERE/'inputs/STAGE04A_R1_DUAL_ARM_REGIONS.tsv'
    check_hash(reg,POL['stage04a_registry_sha256'],'canonical 383 registry')
    candidates=read_tsv(reg)
    if len(candidates)!=383 or len({(r['analysis_id'],r['region_id']) for r in candidates})!=383 or len({r['sequence_id'] for r in candidates})!=174:
        raise RuntimeError('Canonical candidate keys changed')
    stage_manifest=a04/'outputs/STAGE04A_UNIVERSE_MANIFEST.tsv';check_hash(stage_manifest,POL['stage04a_universe_sha256'])
    m=read_tsv(stage_manifest)
    if len(m)!=240 or len({r['sequence_id'] for r in m})!=240:raise RuntimeError('Stage04A universe invalid')
    cache_m=r1a/'outputs/STAGE04B_R1A_DATASETS_CACHE_MANIFEST.tsv';check_hash(cache_m,POL['r1a_manifest_sha256'])
    cm=read_tsv(cache_m)
    if len(cm)!=174 or {r['accession'] for r in cm}!={r['sequence_id'] for r in candidates}:raise RuntimeError('Cache accession set changed')
    r1bg=r1b/'outputs/STAGE04B_R1B_GATE.json';g=read_json(r1bg)
    if g.get('gate_pass') is not True or g.get('n_exact_sequence_equivalent')!=174:raise RuntimeError('Prior R1B gate not PASS 174/174')
    oldtab=r1c/'outputs/stage04b_r1d_freeze/STAGE04B_R1D_NEIGHBORHOOD_PROTEIN_REGISTRY.tsv'
    oldfa=r1c/'work/stage04b_r1d_freeze/STAGE04B_R1D_NEIGHBORHOOD_PROTEINS.faa'
    check_hash(oldtab,POL['r1d0_registry_sha256']);check_hash(oldfa,POL['r1d0_fasta_sha256'])
    old=read_tsv(oldtab);oldseq=fasta_text(oldfa.read_text())
    if len(old)!=9560 or len(oldseq)!=9560:raise RuntimeError('R1D0 reported 9560 records not reproduced')
    for r in old:
        if r['protein_key'] not in oldseq or shatext(oldseq[r['protein_key']])!=r['aa_sha256']:raise RuntimeError('R1D0 FASTA/registry mismatch')
    paths=[reg,stage_manifest,cache_m,r1bg,oldtab,oldfa]
    locks=[dict(path=str(p),bytes=p.stat().st_size,sha256=digest(p)) for p in paths]
    write_tsv(HERE/'outputs/INPUT_LOCK.tsv',locks,['path','bytes','sha256'])
    return candidates,{r['sequence_id']:r for r in m},{r['accession']:r for r in cm},old,oldseq,locks

def analyze(a04,r1a,r1b,r1c):
    candidates,genome_m,cache_m,old,oldseq,locks=freeze_run_inputs(a04.parent,a04,r1a,r1b,r1c)
    out=HERE/'outputs';issues=[];provenance=[];features=[];links=[];availability=[];domains=[]
    reference_index=[];traversal=[];unlinked_inventory=[];source_zips={}
    oldsummary=r1c/'outputs/STAGE04B_R1C_383_CANDIDATE_NEIGHBORHOOD_SUMMARY.tsv'
    old_near=read_tsv(oldsummary) if oldsummary.is_file() else []
    oldkeys={ (r['sequence_id'],r['region_id']):r for r in old_near}
    byacc=defaultdict(list)
    for c in candidates:byacc[c['sequence_id']].append(c)
    genomes={};status={};duplicate_n=0
    for i,(acc,cc) in enumerate(sorted(byacc.items()),1):
        print(f'[{i}/174] {acc}',flush=True)
        mr=genome_m[acc];frozen=(a04/mr['frozen_fasta']).resolve()
        if not frozen.is_relative_to(a04):raise RuntimeError('Frozen input escapes Stage04A root')
        check_hash(frozen,mr['fasta_sha256']);_,genome=single_fasta(frozen)
        zpath=r1a/'work/zips'/f'{acc}.zip'
        # Cache archives were also placed under R1B by the user's cache import.
        if not zpath.is_file():zpath=r1b/'work/zips'/f'{acc}.zip'
        check_hash(zpath,cache_m[acc]['zip_sha256'],'frozen Datasets ZIP')
        source_zips[acc]=zpath
        with zipfile.ZipFile(zpath) as z:
            if z.testzip() is not None:raise RuntimeError('ZIP CRC failed: '+acc)
            names=set(z.namelist());gb=z.read(BASE+'genomic.fna')
            if hashlib.sha256(gb).hexdigest()!=cache_m[acc]['genomic_fna_sha256']:raise RuntimeError('Manifest/ZIP genomic bytes mismatch')
            fna=fasta_text(gb.decode())
            if len(fna)!=1 or acc not in fna or fna[acc]!=genome:raise RuntimeError('Frozen/Datasets exact genomic sequence mismatch: '+acc)
            n=len(genome);genomes[acc]=genome
            if any(int(c['length_bp'])!=n or not 0<=int(c['start0'])<int(c['end0'])<=n for c in cc):raise RuntimeError('Candidate coordinates out of bounds')
            report=parse_report(z.read(BASE+'data_report.jsonl'),acc)
            ann_data=z.read(BASE+'annotation_report.jsonl');ann=parse_report(ann_data,acc)
            reference_index.extend(dict(accession=acc,**r) for r in fasta_reference_inventory(ann))
            cds=fasta_text(z.read(BASE+'cds.fna').decode()) if BASE+'cds.fna' in names else {}
            proteins=fasta_text(z.read(BASE+'protein.faa').decode()) if BASE+'protein.faa' in names else {}
            provenance.append(dict(accession=acc,genome_sha256=shatext(genome),frozen_file_sha256=mr['fasta_sha256'],
                                   zip_sha256=digest(zpath),annotation_json_sha256=hashlib.sha256(ann_data).hexdigest(),
                                   cds_file_sha256=hashlib.sha256(z.read(BASE+'cds.fna')).hexdigest() if cds else '',
                                   protein_file_sha256=hashlib.sha256(z.read(BASE+'protein.faa')).hexdigest() if proteins else ''))
        genes=ann.get('genes',[])
        virus=report.get('virus',{})
        org=get_alias(virus,'organismName','organism_name',default='')
        if not cds and not proteins and not genes:
            status[acc]='GENOME_ONLY_NO_OFFICIAL_CDS'
            if acc not in POL['genome_only_accessions']:issues.append(dict(accession=acc,kind='NEW_GENOME_ONLY_GAP',detail='not the frozen declared gap'))
            availability.append(dict(accession=acc,status=status[acc],n_cds_fasta=0,n_protein_fasta=0,n_report_cds=0,n_unlinked_cds=0,organism=org))
            continue
        status[acc]='OFFICIAL_ANNOTATION_AUDITED'
        seen_cds=set();seen_proteins=set();seen_feature={};n_report=0;acc_issues=0;n_mature=0
        for node,where,feature_type,g,parent_path in iter_peptide_nodes(ann):
                n_report+=1
                n_mature+=int(feature_type=='MATURE_PEPTIDE')
                traversal.append(dict(accession=acc,json_path=where,feature_type=feature_type,parent_json_path=parent_path,cds_seq_id=get_alias(node.get('nucleotide',{}),'seqId','seq_id',default=''),protein_seq_id=get_alias(node.get('protein',{}),'seqId','seq_id',default='')))
                try:
                    nuc=node.get('nucleotide',{});pinfo=node.get('protein',{})
                    cid,nt=select_record(nuc,cds,'CDS')
                    pid,aa=select_record(pinfo,proteins,'protein')
                    parts=parse_ranges(nuc,acc,n);mapped,strand=exact_cds_mapping(genome,nt,parts)
                    seen_cds.add(cid);seen_proteins.add(pid)
                    coding_strand,coding_basis=declared_coding_strand(node,g)
                    st=min(p['start0'] for p in mapped);en=max(p['end0'] for p in mapped)
                    fk=acc+'__'+hashlib.sha256((cid+'|'+pid+'|'+json.dumps(mapped,sort_keys=True)).encode()).hexdigest()[:20]
                    product=node.get('name') or get_alias(pinfo,'title',default='')
                    ps=protein_status(aa)
                    near=[]
                    for c in cc:
                        hit,dist,lo,hi=near_parts(mapped,int(c['start0']),int(c['end0']),10000,n)
                        if hit:near.append((c,dist))
                    if fk in seen_feature:
                        prior=seen_feature[fk]
                        if prior!=(nt,aa):raise ValueError('Conflicting duplicate provider feature')
                        duplicate_n+=1;continue
                    seen_feature[fk]=(nt,aa)
                    if not near:continue
                    if ps!='OFFICIAL_SEQUENCE':raise ValueError(ps)
                    keys=sorted(f"{c['sequence_id']}::{c['region_id']}" for c,dist in near)
                    fr=dict(protein_key=fk,sequence_id=acc,official_cds_fasta_id=cid,official_protein_fasta_id=pid,
                            official_protein_accession=get_alias(pinfo,'accessionVersion','accession_version',default=pid.split(':')[0]),
                            cds_start0=st,cds_end0=en,cds_strand=coding_strand,cds_fasta_payload_strand=strand,coding_strand_basis=coding_basis,
                            feature_type=feature_type,parent_json_path=parent_path,cds_parts_json=json.dumps(mapped,separators=(',',':')),
                            cds_nt_length=len(nt),cds_nt_sha256=shatext(nt),product=product,
                            aa_length=len(aa),aa_sha256=shatext(aa),aa_sequence=aa,n_ambiguous_aa=aa.count('X'),
                            candidate_keys=';'.join(keys),rt_related_text_flag=int(bool(RT.search(product))),
                            annotation_json_path=where,official_protein_primary=True)
                    features.append(fr)
                    for c,dist in near:links.append(dict(candidate_key=f"{acc}::{c['region_id']}",sequence_id=acc,region_id=c['region_id'],
                                                         protein_key=fk,distance_to_candidate_bp=dist,official_protein_accession=fr['official_protein_accession'],
                                                         product=product,rt_related_text_flag=fr['rt_related_text_flag']))
                    for domain in node.get('cdd',[]):
                        domains.append(dict(protein_key=fk,sequence_id=acc,candidate_keys=fr['candidate_keys'],source='NCBI_DATASETS_PREEXISTING_CDD_ANNOTATION',
                                            domain_accession=domain.get('accession',''),domain_name=domain.get('name',''),range_json=json.dumps(domain.get('range',{}),sort_keys=True)))
                except Exception as e:
                    acc_issues+=1;issues.append(dict(accession=acc,kind='CDS_OR_PROTEIN_MAPPING_UNRESOLVED',detail=where+': '+str(e)))
        unlinked=set(cds)-seen_cds;unlinked_p=set(proteins)-seen_proteins
        refs=[x for x in reference_index if x['accession']==acc]
        for kind,unmatched,sequences in [('CDS',unlinked,cds),('PROTEIN',unlinked_p,proteins)]:
            for cid in sorted(unmatched):
                issues.append(dict(accession=acc,kind='UNLINKED_'+kind+'_FASTA_RECORD',detail=cid))
                unlinked_inventory.append(dict(accession=acc,record_type=kind,fasta_id=cid,sequence_length=len(sequences[cid]),sequence_sha256=shatext(sequences[cid]),exact_json_reference_paths=';'.join(x['json_path'] for x in refs if x['value']==cid)))
        if unlinked or unlinked_p or acc_issues:status[acc]='ANNOTATION_LINKAGE_INCOMPLETE'
        availability.append(dict(accession=acc,status=status[acc],n_cds_fasta=len(cds),n_protein_fasta=len(proteins),n_report_cds=n_report-n_mature,n_report_mature_peptides=n_mature,n_report_peptide_nodes=n_report,n_unlinked_cds=len(unlinked),n_unlinked_proteins=len(unlinked_p),organism=org))
    features.sort(key=lambda x:(x['sequence_id'],x['cds_start0'],x['cds_end0'],x['protein_key']))
    bycandidate=defaultdict(list)
    for x in links:bycandidate[x['candidate_key']].append(x)
    summaries=[];changes=[]
    for c in candidates:
        acc=c['sequence_id'];key=f"{acc}::{c['region_id']}";ll=bycandidate[key];n=len(genomes[acc]);state=status[acc]
        n0=max(0,int(c['start0'])-10000);n1=min(n,int(c['end0'])+10000)
        assessed=state=='OFFICIAL_ANNOTATION_AUDITED'
        summary=dict(**c,candidate_key=key,annotation_status=state,neighborhood_start0=n0,neighborhood_end0=n1,
                     n_official_neighborhood_proteins=len(ll),n_rt_related_text_proteins=sum(x['rt_related_text_flag'] for x in ll) if assessed else '',
                     rt_related_text_assessed=assessed,rt_absence_inferred=False,genome_only_orf_rescue_pending=acc in POL['genome_only_accessions'])
        summaries.append(summary)
        prev=oldkeys.get((acc,c['region_id']),{})
        changes.append(dict(candidate_key=key,old_neighborhood_end0=prev.get('neighborhood_end0',''),new_neighborhood_end0=n1,
                            old_nearby_cds=prev.get('n_nearby_cds',''),new_official_proteins=len(ll),new_status=state))
    # Compare all 9,560 former translations to authoritative sequences by
    # accession + locus, not by file order. Strand disagreements remain visible.
    exactloc=defaultdict(list);loose=defaultdict(list)
    for f in features:
        exactloc[(f['sequence_id'],f['cds_start0'],f['cds_end0'],f['cds_strand'])].append(f)
        loose[(f['sequence_id'],f['cds_start0'],f['cds_end0'])].append(f)
    comparison=[]
    for r in old:
        key=(r['sequence_id'],int(r['cds_start0']),int(r['cds_end0']),r['cds_strand'])
        hits=loose.get(key[:3],[]);orientation_mismatch=False
        if len(hits)==1:
            orientation_mismatch=hits[0]['cds_strand'] in ('+','-') and hits[0]['cds_strand']!=r['cds_strand']
        if len(hits)!=1:
            state='UNRESOLVED_LOCUS' if not hits else 'AMBIGUOUS_LOCUS'
            f=None
        else:
            f=hits[0]
            state='AA_IDENTICAL' if r['aa_sha256']==f['aa_sha256'] else 'AA_DIFFERENT'
        comparison.append(dict(old_protein_key=r['protein_key'],sequence_id=r['sequence_id'],status=state,
                               old_aa_length=r['aa_length'],old_aa_sha256=r['aa_sha256'],old_X_count=oldseq[r['protein_key']].count('X'),
                               official_protein_key=f['protein_key'] if f else '',official_protein_accession=f['official_protein_accession'] if f else '',
                               official_aa_length=f['aa_length'] if f else '',official_aa_sha256=f['aa_sha256'] if f else '',
                               aa_difference_class=aa_difference_class(oldseq[r['protein_key']],f['aa_sequence']) if f else 'UNRESOLVED',
                               old_strand_field=r['cds_strand'],official_coding_strand=f['cds_strand'] if f else '',
                               official_cds_payload_strand=f['cds_fasta_payload_strand'] if f else '',coding_strand_differs_from_old_field=orientation_mismatch))
    ffields=['protein_key','sequence_id','official_cds_fasta_id','official_protein_fasta_id','official_protein_accession','cds_start0','cds_end0','cds_strand','cds_fasta_payload_strand','coding_strand_basis','feature_type','parent_json_path','cds_parts_json','cds_nt_length','cds_nt_sha256','product','aa_length','aa_sha256','n_ambiguous_aa','candidate_keys','rt_related_text_flag','annotation_json_path','official_protein_primary']
    write_tsv(out/'OFFICIAL_NEIGHBORHOOD_PROTEIN_REGISTRY.tsv',features,ffields)
    with (out/'OFFICIAL_NEIGHBORHOOD_PROTEINS.faa').open('w') as f:
        for r in features:
            f.write('>'+r['protein_key']+'\n');s=r['aa_sequence']
            for j in range(0,len(s),80):f.write(s[j:j+80]+'\n')
    write_tsv(out/'CANDIDATE_PROTEIN_LINKS.tsv',links,['candidate_key','sequence_id','region_id','protein_key','distance_to_candidate_bp','official_protein_accession','product','rt_related_text_flag'])
    write_tsv(out/'ALL_383_CANDIDATES.tsv',summaries,list(summaries[0]))
    write_tsv(out/'R1D0_VS_OFFICIAL_PROTEIN_AUDIT.tsv',comparison,list(comparison[0]))
    write_tsv(out/'R1C_NEIGHBORHOOD_AUDIT.tsv',changes,list(changes[0]))
    write_tsv(out/'ANNOTATION_COVERAGE.tsv',availability,['accession','status','n_cds_fasta','n_protein_fasta','n_report_cds','n_report_mature_peptides','n_report_peptide_nodes','n_unlinked_cds','n_unlinked_proteins','organism'])
    write_tsv(out/'REPORT_NODE_TRAVERSAL.tsv',traversal,['accession','json_path','feature_type','parent_json_path','cds_seq_id','protein_seq_id'])
    write_tsv(out/'ALL_JSON_FASTA_REFERENCES.tsv',reference_index,['accession','json_path','value'])
    write_tsv(out/'UNLINKED_FASTA_INVENTORY.tsv',unlinked_inventory,['accession','record_type','fasta_id','sequence_length','sequence_sha256','exact_json_reference_paths'])
    write_tsv(out/'SOURCE_PROVENANCE.tsv',provenance,list(provenance[0]))
    write_tsv(out/'ANNOTATION_MAPPING_ISSUES.tsv',issues,['accession','kind','detail'])
    write_tsv(out/'PREEXISTING_CDD_ANNOTATIONS.tsv',domains,['protein_key','sequence_id','candidate_keys','source','domain_accession','domain_name','range_json'])
    rtlinks=[x for x in links if x['rt_related_text_flag']]
    write_tsv(out/'RT_RELATED_TEXT_DETAIL.tsv',rtlinks,['candidate_key','sequence_id','region_id','protein_key','distance_to_candidate_bp','official_protein_accession','product','rt_related_text_flag'])
    # Preserve a full-genome ORF-rescue handoff, with the fixed neighborhood
    # registered separately. No gene predictor or new detector is run here.
    rescue=out/'genome_only_rescue';rescue.mkdir(exist_ok=True)
    orf_rows=[]
    for c in summaries:
        if not c['genome_only_orf_rescue_pending']:continue
        acc=c['sequence_id'];s=genomes[acc];p=rescue/(acc+'.whole_genome.fasta')
        with p.open('w') as f:
            f.write('>'+acc+'\n')
            for j in range(0,len(s),80):f.write(s[j:j+80]+'\n')
        orf_rows.append(dict(candidate_key=c['candidate_key'],genome_file=p.name,genome_sequence_sha256=shatext(s),candidate_start0=c['start0'],candidate_end0=c['end0'],neighborhood_start0=c['neighborhood_start0'],neighborhood_end0=c['neighborhood_end0'],status='ORF_PREDICTION_PENDING_NOT_NEGATIVE'))
    write_tsv(rescue/'ORF_RESCUE_HANDOFF.tsv',orf_rows,['candidate_key','genome_file','genome_sequence_sha256','candidate_start0','candidate_end0','neighborhood_start0','neighborhood_end0','status'])
    # Return exact frozen sources for the three previously observed diagnostics
    # plus any new problem accession. Never request/download current replacements.
    evidence_dir=out/'source_evidence';evidence_dir.mkdir(exist_ok=True)
    evidence_accs=set(POL.get('source_evidence_accessions',[]))|{x['accession'] for x in issues}
    import shutil
    for acc in sorted(evidence_accs):
        if acc in source_zips:shutil.copyfile(source_zips[acc],evidence_dir/(acc+'.zip'))
    for p in [r1c/'outputs/stage04b_r1d_freeze/STAGE04B_R1D_NEIGHBORHOOD_PROTEIN_REGISTRY.tsv',r1c/'work/stage04b_r1d_freeze/STAGE04B_R1D_NEIGHBORHOOD_PROTEINS.faa']:
        if p.exists():shutil.copyfile(p,evidence_dir/p.name)
    for l in locks:check_hash(l['path'],l['sha256'],'source after audit')
    counter=Counter(r['status'] for r in comparison)
    formal=(len(summaries)==383 and len(provenance)==174 and not issues and set(k for k,v in status.items() if v=='GENOME_ONLY_NO_OFFICIAL_CDS')==set(POL['genome_only_accessions']))
    # Old missing records are evidence of a mapping discrepancy, even if the
    # new official table is independently complete. Do not hide this in a PASS.
    unresolved=sum(counter[k] for k in ['UNRESOLVED_LOCUS','AMBIGUOUS_LOCUS'])
    formal=formal and unresolved==0
    man=manifest_tree(out,out/'OUTPUT_MANIFEST.tsv',exclude=['GATE_STAGE04B_R1D1.json'])
    gate=dict(stage='Stage04B-R1D1',implementation_version='v1.0.1',gate='PASS_OFFICIAL_PROTEIN_AUDIT_WITH_DECLARED_GENOME_ONLY_GAP' if formal else 'HOLD_OFFICIAL_PROTEIN_LINKAGE_OR_COVERAGE_AUDIT',
              gate_pass=formal,n_candidates_retained=len(summaries),n_exact_genomes_rechecked=len(provenance),
              n_old_translation_records=len(old),n_official_neighborhood_proteins=len(features),
              n_candidates_with_official_proteins=sum(r['n_official_neighborhood_proteins']>0 for r in summaries),old_vs_official_counts=dict(counter),old_comparison_strand_is_separate_from_aa_status=True,
              n_mapping_issues=len(issues),n_identical_duplicate_report_features=duplicate_n,
              n_neighborhood_mature_peptides=sum(f['feature_type']=='MATURE_PEPTIDE' for f in features),
              n_report_mature_peptides=sum(t['feature_type']=='MATURE_PEPTIDE' for t in traversal),
              n_unknown_coding_strand=sum(f['cds_strand']=='?' for f in features),
              aa_difference_classes=dict(Counter(r['aa_difference_class'] for r in comparison)),
              n_coding_strand_disagreements_with_old=sum(r['coding_strand_differs_from_old_field'] for r in comparison),
              fasta_payload_orientation_not_assumed_coding_direction=True,
              n_candidates_with_rt_related_text=len({x['candidate_key'] for x in rtlinks}),
              genome_only_accessions=POL['genome_only_accessions'],orf_rescue_pending=True,
              official_protein_sequences_used_without_retranslation=True,internal_stops_replaced_with_X=False,
              sources_unchanged=True,mathematical_scores_recomputed=False,candidate_filtering_performed=False,
              full_383_candidate_biological_classification_complete=False,domain_scan_performed=False,
              official_protein_domain_scan_ready=formal,negative_calls_authorized=False,
              output_manifest_sha256=man,runtime=runtime_info())
    js(out/'GATE_STAGE04B_R1D1.json',gate);print(json.dumps(gate,indent=2))
    return 0 if formal else 2

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--project-root',type=Path,required=True)
    ap.add_argument('--stage04a-root',type=Path);ap.add_argument('--r1a-root',type=Path);ap.add_argument('--r1b-root',type=Path);ap.add_argument('--r1c-root',type=Path)
    a=ap.parse_args();root=a.project_root.resolve()
    a04=(a.stage04a_root or root/'ART_MATH_DISCOVERY_STAGE04A_L4_BROAD_UNIVERSE_FROZEN_DISCOVERY_v1_0_0_20260928').resolve()
    r1a=(a.r1a_root or root/'ART_MATH_DISCOVERY_STAGE04B_R1A_NCBI_DATASETS_OFFLINE_CACHE_v1_0_0_20261005').resolve()
    r1b=(a.r1b_root or root/'ART_MATH_DISCOVERY_STAGE04B_R1B_EXACT_SEQUENCE_PROVENANCE_20261005').resolve()
    r1c=(a.r1c_root or root/'ART_MATH_DISCOVERY_STAGE04B_R1C_CANDIDATE_NEIGHBORHOOD_ANNOTATION_v1_0_0_20261005').resolve()
    with lock_dir(HERE/'.run_lock'):
        if (HERE/'outputs/GATE_STAGE04B_R1D1.json').exists():
            stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S_%f')
            old=HERE/'prior_runs'/stamp;old.parent.mkdir(exist_ok=True)
            (HERE/'outputs').rename(old);(HERE/'outputs').mkdir()
        return analyze(a04,r1a,r1b,r1c)
if __name__=='__main__':
    try:sys.exit(main())
    except Exception as e:
        js(HERE/'outputs/GATE_STAGE04B_R1D1.json',dict(stage='Stage04B-R1D1',gate='HOLD_PROTEIN_AUDIT_INPUT_OR_RUNTIME',gate_pass=False,error=str(e),domain_scan_authorized=False))
        print('[HOLD]',e,file=sys.stderr);raise SystemExit(2)
