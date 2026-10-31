#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, gzip, hashlib, json, math, os, sys
from collections import defaultdict, Counter

FA_EXT={'.fa','.faa','.fasta','.fas','.fsa'}

def fasta_iter(p):
    p=Path(p)
    op=gzip.open if p.suffix=='.gz' else open
    h=None; seq=[]
    with op(p,'rt') as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if h is not None:
                    yield h,''.join(seq).upper().replace(' ','')
                h=line[1:].strip(); seq=[]
            else:
                seq.append(line)
        if h is not None:
            yield h,''.join(seq).upper().replace(' ','')

def shaseq(seq):
    return hashlib.sha256(seq.encode()).hexdigest()

def shafile(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):
            h.update(b)
    return h.hexdigest()

def read_tsv(p):
    with open(p,newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(p, rows, fields):
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

def pct(a,b):
    return 100.0*a/b if b else None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--stage04a-root',required=True)
    ap.add_argument('--delta-root',required=True)
    ap.add_argument('--out',required=True)
    a=ap.parse_args()
    A=Path(a.stage04a_root).resolve(); D=Path(a.delta_root).resolve(); O=Path(a.out).resolve(); O.mkdir(parents=True,exist_ok=True)

    satp=A/'outputs/STAGE04A_60G_SATURATION_CENSUS.json'
    famp=A/'outputs/COMBINED_FAMILY_SUMMARY.tsv'
    memp=A/'outputs/COMBINED_FAMILY_MEMBERSHIP.tsv'
    map_p=A/'outputs/GLOBAL_EXACT_AA_MAP.tsv'
    manp=D/'frozen/F3V2_DELTA_SELECTED_ACCESSIONS.tsv'
    faa_dir=D/'outputs/formal_faa'
    required=[satp,famp,memp,map_p,manp]
    for p in required:
        if not p.is_file(): raise RuntimeError('missing required input: '+str(p))
    if not faa_dir.is_dir(): raise RuntimeError('missing formal_faa directory: '+str(faa_dir))

    sat=json.load(open(satp))
    if sat.get('gate_pass') is not True or sat.get('gate')!='PASS_F3V2_STAGE04A_60G_COMBINED_FAMILY_CENSUS_COMPLETE':
        raise RuntimeError('Stage04A saturation gate is not PASS')

    famrows=read_tsv(famp)
    famclass={r['family_rep']:r['family_class'] for r in famrows}
    counts=Counter(famclass.values())
    expected_class_counts={
      'OLD_ONLY':int(sat['combined_old_only_family_n']),
      'DELTA_ONLY':int(sat['combined_delta_only_family_n']),
      'MIXED_OLD_DELTA':int(sat['combined_mixed_old_delta_family_n'])
    }
    for k,v in expected_class_counts.items():
        if counts[k]!=v: raise RuntimeError(f'family class count changed for {k}: {counts[k]} != {v}')

    # sequence SHA -> family mapping from combined family membership.
    sha_to_family={}
    fam_delta_exact=Counter()
    with open(memp,newline='') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            s=r['sequence_sha256']; fr=r['family_rep']; oc=r['origin_class']
            if s in sha_to_family and sha_to_family[s]!=fr:
                raise RuntimeError('same exact AA SHA assigned to multiple combined families: '+s)
            sha_to_family[s]=fr
            if oc in ('DELTA_ONLY','SHARED_EXACT'):
                fam_delta_exact[fr]+=1

    # Global exact map is used to verify which SHA values are represented in delta.
    delta_sha=set(); delta_only_sha=set(); shared_sha=set()
    with open(map_p,newline='') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            c=r['origin_class']; s=r['sequence_sha256']
            if c in ('DELTA_ONLY','SHARED_EXACT'): delta_sha.add(s)
            if c=='DELTA_ONLY': delta_only_sha.add(s)
            elif c=='SHARED_EXACT': shared_sha.add(s)
    if len(delta_sha)!=int(sat['delta_exact_unique_n']):
        raise RuntimeError(f'delta exact SHA universe changed: {len(delta_sha)} != {sat["delta_exact_unique_n"]}')

    mans=read_tsv(manp)
    if len(mans)!=int(sat['delta_accession_n']) if 'delta_accession_n' in sat else 254:
        pass
    if len(mans)!=254:
        raise RuntimeError(f'expected 254 delta accessions, got {len(mans)}')
    if 'accession' not in mans[0]: raise RuntimeError('manifest lacks accession column')
    accs=[r['accession'] for r in mans]
    if len(set(accs))!=len(accs): raise RuntimeError('duplicate accession in delta manifest')
    stratum_cols=[c for c in mans[0] if 'stratum' in c.lower()]
    stratum_col=stratum_cols[0] if len(stratum_cols)==1 else None
    acc_to_stratum={r['accession']:(r.get(stratum_col,'') if stratum_col else '') for r in mans}

    fam_instances=Counter()
    fam_accessions=defaultdict(set)
    fam_strata=defaultdict(set)
    acc_family_unique=defaultdict(set)
    observed_delta_sha=set()
    total_records=0
    per_acc_rows=[]

    for acc in accs:
        fp=faa_dir/f'{acc}.formal.faa'
        if not fp.is_file(): raise RuntimeError('missing per-accession formal FASTA: '+str(fp))
        acc_seen_sha=set(); acc_inst=0
        for h,seq in fasta_iter(fp):
            if not seq: raise RuntimeError(f'empty sequence in {fp}: {h}')
            s=shaseq(seq); total_records+=1; acc_inst+=1
            if s not in delta_sha:
                raise RuntimeError(f'formal sequence SHA absent from Stage04A delta exact universe: accession={acc} sha={s}')
            fr=sha_to_family.get(s)
            if fr is None: raise RuntimeError('delta SHA missing combined family assignment: '+s)
            fam_instances[fr]+=1
            fam_accessions[fr].add(acc)
            if stratum_col and acc_to_stratum[acc]: fam_strata[fr].add(acc_to_stratum[acc])
            acc_family_unique[acc].add(fr)
            acc_seen_sha.add(s); observed_delta_sha.add(s)
        per_acc_rows.append({
          'accession':acc,
          'stratum':acc_to_stratum[acc] if stratum_col else '',
          'formal_instance_n':acc_inst,
          'exact_unique_n':len(acc_seen_sha),
          'combined_family_n':len(acc_family_unique[acc])
        })

    if total_records!=int(sat['delta_formal_rvt_n']):
        raise RuntimeError(f'formal instance count mismatch: {total_records} != {sat["delta_formal_rvt_n"]}')
    if observed_delta_sha!=delta_sha:
        miss=len(delta_sha-observed_delta_sha); extra=len(observed_delta_sha-delta_sha)
        raise RuntimeError(f'per-accession formal FASTAs do not reconstruct delta exact universe: missing={miss} extra={extra}')

    family_rows=[]
    for fr in sorted(famclass):
        cls=famclass[fr]
        if cls=='OLD_ONLY': continue
        accset=sorted(fam_accessions.get(fr,set()))
        strset=sorted(fam_strata.get(fr,set()))
        family_rows.append({
          'family_rep':fr,
          'family_class':cls,
          'n_delta_exact_sequences':fam_delta_exact.get(fr,0),
          'n_formal_instances':fam_instances.get(fr,0),
          'n_accessions':len(accset),
          'n_source_strata':len(strset) if stratum_col else '',
          'accessions':';'.join(accset),
          'source_strata':';'.join(strset) if stratum_col else ''
        })

    # Every delta-represented combined family must be recovered from the per-accession formal FASTAs.
    delta_fams={fr for fr,c in famclass.items() if c in ('DELTA_ONLY','MIXED_OLD_DELTA')}
    obs_fams={r['family_rep'] for r in family_rows if int(r['n_formal_instances'])>0}
    if obs_fams!=delta_fams:
        raise RuntimeError(f'delta family coverage mismatch missing={len(delta_fams-obs_fams)} extra={len(obs_fams-delta_fams)}')

    family_rows.sort(key=lambda r:(0 if r['family_class']=='DELTA_ONLY' else 1,-int(r['n_accessions']),-int(r['n_delta_exact_sequences']),-int(r['n_formal_instances']),r['family_rep']))
    fields=['family_rep','family_class','n_delta_exact_sequences','n_formal_instances','n_accessions','n_source_strata','accessions','source_strata']
    write_tsv(O/'DELTA_REPRESENTED_FAMILY_REPLICATION_CENSUS.tsv',family_rows,fields)
    write_tsv(O/'ACCESSION_FAMILY_BURDEN.tsv',per_acc_rows,['accession','stratum','formal_instance_n','exact_unique_n','combined_family_n'])

    delta_only=[r for r in family_rows if r['family_class']=='DELTA_ONLY']
    mixed=[r for r in family_rows if r['family_class']=='MIXED_OLD_DELTA']
    top=sorted(delta_only,key=lambda r:(-int(r['n_accessions']),-int(r['n_delta_exact_sequences']),-int(r['n_formal_instances']),r['family_rep']))[:1000]
    write_tsv(O/'TOP_REPLICATED_DELTA_ONLY_FAMILIES.tsv',top,fields)

    def n_ge(rows,k): return sum(int(r['n_accessions'])>=k for r in rows)
    def hist(rows):
        c=Counter(int(r['n_accessions']) for r in rows)
        return {str(k):c[k] for k in sorted(c)}
    def strata_ge(rows,k):
        if not stratum_col: return None
        return sum(int(r['n_source_strata'])>=k for r in rows)

    summary={
      'stage':'Stage05F3V2-Stage04B1-delta-family-replication-census',
      'gate':'PASS_F3V2_STAGE04B1_DELTA_FAMILY_REPLICATION_CENSUS_COMPLETE',
      'gate_pass':True,
      'delta_accession_n':len(accs),
      'delta_formal_instance_n':total_records,
      'delta_exact_unique_n':len(delta_sha),
      'delta_only_exact_n':len(delta_only_sha),
      'shared_exact_n':len(shared_sha),
      'combined_delta_only_family_n':len(delta_only),
      'combined_mixed_old_delta_family_n':len(mixed),
      'delta_only_family_accession_ge2_n':n_ge(delta_only,2),
      'delta_only_family_accession_ge3_n':n_ge(delta_only,3),
      'delta_only_family_accession_ge5_n':n_ge(delta_only,5),
      'delta_only_family_accession_ge10_n':n_ge(delta_only,10),
      'mixed_family_accession_ge2_n':n_ge(mixed,2),
      'mixed_family_accession_ge3_n':n_ge(mixed,3),
      'delta_only_family_accession_histogram':hist(delta_only),
      'source_stratum_column':stratum_col,
      'delta_only_family_strata_ge2_n':strata_ge(delta_only,2),
      'delta_only_family_strata_ge3_n':strata_ge(delta_only,3),
      'stage04a_gate_sha256':shafile(satp),
      'combined_family_membership_sha256':shafile(memp),
      'delta_selected_manifest_sha256':shafile(manp),
      'truth_or_annotation_used_for_selection':False,
      'priority_selection_performed':False,
      'structure_prediction_performed':False,
      'retuning_authorized':False,
      'science_change':False,
      'claim_boundary':'Analysis-only replication census of frozen delta-represented combined families. Multi-accession breadth is a replication descriptor, not proof of biological independence or function. No priority threshold, structure, locus architecture, or novelty claim is applied here.'
    }
    with open(O/'STAGE04B1_DELTA_FAMILY_REPLICATION_CENSUS.json','w') as f:
        json.dump(summary,f,indent=2); f.write('\n')

    # deterministic hash manifest for principal outputs
    outs=[O/'DELTA_REPRESENTED_FAMILY_REPLICATION_CENSUS.tsv',O/'TOP_REPLICATED_DELTA_ONLY_FAMILIES.tsv',O/'ACCESSION_FAMILY_BURDEN.tsv',O/'STAGE04B1_DELTA_FAMILY_REPLICATION_CENSUS.json']
    with open(O/'SHA256SUMS.txt','w') as f:
        for p in outs: f.write(f'{shafile(p)}  {p.name}\n')

    print(json.dumps(summary,indent=2))
    print('PASS_STAGE05F3V2_STAGE04B1_DELTA_FAMILY_REPLICATION_CENSUS_COMPLETE')

if __name__=='__main__':
    try: main()
    except Exception as e:
        print('[HOLD]',e,file=sys.stderr)
        raise
