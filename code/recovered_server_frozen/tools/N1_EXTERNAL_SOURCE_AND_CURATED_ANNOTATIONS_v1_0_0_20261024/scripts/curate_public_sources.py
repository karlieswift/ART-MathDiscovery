#!/usr/bin/env python3
"""Frozen, source-backed external ART N1 candidate registry; never assign new primary truth."""
import argparse
import csv
import hashlib
import json
import pathlib
import re

SOURCE_ROWS = [
 ('YOON_ART_2026','ART_PUBLISHED','Yoon et al. ART family','ART genomic RT+partner+repeat arrays reported; RT clusters are not independently verified array loci','PUBLISHED_CLASS_LEVEL','HISTORICALLY_EXPOSED','HOLD', 'Needs original per-locus full-genome array+RT coordinates and newly disjoint families','https://www.biorxiv.org/content/10.64898/2026.10.05.755809v1'),
 ('AI3_L0050','ART_REPRODUCTION','L0050 or closely related assembly','External group recovered published locus fingerprint; identity with original assembly not proven','PUBLISHED_SAMPLE_LEVEL','HISTORICALLY_EXPOSED','HOLD','Known development example and source-assisted reproduction','https://ai3-discovery.onpod.ai/research/art-loci-independent-reproduction'),
 ('AI3_ART34','ART_REPRODUCTION','ERR8055547_17666_12 / ART_34','Authors report alignment of gene sizes and repeat intervals to published ART_34; source is Logan contig ID, not INSDC accession.version','PUBLISHED_SAMPLE_LEVEL','HISTORICALLY_EXPOSED','HOLD','Missing independent array/RT native coordinates; matches illustrated ART','https://ai3-discovery.onpod.ai/research/art-loci-independent-reproduction'),
 ('AI3_ART69','ART_REPRODUCTION','SRR8925777_19762_9 / ART_69','Authors report alignment to published ART_69; source is Logan contig ID','PUBLISHED_SAMPLE_LEVEL','HISTORICALLY_EXPOSED','HOLD','Matches illustrated ART; source is not unique INSDC accession.version','https://ai3-discovery.onpod.ai/research/art-loci-independent-reproduction'),
 ('AI3_ERR10896470','ART_LIKE_AUTHOR_CANDIDATE','ERR10896470_6887_1 / lineage','External group reports short upstream context (232 nt) and read-extended upstream sequence (3597 nt); lineage strong array evidence in multiple samples','AUTHOR_REPORTED_LOCUS_ARCHITECTURE','PREVIOUSLY_REVIEWED_BY_PROJECT','HOLD','No native array/RT coordinates/consensus FASTA published in article; not certified absent from complete ART95 or developer exposure','https://ai3-discovery.onpod.ai/research/art-loci-independent-reproduction'),
 ('AI3_ERR4236129','ART_LINKED_CONTEXT_RECOVERY','ERR4236129_15347_m','External group reports upstream extension 301 nt to 3585 nt enabling array assessment; no per-locus final array call or native coordinates in article','AUTHOR_REPORTED_UPSTREAM_RECONSTRUCTION','PREVIOUSLY_REVIEWED_BY_PROJECT','HOLD','Context recovery alone is not ART true-positive annotation','https://ai3-discovery.onpod.ai/research/art-loci-independent-reproduction'),
 ('MINERVA_UG27','NON_ART_HARD_NEIGHBOR','UG27 RT-ncRNA arrays in prophages','Authors report ~150-nt structured ncRNA repeat units and cDNA synthesis requiring RT and accessory proteins; distinct system class','EXPERIMENTALLY_SUPPORTED_OTHER_CLASS','HISTORICALLY_EXPOSED','HOLD_NON_ART','Related negative-class challenge / positive for generic RT-ncRNA architecture, not ART-positive','https://www.biorxiv.org/content/10.64898/2026.09.22.753630v2.full'),
 ('MINERVA_ZENODO_TABLE','AUTHOR_SUPPLEMENT_PENDING_BYTES','Minerva_SupplementaryTables_preprint.xlsx','Zenodo confirms a 28.4 kB publicly listed supplementary workbook; byte content not inspected in this environment','PUBLIC_FILE_METADATA_ONLY','HISTORICALLY_EXPOSED','HOLD','Download and inspect spreadsheet to recover accession and coordinate content before use','https://zenodo.org/records/22903368'),
 ('INPHARED2_AUG2026','UNLABELED_EXTERNAL_UNIVERSE','INPHARED2 complete phage genomes','Curated ~37387 phage genomes with quality and taxonomy fields; no ART array/RT truth as part of displayed download metadata','GENOME_UNIVERSE_ONLY','EXPOSURE_TO_BE_AUDITED','HOLD_UNLABELED','New phage genomes are not verified negatives or positives; conduct accession+sequence+family exposure audit and blinded locus curation','https://phagebank.uk/inphared2/'),
 ('DRT3_2026','NON_ART_BIOLOGICAL_OUTGROUP','Deng et al. DRT3','Experimentally characterized dual RT / ncRNA producing poly(GT/AC) DNA; distinct defense-associated RT architecture','EXPERIMENTALLY_SUPPORTED_OTHER_CLASS','EXPOSURE_TO_BE_AUDITED','HOLD_NON_ART','Do not interpret DRT3 as array-associated ART without evidence','https://doi.org/10.1126/science.aed1656'),
 ('KORN_JVI_2021','PUBLISHED_GENOME_ANNOTATION','MarsHill / Madawaska / Machias','Korn et al. 2021 originally annotated MarsHill gp206 as a 501-residue retron-like RT and described roughly 1200bp upstream intergenic DNA; MW248466 genome 266637bp. Repeats were not independently established by that study','AUTHOR_GENOME_RECORD_SOURCE','HISTORICALLY_EXPOSED','HOLD','Already studied phage family, not a new blind test; no full array+RT coordinates in source registry','https://pmc.ncbi.nlm.nih.gov/articles/PMC8428398/'),
 ('MARSHILL_MW248466','KNOWN_ART_EXPOSED_NATIVE_RECHECK','MW248466.1 Staphylococcus phage MarsHill','Five native 18-nt anchor sequences and exact coordinates independently re-read from archived full genome; no complete ART array or RT interval newly curated here','NUCLEOTIDE_ANCHOR_INTERVAL_VERIFIED','HISTORICALLY_EXPOSED','HOLD','Use as real evidence-chain positive control for source-QC only, never unseen N1 benchmark',''),
]
FIELDS=['source_id','curation_category','sample_or_cohort','source_backed_fact','annotation_resolution','project_contact_status','n1_primary_gate','missing_evidence_or_reason','evidence_url']


def sha256(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for block in iter(lambda:f.read(1<<20),b''): h.update(block)
 return h.hexdigest()


def tsv_read(path):
 with open(path,newline='',encoding='utf-8') as f: return list(csv.DictReader(f,delimiter='\t'))


def tsv_write(path,fields,rows):
 with open(path,'w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)


def parse_fasta(path):
 lines=path.read_text(encoding='utf-8').splitlines();head=[l for l in lines if l.startswith('>')]
 if len(head)!=1:raise ValueError('Expected one native MarsHill sequence')
 if 'MW248466.1' not in head[0]:raise ValueError('Wrong native accession header')
 seq=''.join(l.strip().upper() for l in lines if l and not l.startswith('>'))
 if set(seq)-set('ACGTN'):raise ValueError('Invalid MarsHill DNA')
 return seq


def audit_mars_native(kit,output):
 fasta=kit/'data/MarsHill_MW248466.1_exposed_reference.fasta'
 intervals=tsv_read(kit/'data/MarsHill_Five_Real_18nt_Annotations.tsv')
 seq=parse_fasta(fasta)
 if len(seq)!=266637: raise AssertionError(f'Unexpected source genome length: {len(seq)}')
 expected_starts=[213994,214207,214404,214614,214787]
 positions=[];rows=[];errors=[]
 for row in intervals:
  a=int(row['start0']);b=int(row['end0']);k=row['native_18nt'].upper()
  observed=seq[a:b];positions.append(a)
  if len(k)!=18 or observed!=k or b-a!=18:errors.append((a,k,observed))
  rows.append(dict(row,independent_native_read=observed,native_match=observed==k,
                   genome_fasta_sha256=sha256(fasta),evidence_class='EXPOSED_TRUE_NUCLEOTIDE_ANCHOR_NOT_NEW_ART_TRUTH'))
 if positions!=expected_starts or errors:raise AssertionError(f'Anchor mismatch {errors}')
 majority=''.join(max('ACGT',key=lambda x:sum(s['native_18nt'][j]==x for s in rows)) for j in range(18))
 if majority!='ATATGAATACGTATAATG':raise AssertionError(f'Majority changed {majority}')
 tsv_write(output/'N1_MARSHILL_5_REAL_NATIVE_ANCHORS.tsv', list(rows[0].keys()), rows)
 return {'n_verified_anchor_copies':len(rows),'genome_length':len(seq),'genome_sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),
         'whole_fasta_sha256':sha256(fasta),'majority_18nt':majority,'starts0':positions,'spacing_bp':[b-a for a,b in zip(positions,positions[1:])],
         'primary_benchmark_admissible':False,'primary_hold_reason':'Previously studied MarsHill; only anchor copies checked, not full array/RT intervals'}


def main():
 p=argparse.ArgumentParser();p.add_argument('--kit',default=str(pathlib.Path(__file__).resolve().parents[1]));p.add_argument('--out',required=True)
 a=p.parse_args();kit=pathlib.Path(a.kit).resolve();out=pathlib.Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True)
 source=[dict(zip(FIELDS,s)) for s in SOURCE_ROWS]
 tsv_write(out/'N1_PUBLIC_SOURCE_CURATION.tsv',FIELDS,source)
 exposed=tsv_read(kit/'data/EXPOSURE_919_AND_HISTORICAL_ACCESSIONS.tsv')
 exposure_set={r['accession_version'] for r in exposed}
 mars=audit_mars_native(kit,out)
 if 'MW248466.1' not in exposure_set:
  # MarsHill exposure is established by project files and published prior handling,
  # regardless of membership in the 919 accession inventory.
  explicit_project_mars_historical=True
 else:explicit_project_mars_historical=False
 # Import template remains strictly empty: source-level facts are not genome-level ART coordinates.
 header=(kit/'data/N1_PRIMARY_EMPTY_INTAKE.tsv').read_text(encoding='utf-8').splitlines()[0]
 (out/'N1_NEW_PRIMARY_INTAKE_STILL_EMPTY.tsv').write_text(header+'\n',encoding='utf-8')
 facts=[]
 for row in source:
  facts.append({'source_id':row['source_id'],'sample_or_cohort':row['sample_or_cohort'],
                'actual_observation':row['source_backed_fact'],'source_url':row['evidence_url'],
                'verified_array_native_coords':'NO (MarsHill has only 5 true individual anchors)',
                'verified_rt_native_coords':'NO','confirmed_unseen':'NO','eligible_new_n1_primary':'NO',
                'evidence_scope':row['annotation_resolution']})
 tsv_write(out/'N1_REAL_SOURCE_FACTS_NOT_BLIND_TRUTH.tsv',list(facts[0].keys()),facts)
 report={'gate':'HOLD_NO_INDEPENDENT_ART_ARRAY_RT_NATIVE_TRUTH','n_public_sources_screened':len(source),
         'n_known_native_motif_copies_verified':mars['n_verified_anchor_copies'],
         'n_new_primary_art_loci_with_verified_array_and_rt_coords':0,'n_admitted_external_primary_genomes':0,
         'n_unique_exposure_accessions':len(exposure_set),'explicit_prior_mars_hill_exposure_outside_919_list':explicit_project_mars_historical,
         'mars_native_audit':mars,'critical_caveat':'Real published source facts and DNA anchor annotations are not full genome-native ART array+RT ground truth.'}
 (out/'N1_PUBLIC_CURATION_GATE.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
 print(json.dumps({k:v for k,v in report.items() if k!='mars_native_audit'},indent=2))

if __name__=='__main__': main()
