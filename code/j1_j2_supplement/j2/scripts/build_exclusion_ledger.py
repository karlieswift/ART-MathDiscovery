#!/usr/bin/env python3
"""Index known exposure metadata without converting absence into independence.

Input is the frozen scientific release. Every entry retains a source file/column.
Protein hashes and local family IDs are kept distinct from whole-DNA hashes.
"""
import argparse,csv,hashlib,json,re,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent))
from j2_pipeline import fasta,digest,write_tsv,write_json,sha

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--release-root',required=True);ap.add_argument('--output',required=True);a=ap.parse_args()
 root=Path(a.release_root).resolve();out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
 ledger={};sources=[]
 def add(kind,value,path,column,cohort):
  value=str(value or '').strip()
  if not value or value in {'.','NA','None','nan','UNRESOLVED'}:return
  key=(kind,value)
  if key not in ledger:ledger[key]={'key_type':kind,'value':value,'first_source':str(path.relative_to(root)),
    'source_column':column,'exposed_cohort':cohort,'evidence_count':0}
  ledger[key]['evidence_count']+=1
 # Explicit source lists cover the formal panel, broad919, UG27 and both
 # discovery cohorts. Additional supplementary accession columns are indexed.
 paths=[root/'code/formal_baseline/config/CASES.tsv']
 paths+=list((root/'supplementary_data').glob('*.tsv'))
 paths+=list((root/'data/raw_supplementary_sources/closure_recovery').rglob('*MANIFEST.tsv'))
 paths+=list((root/'data/raw_supplementary_sources/closure_recovery').rglob('*ACCESSIONS.tsv'))
 paths+=list((root/'data/raw_supplementary_sources/closure_recovery').rglob('*ACQUISITION_LEDGER.tsv'))
 csv.field_size_limit(100_000_000)
 for p in sorted(set(paths)):
  if not p.is_file():continue
  rel=str(p.relative_to(root));n=0;cohort='PREVIOUSLY_EXPOSED_PROJECT_DATA'
  if 'CASES' in rel:cohort='FORMAL_13_INCLUDING_8_PRIMARY_3_CONTROLS_DEVELOPMENT_RETROSPECTIVE'
  elif 'NEW919' in rel or 'broad_genome' in rel:cohort='EXPOSED_BROAD_919'
  elif 'UG27' in rel:cohort='EXPOSED_UG27'
  elif 'independent' in rel or 'F3V2' in rel:cohort='EXPOSED_54G'
  elif 'discovery' in rel:cohort='EXPOSED_HISTORICAL_6_12G'
  with p.open(newline='') as f:
   for r in csv.DictReader(f,delimiter='\t'):
    n+=1
    for k,v in r.items():
     kl=k.lower()
     if kl in {'accession','reference_accession','sequence_accession','source_accession','priority_candidate_accession','native_parent_accession'}:
      add('accession',v,p,k,cohort)
     elif kl in {'native_contig_id','source_contig_id','scoring_id','assembly_id','case_id'}:
      add('contig_id',v,p,k,cohort)
     elif kl in {'family_rep_id','combined_novelty_family_id','delta_operational_family_id','family_id'}:
      add('family_cluster',v,p,k,cohort)
     elif kl in {'dna_sequence_sha256','reverse_complement_sha256'}:
      add('sequence_sha256' if kl.startswith('dna') else kl,v,p,k,cohort)
     elif kl=='sequence_sha256' and ('GENOME_MANIFEST' in p.name or 'broad_genome' in p.name):
      add('sequence_sha256',v,p,k,cohort)
     elif kl in {'aa_sequence_sha256','peptide_sha256'}:
      add('protein_sha256',v,p,k,cohort)
  sources.append({'path':rel,'sha256':sha(p),'rows':n})
 for p in sorted((root/'data/raw_supplementary_sources/benchmark_native_genomes').glob('*.fasta')):
  h,s=fasta(p);add('sequence_sha256',digest(s),p,'actual complete FASTA sequence','FORMAL_HISTORICAL_NATIVE_GENOME')
  add('reverse_complement_sha256',digest(s.translate(str.maketrans('ACGT','TGCA'))[::-1]),p,'actual complete FASTA reverse complement','FORMAL_HISTORICAL_NATIVE_GENOME')
  for token in re.findall(r'(?:[A-Z]{1,4}_?\d{5,9})(?:\.\d+)?',h+' '+p.name):add('accession',token,p,'FASTA header or filename','FORMAL_HISTORICAL_NATIVE_GENOME')
  sources.append({'path':str(p.relative_to(root)),'sha256':sha(p),'rows':1})
 entries=[ledger[k] for k in sorted(ledger)]
 write_tsv(out/'EXCLUSION_LEDGER.tsv',entries,['key_type','value','first_source','source_column','exposed_cohort','evidence_count'])
 write_tsv(out/'EXCLUSION_LEDGER_SOURCES.tsv',sources,['path','sha256','rows'])
 counts={}
 for x in entries:counts[x['key_type']]=counts.get(x['key_type'],0)+1
 report={'status':'KNOWN_EXPOSURE_INDEXED_NOT_A_COMPLETENESS_CERTIFICATE','entry_n':len(entries),'counts':counts,
  'source_file_n':len(sources),'ledger_sha256':sha(out/'EXCLUSION_LEDGER.tsv'),
  'limitations':['not evidence that an unlisted genome was never seen','local operational family IDs do not establish homology between new and old families',
   'not a sequence-similarity search','not all original native discovery contigs are supplied','protein hashes are not DNA hashes',
   'version-stripped accession checks are applied at admission; unresolved aliases require source curator review']}
 write_json(out/'EXCLUSION_LEDGER_STATUS.json',report);print(json.dumps(report,indent=2))
if __name__=='__main__':main()
