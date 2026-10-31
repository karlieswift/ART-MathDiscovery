#!/usr/bin/env python3
import csv, json, hashlib, importlib.util, sys
from pathlib import Path

STARTS={'ATG','GTG','TTG'}
STOPS={'TAA','TAG','TGA'}
PARTNER_SCAN_BP=10000
MIN_AA=30

def loadmod(name,path):
    spec=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''): h.update(b)
    return h.hexdigest()

def wrap(s,w=80): return '\n'.join(s[i:i+w] for i in range(0,len(s),w))+'\n'

def call_forward_orfs(seq,codon,orientation):
    out=[]; n=len(seq)
    for frame in (0,1,2):
        seg_start=frame; i=frame
        while i+3<=n:
            c=seq[i:i+3]
            if c in STOPS:
                st=None
                for j in range(seg_start,i,3):
                    if seq[j:j+3] in STARTS:
                        st=j; break
                if st is not None:
                    aa_len=(i-st)//3
                    if aa_len>=MIN_AA:
                        pep=''.join(codon.get(seq[k:k+3],'X') for k in range(st,i,3))
                        out.append({'orientation':orientation,'frame_zero_based':frame,'local_start0':st,'local_end0_including_stop':i+3,'aa_length':aa_len,'peptide':pep})
                seg_start=i+3
            i+=3
    return out

def map_reverse_orfs(seq,rc_orfs):
    n=len(seq); out=[]
    for r in rc_orfs:
        s=n-int(r['local_end0_including_stop']); e=n-int(r['local_start0'])
        z=dict(r); z['local_start0']=s; z['local_end0_including_stop']=e; out.append(z)
    return out

def main():
    if len(sys.argv)!=8:
        raise SystemExit('usage: partner_orf_scan.py LOCUS_CENSUS GENOMEDB GENOMEDB_INDEX OUTDIR LOCK_JSON VENDOR_PLACEMENT VENDOR_REPEAT')
    locp,gdb,gidx,outdir,lockp,vplace,vrepeat=sys.argv[1:]
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True)
    lock=json.load(open(lockp)); assert lock['gate_pass'] is True
    placement=loadmod('vendor_place',vplace); rpt=loadmod('vendor_repeat',vrepeat)
    rows=list(csv.DictReader(open(locp,encoding='utf-8'),delimiter='\t'))
    if len(rows)!=10598: raise SystemExit(f'FAIL locus census n={len(rows)}')
    unique=[r for r in rows if r['placement_status']=='UNIQUE_EXACT']
    keys={r['parent_contig_key'] for r in unique}; gi=placement.load_index(gidx,keys)

    all_orfs=[]; nearest=[]; fasta=[]
    with open(gdb,'rb',buffering=0) as gfh:
        for ix,r in enumerate(unique,1):
            seq=placement.fetch_db(gfh,*gi[r['parent_contig_key']])
            s=int(r['candidate_start0']); e=int(r['candidate_end0']); strand=r['strand']
            if strand=='+': dn=seq[e:min(len(seq),e+PARTNER_SCAN_BP+3)]
            else: dn=placement.revcomp(seq[max(0,s-PARTNER_SCAN_BP-3):s])
            rt_stop_present=(len(dn)>=3 and dn[:3] in STOPS)
            base=3 if rt_stop_present else 0
            scan=dn[base:base+PARTNER_SCAN_BP]
            fwd=call_forward_orfs(scan,placement.CODON,'SAME_STRAND')
            rcseq=placement.revcomp(scan)
            rev=map_reverse_orfs(scan,call_forward_orfs(rcseq,placement.CODON,'OPPOSITE_STRAND'))
            ors=fwd+rev
            ors.sort(key=lambda z:(int(z['local_start0']),0 if z['orientation']=='SAME_STRAND' else 1,-int(z['aa_length']),int(z['local_end0_including_stop'])))
            for j,z in enumerate(ors,1):
                oid=f"{r['candidate_id']}__ORF{j:03d}"
                q={
                  'orf_id':oid,'family_rep_id':r['family_rep_id'],'candidate_id':r['candidate_id'],'accession':r['accession'],
                  'parent_contig_key':r['parent_contig_key'],'rt_strand':strand,'rt_terminal_stop_present':rt_stop_present,
                  'partner_scan_available_bp':len(scan),'physical_rank':j,'orientation_relative_to_rt':z['orientation'],
                  'gap_after_rt_stop_bp':int(z['local_start0']),'start_from_rt_coding_end_bp':base+int(z['local_start0']),
                  'orf_end_from_rt_stop_bp':int(z['local_end0_including_stop']),'aa_length':z['aa_length'],
                  'local_start0_after_rt_stop':z['local_start0'],'local_end0_including_stop_after_rt_stop':z['local_end0_including_stop'],
                }
                all_orfs.append(q)
                if j==1:
                    q2=dict(q); q2['peptide_sha256']=hashlib.sha256(z['peptide'].encode()).hexdigest(); nearest.append(q2)
                    fasta.append(f'>{oid}\n{wrap(z["peptide"])}')
            if not ors:
                nearest.append({
                  'orf_id':'','family_rep_id':r['family_rep_id'],'candidate_id':r['candidate_id'],'accession':r['accession'],
                  'parent_contig_key':r['parent_contig_key'],'rt_strand':strand,'rt_terminal_stop_present':rt_stop_present,
                  'partner_scan_available_bp':len(scan),'physical_rank':'','orientation_relative_to_rt':'',
                  'gap_after_rt_stop_bp':'','start_from_rt_coding_end_bp':'','orf_end_from_rt_stop_bp':'','aa_length':'',
                  'local_start0_after_rt_stop':'','local_end0_including_stop_after_rt_stop':'','peptide_sha256':'',
                })
            if ix%500==0: print(f'[PARTNER_ORF] {ix}/{len(unique)} nearest_with_orf={sum(bool(x["orf_id"]) for x in nearest)}',flush=True)

    allp=out/'F3B4C2_ALL_DOWNSTREAM_ORF_CANDIDATES.tsv'
    acols=list(all_orfs[0]) if all_orfs else ['orf_id','family_rep_id','candidate_id','accession','parent_contig_key','rt_strand','rt_terminal_stop_present','partner_scan_available_bp','physical_rank','orientation_relative_to_rt','gap_after_rt_stop_bp','start_from_rt_coding_end_bp','orf_end_from_rt_stop_bp','aa_length','local_start0_after_rt_stop','local_end0_including_stop_after_rt_stop']
    with allp.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=acols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(all_orfs)
    nearp=out/'F3B4C2_NEAREST_DOWNSTREAM_ORF_PER_RT_LOCUS.tsv'; ncols=list(nearest[0])
    with nearp.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=ncols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(nearest)
    fa=out/'F3B4C2_NEAREST_DOWNSTREAM_ORF_PROTEINS.faa'; fa.write_text(''.join(fasta),encoding='ascii')
    prep={
      'stage':'Stage05F3B4C2-partner-orf-prep','gate':'PASS_F3B4C2_NEAREST_DOWNSTREAM_ORF_ENUMERATION_COMPLETE','gate_pass':True,
      'input_unique_exact_rt_locus_n':len(unique),'rt_locus_n_with_nearest_complete_orf':sum(bool(x['orf_id']) for x in nearest),
      'partner_scan_bp_after_rt_stop':PARTNER_SCAN_BP,'minimum_complete_orf_aa':MIN_AA,'allowed_start_codons':sorted(STARTS),'stop_codons':sorted(STOPS),
      'same_and_opposite_strands_enumerated':True,'nearest_rule':'smallest physical start after RT stop; same-strand preferred only on exact positional tie; then longer ORF',
      'partner_function_decision_performed':False,'partner_dedicatedness_decision_performed':False,
      'all_orfs_sha256':sha256_file(allp),'nearest_orf_table_sha256':sha256_file(nearp),'nearest_fasta_sha256':sha256_file(fa),
      'claim_boundary':'De novo downstream ORF adjacency enumeration only. Nearest ORF is a partner candidate, not a confirmed dedicated partner gene.'}
    pp=out/'STAGE05F3B4C2_PARTNER_ORF_PREP_GATE.json'; pp.write_text(json.dumps(prep,indent=2)+'\n')
    print(json.dumps(prep,indent=2)); print('PASS_STAGE05F3B4C2_PARTNER_ORF_ENUMERATION_COMPLETE')

if __name__=='__main__': main()
