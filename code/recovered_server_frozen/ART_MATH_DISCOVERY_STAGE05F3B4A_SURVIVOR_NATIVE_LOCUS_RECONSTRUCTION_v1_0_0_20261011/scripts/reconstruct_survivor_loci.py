#!/usr/bin/env python3
import csv, hashlib, json, re, sys
from pathlib import Path

CODON = {
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}
RCMAP = str.maketrans('ACGTRYMKSWBDHVNacgtrymkswbdhvn','TGCAYRKMSWVHDBNtgcayrkmswvhdbn')


def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''): h.update(b)
    return h.hexdigest()

def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()

def clean_entry(b):
    return b.rstrip(b'\x00\r\n')

def revcomp(s):
    return s.translate(RCMAP)[::-1]

def translate_frame(seq, frame):
    out=[]
    n=len(seq)
    for i in range(frame,n-2,3):
        c=seq[i:i+3]
        out.append(CODON.get(c,'X'))
    return ''.join(out)

def read_fasta(path):
    d={}
    name=None; buf=[]
    with open(path,encoding='utf-8') as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if name is not None:
                    d[name]=''.join(buf).replace(' ','').upper()
                name=line[1:].split()[0]
                buf=[]
            else:
                buf.append(line)
        if name is not None:
            d[name]=''.join(buf).replace(' ','').upper()
    return d

def load_index_subset(index_path,wanted):
    out={}
    with open(index_path,encoding='ascii') as f:
        for line in f:
            if not line.strip(): continue
            x=line.split()
            if len(x)<3: raise SystemExit(f'FAIL malformed MMseqs index: {line[:200]!r}')
            k=x[0]
            if k in wanted:
                if k in out: raise SystemExit(f'FAIL duplicate DB key {k}')
                out[k]=(int(x[1]),int(x[2]))
    miss=wanted-set(out)
    if miss: raise SystemExit(f'FAIL keys absent from genomeDB.index n={len(miss)} examples={sorted(miss)[:20]}')
    return out

def fetch_seq(db_fh,off,ln):
    db_fh.seek(off); b=db_fh.read(ln)
    if len(b)!=ln: raise SystemExit(f'FAIL short DB read off={off} len={len(b)} expected={ln}')
    s=clean_entry(b).decode('ascii',errors='strict')
    s=re.sub(r'\s+','',s).upper()
    return s

def exact_peptide_matches(seq,pep):
    hits=[]
    # Forward frames
    for f in (0,1,2):
        aa=translate_frame(seq,f)
        pos=aa.find(pep)
        while pos!=-1:
            s=f+pos*3; e=s+len(pep)*3
            hits.append(('+',f,s,e))
            pos=aa.find(pep,pos+1)
    # Reverse-complement frames
    rc=revcomp(seq)
    n=len(seq)
    for f in (0,1,2):
        aa=translate_frame(rc,f)
        pos=aa.find(pep)
        while pos!=-1:
            rcs=f+pos*3; rce=rcs+len(pep)*3
            s=n-rce; e=n-rcs
            hits.append(('-',f,s,e))
            pos=aa.find(pep,pos+1)
    return hits

def wrap(seq,w=80):
    return '\n'.join(seq[i:i+w] for i in range(0,len(seq),w))+'\n'

def rt_keyword(text):
    t=text.lower()
    pats=['reverse transcriptase','rna-directed dna polymerase','rna directed dna polymerase','group ii intron','group-ii intron','ltra']
    return any(x in t for x in pats)

def main():
    if len(sys.argv)!=11:
        raise SystemExit('usage: reconstruct_survivor_loci.py R4_MAP F3B3D_SUMMARY F3B3D_GATE SHORTLIST_FASTA B3B_AUDIT GENOMEDB GENOMEDB_INDEX OUTDIR LOCK_JSON FLANK_BP')
    r4_map,d_sum,d_gate,short_fa,b_audit,genomedb,genomeidx,outdir,lock_json,flank_s=sys.argv[1:]
    flank=int(flank_s)
    out=Path(outdir); (out/'contexts_genomic').mkdir(parents=True,exist_ok=True); (out/'contexts_oriented').mkdir(parents=True,exist_ok=True)
    lock=json.load(open(lock_json))
    assert lock['gate_pass'] is True and lock['survivor_expected_n']==16 and flank==lock['flank_bp_each_side']==10000
    dg=json.load(open(d_gate)); assert dg['gate_pass'] is True

    with open(d_sum,encoding='utf-8') as f: drows=list(csv.DictReader(f,delimiter='\t'))
    if len(drows)!=32: raise SystemExit(f'FAIL F3B3D summary row_n={len(drows)} != 32')
    survivors=[r for r in drows if r['final_best_tier'] in ('STRONG','SUGGESTIVE')]
    if len(survivors)!=16: raise SystemExit(f'FAIL frozen survivor n={len(survivors)} != 16')
    survivor_ids={r['candidate_id'] for r in survivors}

    pep=read_fasta(short_fa)
    if not survivor_ids.issubset(pep):
        raise SystemExit(f'FAIL shortlist FASTA missing survivors {sorted(survivor_ids-set(pep))}')

    # R4 namespace map
    rm={}
    with open(r4_map,encoding='utf-8') as f:
        rr=csv.DictReader(f,delimiter='\t')
        for r in rr:
            if r['instance_id'] in survivor_ids: rm[r['instance_id']]=r
    if set(rm)!=survivor_ids: raise SystemExit(f'FAIL R4 map missing n={len(survivor_ids-set(rm))}')

    # B3B known-sequence audit is descriptive metadata only.
    ba={}
    with open(b_audit,encoding='utf-8') as f:
        br=csv.DictReader(f,delimiter='\t')
        for r in br: ba[r['candidate_id']]=r
    if not survivor_ids.issubset(ba): raise SystemExit('FAIL B3B audit missing survivor rows')

    wanted={rm[q]['parent_contig_key'] for q in survivor_ids}
    idx=load_index_subset(genomeidx,wanted)
    seq_by_key={}
    with open(genomedb,'rb',buffering=0) as fh:
        for k,(off,ln) in idx.items(): seq_by_key[k]=fetch_seq(fh,off,ln)

    out_rows=[]
    all_gen=[]; all_ori=[]; all_pep=[]
    for srow in sorted(survivors,key=lambda z:int(z['shortlist_rank'])):
        q=srow['candidate_id']; r=rm[q]; p=pep[q].rstrip('*')
        key=r['parent_contig_key']; seq=seq_by_key[key]
        hits=exact_peptide_matches(seq,p)
        status='UNIQUE_EXACT' if len(hits)==1 else ('NO_EXACT' if len(hits)==0 else 'MULTI_EXACT')
        strand=frame=start0=end0=''
        ctx_s=ctx_e=''; gsha=osha=''; gpath=opath=''
        if len(hits)==1:
            strand,frame,start0,end0=hits[0]
            if end0-start0 != len(p)*3: raise SystemExit('internal coordinate length error')
            ctx_s=max(0,start0-flank); ctx_e=min(len(seq),end0+flank)
            gseq=seq[ctx_s:ctx_e]
            oseq=gseq if strand=='+' else revcomp(gseq)
            gpath=out/'contexts_genomic'/f'{q}.plusminus10kb.genomic.fna'
            opath=out/'contexts_oriented'/f'{q}.plusminus10kb.candidate_oriented.fna'
            gh=(f'>{q}|accession={r["rescued_accession"]}|parent_key={key}|contig_len={len(seq)}|candidate={start0+1}-{end0}|strand={strand}|context={ctx_s+1}-{ctx_e}|orientation=genomic\n')
            oh=(f'>{q}|accession={r["rescued_accession"]}|parent_key={key}|contig_len={len(seq)}|candidate={start0+1}-{end0}|strand={strand}|context={ctx_s+1}-{ctx_e}|orientation=candidate_5to3\n')
            gpath.write_text(gh+wrap(gseq),encoding='ascii')
            opath.write_text(oh+wrap(oseq),encoding='ascii')
            gsha=sha256_file(gpath); osha=sha256_file(opath)
            all_gen.append((gh,gseq)); all_ori.append((oh,oseq))
        top_hdr=ba[q].get('top_reference_original_headers','')
        row={
            'shortlist_rank':srow['shortlist_rank'],
            'family_rep_id':srow['family_rep_id'],
            'candidate_id':q,
            'aa_length':len(p),
            'f3b3d_final_best_tier':srow['final_best_tier'],
            'source_accession':r['rescued_accession'],
            'parent_contig_key':key,
            'parent_contig_header':r['parent_contig_header'],
            'parent_contig_length_bp':len(seq),
            'exact_peptide_match_n':len(hits),
            'placement_status':status,
            'candidate_strand':strand,
            'translation_frame_zero_based':frame,
            'candidate_start0':start0,
            'candidate_end0':end0,
            'candidate_start1':(start0+1 if start0!='' else ''),
            'candidate_end1_inclusive':end0,
            'context_start0':ctx_s,
            'context_end0':ctx_e,
            'context_start1':(ctx_s+1 if ctx_s!='' else ''),
            'context_end1_inclusive':ctx_e,
            'context_left_bp':(start0-ctx_s if start0!='' else ''),
            'context_right_bp':(ctx_e-end0 if end0!='' else ''),
            'candidate_peptide_sha256':sha256_bytes(p.encode()),
            'genomic_context_fasta':str(gpath) if gpath else '',
            'genomic_context_fasta_sha256':gsha,
            'oriented_context_fasta':str(opath) if opath else '',
            'oriented_context_fasta_sha256':osha,
            'b3b_sequence_reference_hit_n':ba[q].get('sequence_reference_hit_n',''),
            'b3b_top_reference_id':ba[q].get('top_reference_id',''),
            'b3b_top_reference_header':top_hdr,
            'b3b_top_reference_explicit_rt_keyword':rt_keyword(top_hdr),
            'annotation_used_for_candidate_generation':False,
        }
        out_rows.append(row)
        all_pep.append((f'>{q}|family={srow["family_rep_id"]}|tier={srow["final_best_tier"]}|accession={r["rescued_accession"]}\n',p))

    cols=list(out_rows[0])
    locus_tsv=out/'F3B4A_SURVIVOR_NATIVE_LOCUS_MAP.tsv'
    with locus_tsv.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(out_rows)

    (out/'F3B4A_SURVIVOR_PEPTIDES_16.faa').write_text(''.join(h+wrap(s) for h,s in all_pep),encoding='ascii')
    (out/'F3B4A_SURVIVOR_CONTEXTS_16_GENOMIC.fna').write_text(''.join(h+wrap(s) for h,s in all_gen),encoding='ascii')
    (out/'F3B4A_SURVIVOR_CONTEXTS_16_CANDIDATE_ORIENTED.fna').write_text(''.join(h+wrap(s) for h,s in all_ori),encoding='ascii')

    unique=sum(r['placement_status']=='UNIQUE_EXACT' for r in out_rows)
    nohit=sum(r['placement_status']=='NO_EXACT' for r in out_rows)
    multi=sum(r['placement_status']=='MULTI_EXACT' for r in out_rows)
    # No biological filter: all 16 are retained in the next-stage manifest, even if placement HOLD.
    manifest=out/'F3B4B_FROZEN_INPUT_MANIFEST_16.tsv'
    mcols=['shortlist_rank','family_rep_id','candidate_id','f3b3d_final_best_tier','source_accession','parent_contig_key','candidate_strand','candidate_start0','candidate_end0','oriented_context_fasta','oriented_context_fasta_sha256','placement_status']
    with manifest.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=mcols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows({k:r[k] for k in mcols} for r in out_rows)

    passed=(unique==16 and nohit==0 and multi==0)
    gate={
        'stage':'Stage05F3B4A',
        'gate':'PASS_F3B4A_16_SURVIVOR_NATIVE_LOCI_EXACTLY_RECONSTRUCTED' if passed else 'HOLD_F3B4A_NATIVE_LOCUS_RECONSTRUCTION_REQUIRES_COORDINATE_REVIEW',
        'gate_pass':passed,
        'f3b3d_survivor_expected_n':16,
        'survivor_row_n':len(out_rows),
        'unique_exact_placement_n':unique,
        'no_exact_placement_n':nohit,
        'multi_exact_placement_n':multi,
        'flank_bp_each_side':flank,
        'translation_table':11,
        'translation_replay':'six frame exact peptide replay on canonical parent contig; no annotation coordinates used',
        'annotation_used_for_candidate_generation':False,
        'repeat_or_partner_analysis_performed':False,
        'art_catalogue_membership_decision_performed':False,
        'novelty_decision_performed':False,
        'survivor_membership_changed':False,
        'scientific_threshold_changed':False,
        'locus_map_sha256':sha256_file(locus_tsv),
        'b4b_input_manifest_sha256':sha256_file(manifest),
        'claim_boundary':'Exact native-locus/source-context reconstruction only. This stage does not establish ART architecture, biochemical RT activity, or novelty.'
    }
    with open(out/'STAGE05F3B4A_GATE.json','w') as f: json.dump(gate,f,indent=2); f.write('\n')
    print(json.dumps(gate,indent=2))
    if not passed: raise SystemExit(3)

if __name__=='__main__': main()
