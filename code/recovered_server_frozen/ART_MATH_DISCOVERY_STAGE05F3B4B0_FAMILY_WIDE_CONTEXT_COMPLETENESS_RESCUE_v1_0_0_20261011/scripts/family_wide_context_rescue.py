#!/usr/bin/env python3
import csv, json, hashlib, re, sys
from pathlib import Path
from collections import defaultdict

CODON={
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'R','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
# Correct code-11 standard codon: AAG is K. Fix explicitly to avoid typo propagation.
CODON['AAG']='K'
RCMAP=str.maketrans('ACGTRYMKSWBDHVNacgtrymkswbdhvn','TGCAYRKMSWVHDBNtgcayrkmswvhdbn')
BINS=(500,2000,5000,10000)

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''): h.update(b)
    return h.hexdigest()

def clean_entry(b): return b.rstrip(b'\x00\r\n')
def revcomp(s): return s.translate(RCMAP)[::-1]
def translate(seq,frame): return ''.join(CODON.get(seq[i:i+3],'X') for i in range(frame,len(seq)-2,3))
def fetch_db(dbfh,off,ln):
    dbfh.seek(off); b=dbfh.read(ln)
    if len(b)!=ln: raise SystemExit(f'FAIL short DB read off={off} got={len(b)} expected={ln}')
    return re.sub(r'\s+','',clean_entry(b).decode('ascii','strict')).upper()
def load_index(path,wanted):
    out={}
    with open(path,encoding='ascii') as f:
        for line in f:
            if not line.strip(): continue
            x=line.split()
            if len(x)<3: raise SystemExit(f'FAIL malformed index {line[:120]!r}')
            if x[0] in wanted: out[x[0]]=(int(x[1]),int(x[2]))
    miss=wanted-set(out)
    if miss: raise SystemExit(f'FAIL missing DB keys n={len(miss)} examples={sorted(miss)[:10]}')
    return out

def exact_hits(seq,pep):
    hits=[]
    for f in (0,1,2):
        aa=translate(seq,f); pos=aa.find(pep)
        while pos>=0:
            s=f+3*pos; hits.append(('+',f,s,s+3*len(pep))); pos=aa.find(pep,pos+1)
    rc=revcomp(seq); n=len(seq)
    for f in (0,1,2):
        aa=translate(rc,f); pos=aa.find(pep)
        while pos>=0:
            rs=f+3*pos; re_=rs+3*len(pep); hits.append(('-',f,n-re_,n-rs)); pos=aa.find(pep,pos+1)
    return hits

def wrap(s,w=80): return '\n'.join(s[i:i+w] for i in range(0,len(s),w))+'\n'

def main():
    if len(sys.argv)!=10:
        raise SystemExit('usage: family_wide_context_rescue.py B4A_LOCUS_MAP R4_FAMILY_MAP R4_INSTANCE_MAP ORFDB ORFDB_INDEX GENOMEDB GENOMEDB_INDEX OUTDIR LOCK_JSON')
    b4a,fmap_p,imap_p,orfdb,orfidx,gdb,gidx,outdir,lockp=sys.argv[1:]
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True); (out/'best_contexts').mkdir(exist_ok=True)
    lock=json.load(open(lockp)); assert lock['gate_pass'] is True

    with open(b4a,encoding='utf-8') as f: b4=list(csv.DictReader(f,delimiter='\t'))
    fam_order=[]; tier={}; priority_candidate={}
    for r in sorted(b4,key=lambda x:int(x['shortlist_rank'])):
        fam=r['family_rep_id']
        if fam not in tier:
            fam_order.append(fam); tier[fam]=r['f3b3d_final_best_tier']; priority_candidate[fam]=r['candidate_id']
    if len(fam_order)!=16: raise SystemExit(f'FAIL survivor family n={len(fam_order)} !=16')
    famset=set(fam_order)

    members=[]
    with open(fmap_p,encoding='utf-8') as f:
        rr=csv.DictReader(f,delimiter='\t')
        req={'candidate_id','family_rep_id','accession'}
        if not req.issubset(rr.fieldnames or []): raise SystemExit(f'FAIL R4 family map missing {sorted(req-set(rr.fieldnames or []))}')
        for r in rr:
            if r['family_rep_id'] in famset: members.append(r)
    byfam=defaultdict(list)
    for r in members: byfam[r['family_rep_id']].append(r)
    missing=famset-set(byfam)
    if missing: raise SystemExit(f'FAIL no members for survivor families {sorted(missing)}')

    im={}
    with open(imap_p,encoding='utf-8') as f:
        rr=csv.DictReader(f,delimiter='\t')
        req={'instance_id','orfdb_key','parent_contig_key','rescued_accession','parent_contig_header'}
        if not req.issubset(rr.fieldnames or []): raise SystemExit(f'FAIL R4 instance map missing {sorted(req-set(rr.fieldnames or []))}')
        wanted={r['candidate_id'] for r in members}
        for r in rr:
            if r['instance_id'] in wanted: im[r['instance_id']]=r
    miss={r['candidate_id'] for r in members}-set(im)
    if miss: raise SystemExit(f'FAIL instance provenance missing n={len(miss)} examples={sorted(miss)[:10]}')

    orfkeys={im[r['candidate_id']]['orfdb_key'] for r in members}
    contigkeys={im[r['candidate_id']]['parent_contig_key'] for r in members}
    oi=load_index(orfidx,orfkeys); gi=load_index(gidx,contigkeys)
    peps={}; seqs={}
    with open(orfdb,'rb',buffering=0) as fh:
        for k,(o,l) in oi.items(): peps[k]=fetch_db(fh,o,l).rstrip('*')
    with open(gdb,'rb',buffering=0) as fh:
        for k,(o,l) in gi.items(): seqs[k]=fetch_db(fh,o,l)

    rows=[]
    for m in members:
        q=m['candidate_id']; p=im[q]; pep=peps[p['orfdb_key']]; seq=seqs[p['parent_contig_key']]
        hits=exact_hits(seq,pep)
        status='UNIQUE_EXACT' if len(hits)==1 else ('NO_EXACT' if not hits else 'MULTI_EXACT')
        strand=frame=start=end=''; left=right=minflank=totflank=''
        if len(hits)==1:
            strand,frame,start,end=hits[0]; left=start; right=len(seq)-end; minflank=min(left,right); totflank=left+right
        row={
            'family_rep_id':m['family_rep_id'],'survivor_tier':tier[m['family_rep_id']],
            'priority_candidate_id':priority_candidate[m['family_rep_id']], 'candidate_id':q,
            'accession':p['rescued_accession'],'parent_contig_key':p['parent_contig_key'],
            'parent_contig_header':p['parent_contig_header'],'parent_contig_length_bp':len(seq),
            'aa_length':len(pep),'exact_peptide_match_n':len(hits),'placement_status':status,
            'strand':strand,'frame_zero_based':frame,'candidate_start0':start,'candidate_end0':end,
            'left_available_bp':left,'right_available_bp':right,'min_flank_bp':minflank,'total_flank_bp':totflank,
        }
        for b in BINS: row[f'complete_plusminus_{b}bp']=bool(len(hits)==1 and minflank>=b)
        rows.append(row)

    cols=list(rows[0])
    allp=out/'F3B4B0_ALL_SURVIVOR_FAMILY_MEMBER_LOCUS_COMPLETENESS.tsv'
    with allp.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)

    # Technical best-context representative: unique placement, max min flank, then total flank, contig length, candidate_id.
    bestrows=[]; summaries=[]
    for fam in fam_order:
        fr=[r for r in rows if r['family_rep_id']==fam]
        uniq=[r for r in fr if r['placement_status']=='UNIQUE_EXACT']
        best=None
        if uniq:
            best=sorted(uniq,key=lambda r:(-int(r['min_flank_bp']),-int(r['total_flank_bp']),-int(r['parent_contig_length_bp']),r['candidate_id']))[0]
            seq=seqs[best['parent_contig_key']]; s=int(best['candidate_start0']); e=int(best['candidate_end0'])
            cs=max(0,s-10000); ce=min(len(seq),e+10000); ctx=seq[cs:ce]
            if best['strand']=='-': ctx=revcomp(ctx)
            fp=out/'best_contexts'/f'{fam}__{best["candidate_id"]}.best_observed_context.fna'
            fp.write_text(f'>{fam}|candidate={best["candidate_id"]}|accession={best["accession"]}|min_flank={best["min_flank_bp"]}|contig_len={best["parent_contig_length_bp"]}|orientation=candidate_5to3\n'+wrap(ctx),encoding='ascii')
            br=dict(best); br['best_context_fasta']=str(fp); br['best_context_fasta_sha256']=sha256_file(fp); bestrows.append(br)
        accn=len({r['accession'] for r in fr})
        sm={
            'family_rep_id':fam,'survivor_tier':tier[fam],'priority_candidate_id':priority_candidate[fam],
            'family_candidate_instance_n':len(fr),'family_accession_n':accn,
            'unique_exact_placement_n':sum(r['placement_status']=='UNIQUE_EXACT' for r in fr),
            'no_exact_placement_n':sum(r['placement_status']=='NO_EXACT' for r in fr),
            'multi_exact_placement_n':sum(r['placement_status']=='MULTI_EXACT' for r in fr),
            'max_min_flank_bp':max([int(r['min_flank_bp']) for r in uniq],default=0),
            'best_context_candidate_id':best['candidate_id'] if best else '',
            'best_context_accession':best['accession'] if best else '',
            'best_context_parent_contig_length_bp':best['parent_contig_length_bp'] if best else '',
        }
        for b in BINS:
            sm[f'instance_n_complete_plusminus_{b}bp']=sum(str(r[f'complete_plusminus_{b}bp']).lower()=='true' for r in fr)
            sm[f'accession_n_with_complete_plusminus_{b}bp']=len({r['accession'] for r in fr if str(r[f'complete_plusminus_{b}bp']).lower()=='true'})
        summaries.append(sm)

    bestp=out/'F3B4B0_TECHNICAL_BEST_CONTEXT_PER_SURVIVOR_FAMILY.tsv'
    bcols=list(bestrows[0]) if bestrows else cols+['best_context_fasta','best_context_fasta_sha256']
    with bestp.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=bcols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(bestrows)
    sump=out/'F3B4B0_SURVIVOR_FAMILY_CONTEXT_COMPLETENESS_SUMMARY.tsv'
    scols=list(summaries[0])
    with sump.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=scols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(summaries)

    gate={
      'stage':'Stage05F3B4B0','gate':'PASS_F3B4B0_FAMILY_WIDE_CONTEXT_COMPLETENESS_RESCUE_COMPLETE','gate_pass':True,
      'survivor_family_n':16,'selected_family_candidate_instance_n':len(rows),
      'selected_family_accession_n':len({r['accession'] for r in rows}),
      'unique_exact_placement_n':sum(r['placement_status']=='UNIQUE_EXACT' for r in rows),
      'no_exact_placement_n':sum(r['placement_status']=='NO_EXACT' for r in rows),
      'multi_exact_placement_n':sum(r['placement_status']=='MULTI_EXACT' for r in rows),
      'family_n_with_at_least_one_complete_plusminus_500bp':sum(s['instance_n_complete_plusminus_500bp']>0 for s in summaries),
      'family_n_with_at_least_one_complete_plusminus_2000bp':sum(s['instance_n_complete_plusminus_2000bp']>0 for s in summaries),
      'family_n_with_at_least_one_complete_plusminus_5000bp':sum(s['instance_n_complete_plusminus_5000bp']>0 for s in summaries),
      'family_n_with_at_least_one_complete_plusminus_10000bp':sum(s['instance_n_complete_plusminus_10000bp']>0 for s in summaries),
      'technical_best_context_rule':'unique exact placement only; maximize min(left,right), then total flank, then parent contig length, then candidate_id lexical; no repeat/annotation/novelty information used',
      'context_completeness_bins_bp_each_side':[500,2000,5000,10000],
      'repeat_analysis_performed':False,'partner_analysis_performed':False,'novelty_decision_performed':False,
      'architecture_absence_claim_authorized':False,
      'all_instances_table_sha256':sha256_file(allp),'family_summary_sha256':sha256_file(sump),'best_context_table_sha256':sha256_file(bestp),
      'claim_boundary':'Family-wide technical context rescue only. Completeness is a sequencing/assembly property; no ART/repeat/partner or novelty claim.'}
    gp=out/'STAGE05F3B4B0_GATE.json'; gp.write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate,indent=2))
    print('\n=== FAMILY CONTEXT COMPLETENESS 16 ===')
    print('\t'.join(scols))
    for s in summaries: print('\t'.join(str(s[c]) for c in scols))
    print('PASS_STAGE05F3B4B0_FAMILY_WIDE_CONTEXT_COMPLETENESS_RESCUE_COMPLETE')

if __name__=='__main__': main()
