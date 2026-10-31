#!/usr/bin/env python3
import csv, json, hashlib, statistics, sys
from pathlib import Path
from collections import defaultdict, Counter, deque


def boolv(x): return str(x).lower()=='true'

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''): h.update(b)
    return h.hexdigest()

def read_fasta(path):
    out={}; h=None; seq=[]
    with open(path,encoding='utf-8',errors='ignore') as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if h is not None: out[h]=''.join(seq).upper()
                h=line[1:].split()[0]; seq=[]
            else: seq.append(line)
        if h is not None: out[h]=''.join(seq).upper()
    return out

def kmers(s,k=31):
    if len(s)<k: return {s}
    return {s[i:i+k] for i in range(len(s)-k+1)}

def jacc(a,b):
    if not a and not b: return 1.0
    return len(a&b)/len(a|b) if (a|b) else 0.0

def lev(a,b):
    if len(a)>len(b): a,b=b,a
    prev=list(range(len(a)+1))
    for i,cb in enumerate(b,1):
        cur=[i]
        for j,ca in enumerate(a,1):
            cur.append(min(cur[-1]+1,prev[j]+1,prev[j-1]+(ca!=cb)))
        prev=cur
    return prev[-1]

def components(nodes,edges):
    adj={n:set() for n in nodes}
    for a,b in edges:
        adj[a].add(b); adj[b].add(a)
    comps=[]; seen=set()
    for n in sorted(nodes):
        if n in seen: continue
        q=[n]; seen.add(n); c=[]
        while q:
            x=q.pop(); c.append(x)
            for y in sorted(adj[x]):
                if y not in seen: seen.add(y); q.append(y)
        comps.append(sorted(c))
    return comps


def nw_identity(a,b):
    # Deterministic global alignment: match +1, mismatch -1, gap -1.
    # State stores (score, matches, aligned_columns); ties prefer more matches, then fewer columns.
    m=len(a); n=len(b)
    prev=[(-j,0,j) for j in range(n+1)]
    for i in range(1,m+1):
        cur=[(-i,0,i)]
        for j in range(1,n+1):
            opts=[]
            sc,ma,al=prev[j-1]; opts.append((sc+(1 if a[i-1]==b[j-1] else -1),ma+(a[i-1]==b[j-1]),al+1))
            sc,ma,al=prev[j]; opts.append((sc-1,ma,al+1))
            sc,ma,al=cur[j-1]; opts.append((sc-1,ma,al+1))
            best=max(opts,key=lambda z:(z[0],z[1],-z[2]))
            cur.append(best)
        prev=cur
    sc,ma,al=prev[n]
    return ma/al if al else 1.0

def median_or_blank(xs):
    return statistics.median(xs) if xs else ''

def cv(xs):
    if len(xs)<2: return 0.0 if xs else ''
    m=statistics.mean(xs)
    return statistics.pstdev(xs)/m if m else ''

def write_tsv(path,rows,cols=None):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    if cols is None: cols=list(rows[0]) if rows else []
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)

def main():
    if len(sys.argv)!=9:
        raise SystemExit('usage: analyze_abc.py CORE_MANIFEST CORE_CONTEXT_FASTA B4C_ARRAYS B4C_FAMILY_SUMMARY CORE_PARTNER_FASTA TARGET7 OUTDIR LOCK')
    core_p,ctx_fa,arrays_p,famsum_p,partner_fa,target7_p,outdir,lockp=sys.argv[1:]
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True)
    lock=json.load(open(lockp)); assert lock['gate_pass'] is True
    core=list(csv.DictReader(open(core_p,encoding='utf-8'),delimiter='\t'))
    if len(core)!=16: raise SystemExit(f'FAIL core n={len(core)}')
    ctx=read_fasta(ctx_fa); pseqs=read_fasta(partner_fa)
    # FASTA headers start with candidate/orf IDs.
    if len(ctx)!=16 or len(pseqs)!=16: raise SystemExit(f'FAIL fasta dims ctx={len(ctx)} partner={len(pseqs)}')
    fam_order=[r['family_rep_id'] for r in csv.DictReader(open(target7_p,encoding='utf-8'),delimiter='\t')]
    famset=set(fam_order)
    if len(famset)!=7: raise SystemExit('FAIL target family n')

    # A: exact and near-context dereplication within family.
    byfam=defaultdict(list)
    for r in core: byfam[r['family_rep_id']].append(r)
    pairrows=[]; locusA=[]; famA=[]
    for fam in fam_order:
        xs=byfam[fam]
        seq={r['candidate_id']:ctx[r['candidate_id']] for r in xs}
        ks={cid:kmers(s,31) for cid,s in seq.items()}
        edges=[]
        for i,a in enumerate(sorted(seq)):
            for b in sorted(seq)[i+1:]:
                la,lb=len(seq[a]),len(seq[b]); lr=min(la,lb)/max(la,lb)
                jac=jacc(ks[a],ks[b])
                exact=seq[a]==seq[b]
                near=(lr>=0.95 and jac>=0.95)
                pairrows.append({'family_rep_id':fam,'candidate_a':a,'candidate_b':b,'length_a':la,'length_b':lb,'length_ratio':lr,'k31_jaccard':jac,'exact_architecture_core':exact,'near_context_edge':near})
                if near: edges.append((a,b))
        comps=components(seq.keys(),edges)
        compid={cid:f'{fam}__NC{j:02d}' for j,c in enumerate(comps,1) for cid in c}
        hash_groups=defaultdict(list)
        for cid,s in seq.items(): hash_groups[hashlib.sha256(s.encode()).hexdigest()].append(cid)
        exactid={cid:f'{fam}__EX{j:02d}' for j,(h,cids) in enumerate(sorted(hash_groups.items()),1) for cid in cids}
        for r in xs:
            z=dict(r); z['exact_context_cluster_id']=exactid[r['candidate_id']]; z['near_context_cluster_id']=compid[r['candidate_id']]; locusA.append(z)
        famA.append({
          'family_rep_id':fam,'core_architecture_locus_n':len(xs),'accession_n':len({r['accession'] for r in xs}),
          'exact_architecture_context_cluster_n':len(hash_groups),'near_context_cluster_n':len(comps),
          'near_context_dereplicated_accession_lower_bound':len(comps),
          'all_core_contexts_exactly_identical':len(hash_groups)==1,
          'all_core_contexts_near_identical_by_frozen_rule':len(comps)==1,
          'lineage_independence_confirmed':False,
        })
    write_tsv(out/'F3B4D_A_CORE_LOCUS_LINEAGE_DEREP_16.tsv',locusA)
    write_tsv(out/'F3B4D_A_PAIRWISE_CONTEXT_SIMILARITY.tsv',pairrows)
    write_tsv(out/'F3B4D_A_FAMILY_CONTEXT_DEREP_SUMMARY_7.tsv',famA)
    gateA={
      'stage':'Stage05F3B4D-A','gate':'PASS_F3B4D_A_NEAR_CONTEXT_DEREPLICATION_COMPLETE','gate_pass':True,
      'target_family_n':7,'core_architecture_locus_n':16,
      'exact_context_rule':'SHA256 of candidate-oriented architecture core from leftmost strict repeat start through recurrent nearest-partner ORF stop',
      'near_context_rule':'within-family 31-mer Jaccard>=0.95 AND architecture-core length ratio>=0.95; connected components used only as near-context dereplication clusters',
      'evolutionary_lineage_independence_confirmed':False,
      'family_summary_sha256':sha256_file(out/'F3B4D_A_FAMILY_CONTEXT_DEREP_SUMMARY_7.tsv'),
      'claim_boundary':'Near-context dereplication reduces obvious duplicate/near-duplicate SRA loci. It does not prove independent evolutionary lineages.'}
    (out/'STAGE05F3B4D_A_GATE.json').write_text(json.dumps(gateA,indent=2)+'\n')

    # B: strict-array repeat-core conservation across all strict loci in target 7 families.
    arrays=[r for r in csv.DictReader(open(arrays_p,encoding='utf-8'),delimiter='\t') if r['family_rep_id'] in famset and boolv(r['strict_published_art_geometry'])]
    if not arrays: raise SystemExit('FAIL no strict arrays for target families')
    motifrows=[]; famB=[]
    for fam in fam_order:
        ar=[r for r in arrays if r['family_rep_id']==fam]
        motifs=sorted(set(r['motif'] for r in ar))
        medges=[]
        for i,a in enumerate(motifs):
            for b in motifs[i+1:]:
                d=lev(a,b)
                if d<=1: medges.append((a,b))
        mcomps=components(motifs,medges)
        mid={m:f'{fam}__MC{j:02d}' for j,c in enumerate(mcomps,1) for m in c}
        # rank motif neighborhoods by independent accessions, array count, copy sum.
        stats=[]
        for j,c in enumerate(mcomps,1):
            xs=[r for r in ar if r['motif'] in c]
            stats.append((len({r['accession'] for r in xs}),len(xs),sum(int(r['copy_n']) for r in xs),sorted(c)[0],j,c,xs))
        stats.sort(key=lambda z:(-z[0],-z[1],-z[2],z[3]))
        dom=stats[0]
        allsp=[]; allcopies=[]; arrdist=[]; patterns=Counter(); motif_count=Counter()
        for r in ar:
            sps=[int(x) for x in r['spacer_lengths'].split(',') if x!='']
            sts=[int(x) for x in r['candidate_relative_repeat_starts'].split(',') if x!='']
            allsp += sps; allcopies.append(int(r['copy_n'])); motif_count[r['motif']]+=1
            if sts: arrdist.append(-max(sts))
            patterns[(r['motif'],r['copy_n'],r['spacer_lengths'])]+=1
            motifrows.append({
              'family_rep_id':fam,'candidate_id':r['candidate_id'],'accession':r['accession'],'motif':r['motif'],'motif_len':r['motif_len'],
              'motif_neighborhood_cluster_id':mid[r['motif']],'copy_n':r['copy_n'],'spacer_lengths':r['spacer_lengths'],'median_spacer_bp':r['median_spacer_bp'],
              'candidate_relative_repeat_starts':r['candidate_relative_repeat_starts'],'array_span_bp':r['array_span_bp']})
        modal_motif=sorted(motif_count.items(),key=lambda x:(-x[1],x[0]))[0][0]
        modal_pattern=sorted(patterns.items(),key=lambda x:(-x[1],x[0]))[0]
        famB.append({
          'family_rep_id':fam,'strict_array_n':len(ar),'strict_array_locus_n':len({r['candidate_id'] for r in ar}),'strict_array_accession_n':len({r['accession'] for r in ar}),
          'exact_motif_core_n':len(motifs),'one_edit_motif_neighborhood_n':len(mcomps),'dominant_motif_neighborhood_id':f'{fam}__MC{dom[4]:02d}',
          'dominant_motif_neighborhood_accession_n':dom[0],'dominant_motif_neighborhood_array_n':dom[1],'dominant_motif_neighborhood_members':','.join(sorted(dom[5])),
          'modal_exact_motif_core':modal_motif,'modal_exact_motif_array_n':motif_count[modal_motif],
          'copy_n_min':min(allcopies),'copy_n_median':statistics.median(allcopies),'copy_n_max':max(allcopies),
          'spacer_bp_median':median_or_blank(allsp),'spacer_bp_cv':cv(allsp),'spacer_bp_min':min(allsp) if allsp else '','spacer_bp_max':max(allsp) if allsp else '',
          'array_to_rt_nearest_repeat_distance_bp_median':median_or_blank(arrdist),
          'modal_exact_array_signature':f'{modal_pattern[0][0]}|copy={modal_pattern[0][1]}|spacers={modal_pattern[0][2]}','modal_exact_array_signature_n':modal_pattern[1],
          'full_biological_repeat_unit_consensus_inferred':False,
        })
    write_tsv(out/'F3B4D_B_STRICT_ARRAY_RECORDS_7FAM.tsv',motifrows)
    write_tsv(out/'F3B4D_B_REPEAT_CORE_CONSERVATION_SUMMARY_7.tsv',famB)
    gateB={
      'stage':'Stage05F3B4D-B','gate':'PASS_F3B4D_B_REPEAT_CORE_CONSERVATION_COMPLETE','gate_pass':True,
      'target_family_n':7,'strict_array_n':len(arrays),'strict_array_locus_n':len({r['candidate_id'] for r in arrays}),
      'motif_core_neighborhood_rule':'exact detector motif cores linked when Levenshtein distance<=1; no repeat-detector threshold retuning',
      'full_repeat_unit_consensus_claimed':False,'repeat_detector_changed':False,
      'family_summary_sha256':sha256_file(out/'F3B4D_B_REPEAT_CORE_CONSERVATION_SUMMARY_7.tsv'),
      'claim_boundary':'Conservation of detector-defined exact motif cores, copy number and spacer geometry only; not a full RNA repeat-unit sequence model.'}
    (out/'STAGE05F3B4D_B_GATE.json').write_text(json.dumps(gateB,indent=2)+'\n')

    # C base partner conservation for the 16 recurrent-partner strict loci.
    c_byfam=defaultdict(list)
    for r in core: c_byfam[r['family_rep_id']].append(r)
    partner_rows=[]; famC=[]
    for fam in fam_order:
        xs=c_byfam[fam]
        seqs=[]
        for r in xs:
            oid=r['partner_orf_id']; s=pseqs[oid]; seqs.append((r,s))
            partner_rows.append({'family_rep_id':fam,'candidate_id':r['candidate_id'],'accession':r['accession'],'partner_orf_id':oid,'partner_aa_length':len(s),'partner_peptide_sha256':hashlib.sha256(s.encode()).hexdigest(),'orientation':r['partner_orientation'],'gap_bp':r['partner_gap_bp']})
        hashes=Counter(hashlib.sha256(s.encode()).hexdigest() for _,s in seqs)
        lengths=[len(s) for _,s in seqs]; gaps=[int(r['partner_gap_bp']) for r,_ in seqs]
        pis=[]
        for i in range(len(seqs)):
            for j in range(i+1,len(seqs)):
                pis.append(nw_identity(seqs[i][1],seqs[j][1]))
        famC.append({
          'family_rep_id':fam,'core_partner_locus_n':len(xs),'partner_accession_n':len({r['accession'] for r,_ in seqs}),
          'exact_partner_aa_sequence_n':len(hashes),'most_common_exact_partner_aa_locus_n':max(hashes.values()),
          'partner_pairwise_global_identity_min':min(pis) if pis else 1.0,'partner_pairwise_global_identity_median':statistics.median(pis) if pis else 1.0,
          'partner_aa_length_min':min(lengths),'partner_aa_length_median':statistics.median(lengths),'partner_aa_length_max':max(lengths),
          'same_strand_n':sum(r['partner_orientation']=='SAME_STRAND' for r,_ in seqs),'opposite_strand_n':sum(r['partner_orientation']=='OPPOSITE_STRAND' for r,_ in seqs),
          'rt_partner_gap_bp_min':min(gaps),'rt_partner_gap_bp_median':statistics.median(gaps),'rt_partner_gap_bp_max':max(gaps),
          'partner_function_confirmed':False,'dedicated_partner_gene_confirmed':False,
        })
    write_tsv(out/'F3B4D_C_CORE_PARTNER_CONSERVATION_16.tsv',partner_rows)
    write_tsv(out/'F3B4D_C_PARTNER_CONSERVATION_SUMMARY_7.tsv',famC)
    gateC={
      'stage':'Stage05F3B4D-C','gate':'PASS_F3B4D_C_RECURRENT_PARTNER_CONSERVATION_BASE_COMPLETE','gate_pass':True,
      'target_family_n':7,'core_partner_locus_n':16,
      'conservation_dimensions':['exact AA sequence recurrence','AA length','orientation relative to RT','RT-to-partner gap'],
      'catalog_annotation_integrated':False,'partner_function_confirmed':False,'dedicated_partner_gene_confirmed':False,
      'family_summary_sha256':sha256_file(out/'F3B4D_C_PARTNER_CONSERVATION_SUMMARY_7.tsv'),
      'claim_boundary':'Recurrent nearest-downstream partner-candidate sequence/topology conservation only; function and dedicatedness require independent evidence.'}
    (out/'STAGE05F3B4D_C_BASE_GATE.json').write_text(json.dumps(gateC,indent=2)+'\n')

    print(json.dumps(gateA,indent=2)); print(json.dumps(gateB,indent=2)); print(json.dumps(gateC,indent=2))
    print('\n=== A FAMILY CONTEXT DEREP 7 ===')
    with open(out/'F3B4D_A_FAMILY_CONTEXT_DEREP_SUMMARY_7.tsv') as f: print(f.read().rstrip())
    print('\n=== B REPEAT CORE SUMMARY 7 ===')
    with open(out/'F3B4D_B_REPEAT_CORE_CONSERVATION_SUMMARY_7.tsv') as f: print(f.read().rstrip())
    print('\n=== C PARTNER CONSERVATION SUMMARY 7 ===')
    with open(out/'F3B4D_C_PARTNER_CONSERVATION_SUMMARY_7.tsv') as f: print(f.read().rstrip())
    print('PASS_F3B4D_ABC_BASE_ANALYSIS_COMPLETE')

if __name__=='__main__': main()
