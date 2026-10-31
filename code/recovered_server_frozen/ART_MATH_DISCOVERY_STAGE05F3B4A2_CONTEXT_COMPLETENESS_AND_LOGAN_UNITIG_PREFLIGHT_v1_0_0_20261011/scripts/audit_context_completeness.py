#!/usr/bin/env python3
import argparse,csv,hashlib,json
from pathlib import Path

def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def read_one_fasta(path):
    name=None; seq=[]
    with open(path) as f:
        for ln in f:
            ln=ln.strip()
            if not ln: continue
            if ln.startswith('>'):
                if name is not None: raise RuntimeError(f'multi-record FASTA: {path}')
                name=ln[1:]
            else: seq.append(ln.upper())
    if name is None: raise RuntimeError(f'empty FASTA: {path}')
    return name,''.join(seq)

ap=argparse.ArgumentParser()
ap.add_argument('--b4a',required=True)
ap.add_argument('--out',required=True)
a=ap.parse_args()
b4a=Path(a.b4a); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
gate=json.load(open(b4a/'outputs/STAGE05F3B4A_GATE.json'))
assert gate['gate_pass'] and gate['survivor_row_n']==16 and gate['unique_exact_placement_n']==16, gate
mp=b4a/'outputs/F3B4A_SURVIVOR_NATIVE_LOCUS_MAP.tsv'
rows=list(csv.DictReader(open(mp),delimiter='\t'))
assert len(rows)==16
assert len({r['candidate_id'] for r in rows})==16
want=10000
outrows=[]; anchors=[]
for r in rows:
    clen=int(r['parent_contig_length_bp']); s=int(r['candidate_start0']); e=int(r['candidate_end0'])
    left=int(r['context_left_bp']); right=int(r['context_right_bp'])
    left_missing=max(0,want-left); right_missing=max(0,want-right)
    edge=min(s,clen-e)
    cstatus=('BOTH_COMPLETE' if left>=want and right>=want else
             'LEFT_TRUNCATED' if left<want and right>=want else
             'RIGHT_TRUNCATED' if left>=want and right<want else 'BOTH_TRUNCATED')
    potential_partial=edge<=2
    ctx=Path(r['genomic_context_fasta'])
    if not ctx.exists():
        # collected result archive path fallback
        ctx=b4a/'outputs/contexts_genomic'/f"{r['candidate_id']}.plusminus10kb.genomic.fna"
    _,seq=read_one_fasta(ctx)
    assert len(seq)==int(r['context_end0'])-int(r['context_start0']), (r['candidate_id'],len(seq),r['context_start0'],r['context_end0'])
    rr=dict(r)
    rr.update({
        'context_completeness_status':cstatus,
        'desired_flank_each_side_bp':want,
        'left_missing_bp':left_missing,
        'right_missing_bp':right_missing,
        'candidate_edge_distance_bp':edge,
        'potential_contig_edge_truncated_orf':str(potential_partial),
    })
    outrows.append(rr)
    if left<want and int(r['context_start0'])==0 and len(seq)>=31:
        anchors.append({'candidate_id':r['candidate_id'],'source_accession':r['source_accession'],'side':'LEFT','anchor31':seq[:31],'missing_bp':left_missing,'parent_contig_length_bp':clen})
    if right<want and int(r['context_end0'])==clen and len(seq)>=31:
        anchors.append({'candidate_id':r['candidate_id'],'source_accession':r['source_accession'],'side':'RIGHT','anchor31':seq[-31:],'missing_bp':right_missing,'parent_contig_length_bp':clen})

p=out/'F3B4A2_CONTEXT_COMPLETENESS_16.tsv'
with open(p,'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(outrows[0]),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(outrows)
apath=out/'F3B4A2_BOUNDARY_ANCHORS_31MER.tsv'
with open(apath,'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(anchors[0]),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(anchors)
acc=sorted({r['source_accession'] for r in rows})
(out/'F3B4A2_SOURCE_ACCESSIONS_12.txt').write_text('\n'.join(acc)+'\n')
summary={
 'stage':'Stage05F3B4A2-context-audit',
 'gate':'PASS_F3B4A2_CONTEXT_COMPLETENESS_AUDITED_PRE_ARCHITECTURE',
 'gate_pass':True,
 'survivor_n':len(rows),
 'source_accession_n':len(acc),
 'desired_flank_each_side_bp':want,
 'full_plusminus10kb_n':sum(int(r['context_left_bp'])>=want and int(r['context_right_bp'])>=want for r in rows),
 'both_sides_truncated_n':sum(r['context_completeness_status']=='BOTH_TRUNCATED' for r in outrows),
 'one_side_truncated_n':sum(r['context_completeness_status'] in {'LEFT_TRUNCATED','RIGHT_TRUNCATED'} for r in outrows),
 'candidate_within_2bp_of_contig_edge_n':sum(int(r['candidate_edge_distance_bp'])<=2 for r in outrows),
 'candidate_within_100bp_of_contig_edge_n':sum(int(r['candidate_edge_distance_bp'])<=100 for r in outrows),
 'parent_contig_lt_2kb_n':sum(int(r['parent_contig_length_bp'])<2000 for r in rows),
 'parent_contig_lt_5kb_n':sum(int(r['parent_contig_length_bp'])<5000 for r in rows),
 'parent_contig_ge_10kb_n':sum(int(r['parent_contig_length_bp'])>=10000 for r in rows),
 'architecture_absence_claim_authorized':False,
 'wetlab_freeze_authorized_from_current_contigs':False,
 'next_required_step':'source-context recovery before array/partner absence adjudication; all 16 retained',
 'b4a_gate_sha256':sha256(b4a/'outputs/STAGE05F3B4A_GATE.json'),
 'b4a_locus_map_sha256':sha256(mp),
 'context_table_sha256':sha256(p),
 'anchor_table_sha256':sha256(apath),
 'claim_boundary':'Context completeness audit only. Truncated source contigs cannot support absence-of-array/partner claims; potential contig-edge ORF truncation is technical, not biological negative evidence.'
}
json.dump(summary,open(out/'STAGE05F3B4A2_CONTEXT_AUDIT_GATE.json','w'),indent=2); open(out/'STAGE05F3B4A2_CONTEXT_AUDIT_GATE.json','a').write('\n')
print(json.dumps(summary,indent=2))
