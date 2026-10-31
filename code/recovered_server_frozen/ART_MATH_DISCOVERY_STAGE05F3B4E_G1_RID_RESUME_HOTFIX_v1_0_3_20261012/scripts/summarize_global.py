#!/usr/bin/env python3
import csv,json,sys,re,datetime
from pathlib import Path
from common import *

def load_hits(p):
    if not Path(p).exists() or Path(p).stat().st_size==0: return []
    fields=['qseqid','saccver','pident','length','qstart','qend','sstart','send','evalue','bitscore','qlen','slen','qcovs','stitle','sscinames']
    rows=[]
    with open(p,encoding='utf-8',errors='ignore') as f:
        for line in f:
            z=line.rstrip('\n').split('\t')
            if len(z)<13: continue
            if len(z)<len(fields): z += ['']*(len(fields)-len(z))
            r=dict(zip(fields,z[:len(fields)])); rows.append(r)
    return rows

def summarize_protein(hits,qids):
    by={q:[] for q in qids}
    for r in hits:
        if r['qseqid'] in by: by[r['qseqid']].append(r)
    out={}
    for q,xs in by.items():
        for r in xs:
            r['_pid']=fnum(r['pident'])/100.; r['_qcov']=fnum(r['qcovs'])/100.; r['_tcov']=fnum(r['length'])/max(1,fnum(r['slen'])); r['_bits']=fnum(r['bitscore']); r['_e']=fnum(r['evalue'],1e99)
        xs=sorted(xs,key=lambda r:(r['_e'],-r['_bits'],-r['_qcov'],-r['_pid'],r['saccver']))
        top=xs[0] if xs else None
        exact=any(r['_pid']>=0.99999 and r['_qcov']>=0.99 and r['_tcov']>=0.99 for r in xs)
        high=any(r['_pid']>=0.90 and r['_qcov']>=0.80 and r['_tcov']>=0.80 for r in xs)
        explicit_rt=any(re.search(r'reverse transcriptase|rna-directed dna polymerase|retrotranspos|group ii intron|retron',r.get('stitle',''),re.I) for r in xs[:20])
        explicit_art=any(re.search(r'array[- ]associated reverse transcriptase',r.get('stitle',''),re.I) for r in xs[:50])
        out[q]={'hit_n':len(xs),'exact':exact,'high':high,'explicit_rt_top20':explicit_rt,'explicit_art_top50':explicit_art,
                'top_saccver':top['saccver'] if top else '','top_pident':top['_pid'] if top else '','top_qcov':top['_qcov'] if top else '','top_tcov':top['_tcov'] if top else '',
                'top_evalue':top['evalue'] if top else '','top_bitscore':top['bitscore'] if top else '','top_title':top['stitle'] if top else '','top_species':top['sscinames'] if top else ''}
    return out

def summarize_arch(hits,qids):
    by={q:[] for q in qids}
    for r in hits:
        if r['qseqid'] in by: by[r['qseqid']].append(r)
    out={}
    for q,xs in by.items():
        for r in xs:
            r['_pid']=fnum(r['pident'])/100.; r['_qcov']=fnum(r['qcovs'])/100.; r['_bits']=fnum(r['bitscore']); r['_e']=fnum(r['evalue'],1e99)
        xs=sorted(xs,key=lambda r:(r['_e'],-r['_bits'],-r['_qcov'],-r['_pid'],r['saccver']))
        top=xs[0] if xs else None
        exact=any(r['_pid']>=0.999 and r['_qcov']>=0.95 for r in xs)
        near=any(r['_pid']>=0.95 and r['_qcov']>=0.80 for r in xs)
        out[q]={'hit_n':len(xs),'exact_like':exact,'near_full':near,'top_saccver':top['saccver'] if top else '','top_pident':top['_pid'] if top else '','top_qcov':top['_qcov'] if top else '',
                'top_evalue':top['evalue'] if top else '','top_bitscore':top['bitscore'] if top else '','top_title':top['stitle'] if top else '','top_species':top['sscinames'] if top else ''}
    return out

def main():
    if len(sys.argv)!=8: raise SystemExit('usage summarize_global.py targets rt_hits partner_hits arch_hits out lock metadata_json')
    targ,rtp,pp,ap,outp,lockp,metap=sys.argv[1:]
    targets=read_tsv(targ); fam_by_rt={r['rt_candidate_id']:r['family_rep_id'] for r in targets}; fam_by_partner={r['partner_orf_id']:r['family_rep_id'] for r in targets}
    rt=summarize_protein(load_hits(rtp),list(fam_by_rt)); pr=summarize_protein(load_hits(pp),list(fam_by_partner)); ar=summarize_arch(load_hits(ap),list(fam_by_rt))
    rows=[]
    for t in targets:
        fam=t['family_rep_id']; q=t['rt_candidate_id']; pq=t['partner_orf_id']; R=rt[q]; P=pr[pq]; A=ar[q]
        state='PUBLIC_HIGH_ID_SEQUENCE_OR_LOCUS_PRESENT' if (R['high'] or A['near_full']) else 'NO_HIGH_ID_NEAR_FULL_MATCH_IN_NCBI_NR_NT_AUDIT'
        z={'family_rep_id':fam,'rt_candidate_id':q,'global_public_audit_state':state,
           'rt_nr_hit_n':R['hit_n'],'rt_nr_exact_full_length_match':R['exact'],'rt_nr_high_identity_near_full_match':R['high'],'rt_nr_explicit_rt_annotation_top20':R['explicit_rt_top20'],'rt_nr_explicit_art_annotation_top50':R['explicit_art_top50'],
           'rt_nr_top_saccver':R['top_saccver'],'rt_nr_top_pident':R['top_pident'],'rt_nr_top_qcov':R['top_qcov'],'rt_nr_top_tcov':R['top_tcov'],'rt_nr_top_title':R['top_title'],'rt_nr_top_species':R['top_species'],
           'partner_nr_hit_n':P['hit_n'],'partner_nr_exact_full_length_match':P['exact'],'partner_nr_high_identity_near_full_match':P['high'],'partner_nr_top_saccver':P['top_saccver'],'partner_nr_top_pident':P['top_pident'],'partner_nr_top_qcov':P['top_qcov'],'partner_nr_top_tcov':P['top_tcov'],'partner_nr_top_title':P['top_title'],'partner_nr_top_species':P['top_species'],
           'architecture_nt_hit_n':A['hit_n'],'architecture_nt_exact_like_match':A['exact_like'],'architecture_nt_near_full_match':A['near_full'],'architecture_nt_top_saccver':A['top_saccver'],'architecture_nt_top_pident':A['top_pident'],'architecture_nt_top_qcov':A['top_qcov'],'architecture_nt_top_title':A['top_title'],'architecture_nt_top_species':A['top_species']}
        rows.append(z)
    write_tsv(outp,rows)
    meta=json.load(open(metap)) if Path(metap).exists() else {}
    gate={'stage':'Stage05F3B4E-G1','gate':'PASS_F3B4E_DATED_NCBI_NR_NT_GLOBAL_PUBLIC_AUDIT_COMPLETE','gate_pass':True,'target_family_n':5,
          'nr_rt_high_identity_near_full_family_n':sum(b(r['rt_nr_high_identity_near_full_match']) for r in rows),'nt_architecture_near_full_family_n':sum(b(r['architecture_nt_near_full_match']) for r in rows),
          'absolute_global_novelty_claim_authorized':False,'database_scope':'NCBI nr/nt remote BLAST; public database audit, not every metagenomic/proprietary sequence corpus',
          'runtime_metadata':meta,'summary_sha256':sha256_file(outp),'claim_boundary':'No high-identity match in nr/nt is a dated negative search result, not proof of absolute global novelty.'}
    gp=Path(outp).parent/'STAGE05F3B4E_G1_GATE.json'; gp.write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate,indent=2)); print('\n=== GLOBAL PUBLIC AUDIT 5 ==='); print(Path(outp).read_text().rstrip())

if __name__=='__main__': main()
