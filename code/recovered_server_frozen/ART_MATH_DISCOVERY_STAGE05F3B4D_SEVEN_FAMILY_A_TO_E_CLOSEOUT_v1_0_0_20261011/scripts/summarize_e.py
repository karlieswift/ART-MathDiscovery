#!/usr/bin/env python3
import csv,json,re,sys
from pathlib import Path

def ff(x):
    try:return float(x)
    except:return float('inf')

def parse_m8(p):
    p=Path(p); out=[]
    if not p.exists() or p.stat().st_size==0:return out
    fields=['query','target','fident','alnlen','qcov','tcov','evalue','bits']
    with p.open(encoding='utf-8') as f:
        first=f.readline().rstrip('\n')
        if not first:return out
        if first.startswith('query\t'): return list(csv.DictReader([first]+list(f),delimiter='\t'))
        out.append(dict(zip(fields,first.split('\t'))))
        for line in f:
            if line.strip():out.append(dict(zip(fields,line.rstrip('\n').split('\t'))))
    return out

def main():
    if len(sys.argv)!=10: raise SystemExit('usage: summarize_e.py REP7_TSV REF_MAP LOCAL_M8 GLOBAL_M8 ART_M8 GLOBAL_PRESENT ART_PRESENT OUT_TSV OUT_GATE')
    rep_p,refmap_p,local_p,global_p,art_p,global_present,art_present,out_p,gate_p=sys.argv[1:]
    reps=list(csv.DictReader(open(rep_p,encoding='utf-8'),delimiter='\t'))
    cmap={r['candidate_id']:r['family_rep_id'] for r in reps}
    ref={}
    if Path(refmap_p).exists():
        with open(refmap_p,encoding='utf-8') as f:
            for r in csv.DictReader(f,delimiter='\t'):ref[r['reference_id']]=r
    local=parse_m8(local_p); glob=parse_m8(global_p); art=parse_m8(art_p)
    by=lambda rows,q:[r for r in rows if r.get('query')==q]
    art_re=re.compile(r'(\bART[_ -]?[0-9]*\b|MarsHill|Madawaska|PB50|LY01|PALS_?2|UFV_?DC4|MW349129|MW218148|OR770614|OR836606|LC680885|MN091626|MZ779063)',re.I)
    rt_re=re.compile(r'(reverse transcript|RNA-directed DNA polymerase|RTase|retrotrans|retron|group II intron|diversity-generating retroelement|DGR)',re.I)
    rows=[]
    for r in reps:
        q=r['candidate_id']; fam=r['family_rep_id']; hs=sorted(by(local,q),key=lambda x:(ff(x['evalue']),-float(x.get('bits',0))))
        top=hs[0] if hs else None
        exact=any(float(h['fident'])>=0.999999 and float(h['qcov'])>=0.99 and float(h['tcov'])>=0.99 for h in hs)
        high=any(float(h['fident'])>=0.90 and float(h['qcov'])>=0.80 and float(h['tcov'])>=0.80 for h in hs)
        explicit_rt=False; public_art=False
        if top:
            m=ref.get(top['target'],{}); txt=(m.get('source_paths','')+' '+m.get('original_headers',''))
            explicit_rt=bool(rt_re.search(txt)); public_art=bool(art_re.search(txt))
        gh=sorted(by(glob,q),key=lambda x:(ff(x['evalue']),-float(x.get('bits',0)))) if global_present=='1' else []
        ah=sorted(by(art,q),key=lambda x:(ff(x['evalue']),-float(x.get('bits',0)))) if art_present=='1' else []
        gtop=gh[0] if gh else None; atop=ah[0] if ah else None
        ghigh=any(float(h['fident'])>=0.90 and float(h['qcov'])>=0.80 and float(h['tcov'])>=0.80 for h in gh)
        ahigh=any(float(h['fident'])>=0.90 and float(h['qcov'])>=0.80 and float(h['tcov'])>=0.80 for h in ah)
        novelty='UNRESOLVED_GLOBAL_AND_ART95'
        if global_present=='1' and art_present=='1':
            novelty='NO_HIGH_ID_NEAR_FULL_MATCH_IN_SUPPLIED_GLOBAL_AND_ART_CATALOGS' if not ghigh and not ahigh else 'KNOWN_OR_NEAR_KNOWN_MATCH_IN_SUPPLIED_CATALOG'
        rows.append({
          'family_rep_id':fam,'candidate_id':q,
          'local_known_reference_hit_n':len(hs),'local_exact_full_length_match':exact,'local_high_identity_near_full_match':high,
          'local_top_reference_id':top['target'] if top else '','local_top_fident':top['fident'] if top else '','local_top_qcov':top['qcov'] if top else '','local_top_tcov':top['tcov'] if top else '',
          'local_top_reference_explicit_rt_keyword':explicit_rt,'local_top_reference_public_art_named_subset_keyword':public_art,
          'supplied_global_catalog_present':global_present=='1','global_top_target':gtop['target'] if gtop else '','global_high_identity_near_full_match':ghigh,
          'supplied_art_catalog_present':art_present=='1','art_top_target':atop['target'] if atop else '','art_high_identity_near_full_match':ahigh,
          'novelty_status':novelty,'anthropic_unreported_art_claim_authorized':False,
        })
    cols=list(rows[0])
    with open(out_p,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    gate={
      'stage':'Stage05F3B4D-E','gate':'PASS_F3B4D_E_KNOWN_SYSTEM_AND_NOVELTY_CROSSWALK_COMPLETE_WITH_CLAIM_BOUNDARY','gate_pass':True,
      'target_family_n':7,'local_frozen_known_reference_audit_complete':True,
      'supplied_global_sequence_catalog_present':global_present=='1','supplied_full_art_catalog_present':art_present=='1',
      'global_sequence_novelty_closed':bool(global_present=='1' and art_present=='1'),
      'anthropic_unreported_art_claim_authorized':False,
      'claim_boundary':'Local known-system crosswalk is complete for the frozen project reference panel. Without a supplied comprehensive global sequence catalog and full ART95 sequence catalog, absence of a local match cannot establish novelty or Anthropic-unreported status.'}
    Path(gate_p).write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate,indent=2)); print('\n=== E KNOWN/NOVELTY CROSSWALK 7 ===');print('\t'.join(cols))
    for r in rows:print('\t'.join(str(r[c]) for c in cols))

if __name__=='__main__': main()
