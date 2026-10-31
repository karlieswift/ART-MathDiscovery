#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, json, hashlib, re
from collections import defaultdict, Counter

FORMAL_AA = {
"42a9674ed0f1bce773c92ae6d0ba7273acc94fa323dd0eb8947a6b8299a2804b",
"4330ef5e0b72ac1244f825df7df03355cbbba35736df7593a9c54afdbf977ff",
"8a5b3b69d22ae7b1575029d0f74ea022798f5f78b04c72ea79154d6b3c2efa00",
"955affef28083a9cde889643e7b90a24c752f0cd02d0948dd5626921ba335c70",
"d7186dbba2a9ae64d8b3a1c874b9303b8b65500c4ddcd0db4bd62f1be6921cd",
}

def read_tsv(p):
    with open(p, newline='') as f: return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(p, rows, fields):
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)

def sha256_bytes(b): return hashlib.sha256(b).hexdigest()
def sha256_file(p): return sha256_bytes(Path(p).read_bytes())

def fasta(path):
    name=None; seq=[]
    with open(path) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if name is not None: yield name,''.join(seq)
                name=line[1:].split()[0]; seq=[]
            else: seq.append(line)
        if name is not None: yield name,''.join(seq)

def pick(row, *names, default=''):
    for n in names:
        if n in row and row[n] not in (None,''): return row[n]
    return default

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('d1'); ap.add_argument('--out',default='outputs'); a=ap.parse_args()
    d1=Path(a.d1); out=Path(a.out); out.mkdir(exist_ok=True)
    candp=d1/'config/NOVELTY_FIRST_999_CANDIDATES.tsv'
    linkp=d1/'config/NOVELTY_FIRST_CANDIDATE_ORF_LINKS.tsv'
    orfp=d1/'config/NOVELTY_FIRST_24020_ORF_REGISTRY.tsv'
    fanp=d1/'config/NOVELTY_FIRST_PROTEIN_FANOUT.tsv'
    faap=d1/'inputs/NOVELTY_FIRST_UNIQUE_PROTEINS.faa'
    for p in [candp,linkp,orfp,fanp,faap]:
        if not p.exists(): raise SystemExit(f'MISSING {p}')
    cand=read_tsv(candp); links=read_tsv(linkp); orfs=read_tsv(orfp)
    seqs=dict(fasta(faap))
    # candidate key and stratum autodetect
    ckey_name = next((x for x in ['candidate_key','candidate_id','region_key'] if x in cand[0]), None)
    strat_name = next((x for x in ['novelty_architecture_stratum','novelty_stratum','stratum','architecture_stratum'] if x in cand[0]), None)
    if not ckey_name:
        # fallback: candidate keys look like ACCESSION::AR####
        ckey_name = next((k for k in cand[0] if any('::AR' in str(r.get(k,'')) for r in cand[:20])), None)
    if not strat_name:
        # fallback: N1/N2/N3/N4 architecture labels
        strat_name = next((k for k in cand[0] if any(re.match(r'^N[1-4](?:_|$)', str(r.get(k,''))) for r in cand[:50])), None)
    if not ckey_name or not strat_name: raise SystemExit(f'Cannot detect candidate key/stratum columns: {cand[0].keys()}')
    n12=[r for r in cand if str(r[strat_name]).startswith('N1') or str(r[strat_name]).startswith('N2')]
    counts=Counter('N1' if str(r[strat_name]).startswith('N1') else 'N2' for r in n12)
    if counts['N1']!=262 or counts['N2']!=340 or len(n12)!=602:
        raise SystemExit(f'FAIL N1/N2 counts {counts} total={len(n12)} expected 262/340/602')
    keys={r[ckey_name] for r in n12}
    # link columns
    if not links: raise SystemExit('empty links')
    lckey=next((x for x in ['candidate_key','candidate_id','region_key'] if x in links[0]),None)
    if not lckey:
        lckey=next((k for k in links[0] if sum(str(r.get(k,'')) in keys for r in links[:500])>0),None)
    # ORF registry key + aa id; infer before link ORF column so schemas can vary safely
    okey=next((x for x in ['orf_record_key','record_key','orf_id','protein_record_key'] if x in orfs[0]),None)
    if not okey:
        okey=next((k for k in orfs[0] if any('__PRODIGAL_' in str(r.get(k,'')) for r in orfs[:100])),None)
    aacol=next((x for x in ['aa_id','aa_sha256','protein_sha256'] if x in orfs[0]),None)
    if not aacol:
        hex64=re.compile(r'^(?:AA_)?[0-9a-fA-F]{64}$')
        aacol=next((k for k in orfs[0] if any(hex64.match(str(r.get(k,''))) for r in orfs[:100])),None)
    if not okey or not aacol: raise SystemExit(f'Cannot detect ORF columns: {orfs[0].keys()}')
    registry_keys={r[okey] for r in orfs}
    lorf=next((x for x in ['orf_record_key','record_key','orf_id','protein_record_key'] if x in links[0]),None)
    if not lorf:
        lorf=next((k for k in links[0] if sum(str(r.get(k,'')) in registry_keys for r in links[:500])>0),None)
    if not lckey or not lorf: raise SystemExit(f'Cannot detect link columns: {links[0].keys()}')
    sel_links=[r for r in links if r.get(lckey) in keys]
    orfkeys={r[lorf] for r in sel_links}
    omap={r[okey]:r for r in orfs}
    missing=sorted(orfkeys-set(omap))
    if missing: raise SystemExit(f'Missing {len(missing)} linked ORFs in registry; first={missing[:3]}')
    aa_to_orfs=defaultdict(list); aa_to_candidates=defaultdict(set); aa_to_strata=defaultdict(set)
    cstr={r[ckey_name]: ('N1' if str(r[strat_name]).startswith('N1') else 'N2') for r in n12}
    for lr in sel_links:
        ok=lr[lorf]; aa=omap[ok][aacol]; aa=aa.replace('AA_','') if aa.startswith('AA_') else aa
        aa_to_orfs[aa].append(ok); aa_to_candidates[aa].add(lr[lckey]); aa_to_strata[aa].add(cstr[lr[lckey]])
    # sequence lookup supports AA_<sha> or bare sha headers
    def getseq(aa):
        for k in [f'AA_{aa}',aa]:
            if k in seqs: return seqs[k]
        return ''
    rows=[]; hneg=[]
    outfa=out/'STAGE05D4_N1N2_UNIQUE_PROTEINS.faa'
    outfa_hneg=out/'STAGE05D4_N1N2_NONFORMAL_BROAD_RT_PROTEINS.faa'
    with open(outfa,'w') as f_all, open(outfa_hneg,'w') as f_neg:
        for aa in sorted(aa_to_orfs):
            s=getseq(aa)
            if not s: raise SystemExit(f'AA sequence missing for {aa}')
            formal=aa in FORMAL_AA
            # descriptive QC only
            acidic=sum(c in 'DE' for c in s)/len(s) if s else 0
            glypro=sum(c in 'GP' for c in s)/len(s) if s else 0
            motif_yxdd=bool(re.search(r'Y.DD',s))
            motif_dxd=bool(re.search(r'D.{0,3}D',s))
            row={
              'aa_sha256':aa,'aa_length':len(s),'n_orf_instances':len(aa_to_orfs[aa]),
              'n_candidate_links':len(aa_to_candidates[aa]),'candidate_keys':';'.join(sorted(aa_to_candidates[aa])),
              'strata':';'.join(sorted(aa_to_strata[aa])),'formal_stage05d1_broad_rt_positive':str(formal),
              'acidic_fraction':f'{acidic:.6f}','gly_pro_fraction':f'{glypro:.6f}',
              'YxDD_present':str(motif_yxdd),'D_x0_3_D_present':str(motif_dxd)
            }
            rows.append(row)
            hdr=f'>AA_{aa} n_orf={row["n_orf_instances"]} n_candidates={row["n_candidate_links"]} strata={row["strata"]} formal_broad_rt={formal}'
            f_all.write(hdr+'\n'+s+'\n')
            if not formal:
                f_neg.write(hdr+'\n'+s+'\n'); hneg.append(row)
    fields=list(rows[0])
    write_tsv(out/'STAGE05D4_N1N2_PROTEIN_UNIVERSE.tsv',rows,fields)
    # fanout table from exact candidate-link mapping
    fan=[]
    for aa in sorted(aa_to_orfs):
        for ok in sorted(set(aa_to_orfs[aa])):
            fan.append({'aa_sha256':aa,'orf_record_key':ok,'candidate_keys':';'.join(sorted({r[lckey] for r in sel_links if r[lorf]==ok}))})
    write_tsv(out/'STAGE05D4_EXACT_AA_FANOUT.tsv',fan,['aa_sha256','orf_record_key','candidate_keys'])
    # priority bins are descriptive only; no candidate filtering
    def lb(n):
        if n<80:return '<80'
        if n<150:return '80-149'
        if n<250:return '150-249'
        if n<400:return '250-399'
        return '>=400'
    summ={
      'stage':'Stage05D4','gate':'PASS_N1N2_REMOTE_HOMOLOGY_UNIVERSE_FROZEN','gate_pass':True,
      'n_N1_candidates':counts['N1'],'n_N2_candidates':counts['N2'],'n_N1N2_candidates':len(n12),
      'n_candidate_orf_links':len(sel_links),'n_linked_orf_instances':len(orfkeys),'n_unique_aa':len(rows),
      'n_formal_stage05d1_broad_rt_positive_unique_aa':sum(r['formal_stage05d1_broad_rt_positive']=='True' for r in rows),
      'n_nonformal_broad_rt_unique_aa':len(hneg),
      'aa_length_bins':dict(Counter(lb(int(r['aa_length'])) for r in rows)),
      'dna_score_recomputed':False,'dna_rank_recomputed':False,'candidate_filtering_performed':False,
      'formal_broad_rt_used_as_negative_filter':False,
      'claim_boundary':'Stage05D4 freezes the full N1/N2 candidate-linked protein universe for later remote-homology/structure analysis. No RT or ART identity is assigned here.'
    }
    (out/'GATE_STAGE05D4.json').write_text(json.dumps(summ,indent=2)+'\n')
    manifest=[]
    for p in sorted(out.iterdir()):
        if p.is_file(): manifest.append({'file':str(p),'bytes':p.stat().st_size,'sha256':sha256_file(p)})
    write_tsv(out/'STAGE05D4_OUTPUT_MANIFEST.tsv',manifest,['file','bytes','sha256'])
    print(json.dumps(summ,indent=2))

if __name__=='__main__': main()
