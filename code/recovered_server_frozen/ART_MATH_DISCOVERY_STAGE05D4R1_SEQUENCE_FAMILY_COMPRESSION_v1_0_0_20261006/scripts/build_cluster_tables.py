#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, json, hashlib
from Bio import SeqIO

def read_tsv(p):
    with open(p,newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):
            h.update(b)
    return h.hexdigest()

def parse_cluster_tsv(p):
    rows=[]
    with open(p) as f:
        for line in f:
            if not line.strip():
                continue
            rep,member=line.rstrip('\n').split('\t')[:2]
            rows.append((rep,member))
    return rows

def raw_sha(aa_id):
    return aa_id[3:] if aa_id.startswith('AA_') else aa_id

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--d4',required=True)
    ap.add_argument('--work',required=True)
    ap.add_argument('--out',required=True)
    a=ap.parse_args()
    d4=Path(a.d4); work=Path(a.work); out=Path(a.out)
    out.mkdir(parents=True,exist_ok=True)

    gate=json.loads((d4/'outputs/GATE_STAGE05D4.json').read_text())
    assert gate['gate']=='PASS_N1N2_REMOTE_HOMOLOGY_UNIVERSE_FROZEN'
    assert gate['n_N1N2_candidates']==602 and gate['n_nonformal_broad_rt_unique_aa']==13062

    rows=read_tsv(d4/'outputs/STAGE05D4_N1N2_PROTEIN_UNIVERSE.tsv')
    if not rows:
        raise SystemExit('Empty Stage05D4 protein universe')
    if 'aa_sha256' not in rows[0]:
        raise SystemExit(f"aa_sha256 missing from universe: {rows[0].keys()}")
    universe={r['aa_sha256']:r for r in rows}

    seqs={r.id:str(r.seq) for r in SeqIO.parse(d4/'outputs/STAGE05D4_N1N2_NONFORMAL_BROAD_RT_PROTEINS.faa','fasta')}
    assert len(seqs)==13062

    # MMseqs preserves FASTA record IDs, which are AA_<sha256>, whereas
    # STAGE05D4_N1N2_PROTEIN_UNIVERSE.tsv stores bare aa_sha256.  Normalize
    # only at this join boundary; scientific inputs, clusters and thresholds
    # are unchanged.
    missing=[x for x in seqs if raw_sha(x) not in universe]
    if missing:
        raise SystemExit(f'FASTA IDs not recoverable in universe after AA_ normalization: n={len(missing)} first={missing[:3]}')

    for tier in ['nr90','fam50']:
        ctsv=work/tier/f'{tier}_cluster.tsv'
        pairs=parse_cluster_tsv(ctsv)
        byrep={}
        for rep,mem in pairs:
            byrep.setdefault(rep,[]).append(mem)
        members=set(m for _,m in pairs)
        assert members==set(seqs), (tier,len(members),len(seqs),len(set(seqs)-members))

        tab=out/f'STAGE05D4R1_{tier.upper()}_CLUSTERS.tsv'
        with tab.open('w',newline='') as f:
            cols=['cluster_id','representative_aa_sha256','member_aa_sha256','member_aa_length','member_candidate_keys','member_strata','member_YxDD_present','member_D_x0_3_D_present']
            w=csv.DictWriter(f,fieldnames=cols,delimiter='\t')
            w.writeheader()
            for i,rep in enumerate(sorted(byrep),1):
                cid=f'{tier.upper()}_{i:05d}'
                rep_sha=raw_sha(rep)
                for mem in sorted(byrep[rep]):
                    mem_sha=raw_sha(mem)
                    u=universe[mem_sha]
                    w.writerow(dict(
                        cluster_id=cid,
                        representative_aa_sha256=rep_sha,
                        member_aa_sha256=mem_sha,
                        member_aa_length=u['aa_length'],
                        member_candidate_keys=u['candidate_keys'],
                        member_strata=u['strata'],
                        member_YxDD_present=u['YxDD_present'],
                        member_D_x0_3_D_present=u['D_x0_3_D_present']))

        repfa=out/f'STAGE05D4R1_{tier.upper()}_REPRESENTATIVES.faa'
        with repfa.open('w') as f:
            for rep in sorted(byrep):
                f.write(f'>{rep}\n{seqs[rep]}\n')

        summary={
            'tier':tier,
            'n_input_unique_aa':len(seqs),
            'n_clusters':len(byrep),
            'compression_ratio':len(seqs)/len(byrep),
            'singleton_clusters':sum(len(v)==1 for v in byrep.values()),
            'largest_cluster_size':max(map(len,byrep.values())),
            'cluster_table_sha256':sha(tab),
            'representative_fasta_sha256':sha(repfa),
            'id_namespace_hotfix':'MMseqs AA_<sha256> normalized to bare aa_sha256 for universe join only'
        }
        (out/f'STAGE05D4R1_{tier.upper()}_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
        print(json.dumps(summary,indent=2))

    g={
      'stage':'Stage05D4-R1',
      'gate':'PASS_N1N2_SEQUENCE_FAMILY_COMPRESSION_FROZEN',
      'gate_pass':True,
      'source_stage05d4_gate':gate['gate'],
      'n_N1N2_candidates':602,
      'n_input_nonformal_unique_aa':13062,
      'dna_score_recomputed':False,
      'dna_rank_recomputed':False,
      'candidate_filtering_performed':False,
      'nr90_used_for_redundancy_compression_only':True,
      'fam50_used_as_descriptive_family_map_only':True,
      'implementation_hotfix':'v1.0.1 AA_ prefix namespace normalization only; clustering thresholds and MMseqs outputs unchanged',
      'claim_boundary':'Sequence clustering compresses redundant protein queries and maps broad sequence families; it does not assign RT, ART, or novelty identity.'
    }
    (out/'GATE_STAGE05D4R1.json').write_text(json.dumps(g,indent=2)+'\n')

if __name__=='__main__':
    main()
