#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, datetime, json, shutil
from common import parse_fasta, safe_id, accession_root, sha256_file, write_json

EXCLUDE_ROOTS={'MW248466','MW349129','MW349128','MW218148'}
ALLOWED_ARMS={'L4_CURATED','L4_PROSPECTIVE','L4_METAVIROME'}
CUTOFF=datetime.date(2026,9,29)

def parse_date(s):
    if not s: return None
    x=str(s).strip()[:10]
    try: return datetime.date.fromisoformat(x)
    except Exception: raise ValueError(f'invalid ISO release_date {s!r}; expected YYYY-MM-DD')

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); ap.add_argument('--cohort-tsv',required=True); a=ap.parse_args()
    root=Path(a.root).resolve(); engkids=[p for p in (root/'runtime_engine').iterdir() if p.is_dir()];
    if len(engkids)!=1: raise SystemExit('run preflight first')
    eng=engkids[0]; (root/'inputs').mkdir(exist_ok=True); (eng/'inputs').mkdir(exist_ok=True)
    with open(a.cohort_tsv) as f: rows=list(csv.DictReader(f,delimiter='\t'))
    if not rows: raise SystemExit('empty Stage04 cohort')
    required={'sequence_id','arm','source','release_date','fasta_path'}
    if not required.issubset(rows[0]): raise SystemExit('cohort missing required columns: '+','.join(sorted(required-set(rows[0]))))
    out=[]; seen=set()
    for r in rows:
        sid=r['sequence_id'].strip(); arm=r['arm'].strip();
        if arm not in ALLOWED_ARMS: raise SystemExit(f'bad arm {arm} for {sid}')
        if accession_root(sid) in EXCLUDE_ROOTS: raise SystemExit(f'known development/validation accession forbidden in Stage04 primary cohort: {sid}')
        rd=parse_date(r.get('release_date',''))
        if arm=='L4_PROSPECTIVE' and (rd is None or rd<CUTOFF): raise SystemExit(f'prospective row {sid} has release_date {r.get("release_date","")} before 2026-09-29')
        aid=safe_id(sid)
        if aid in seen: raise SystemExit(f'duplicate analysis_id after sanitization: {aid}')
        seen.add(aid); src=Path(r['fasta_path']).expanduser().resolve(); header,seq=parse_fasta(src)
        if len(seq)<4096: raise SystemExit(f'{sid}: length <4096')
        extra=set(seq)-set('ACGT')
        if extra: raise SystemExit(f'{sid}: non-ACGT characters {sorted(extra)}; frozen engine has no N policy')
        dst=root/'inputs'/f'{aid}.frozen.fasta'; edst=eng/'inputs'/f'{aid}.frozen.fasta'
        txt='>'+header+'\n'+'\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n'; dst.write_text(txt); shutil.copy2(dst,edst)
        out.append({'analysis_id':aid,'sequence_id':sid,'arm':arm,'source':r.get('source',''),'release_date':r.get('release_date',''),'cluster_id':r.get('cluster_id',''),'notes':r.get('notes',''),'length_bp':len(seq),'fasta_sha256':sha256_file(dst),'fasta_header':header,'frozen_fasta':str(dst.relative_to(root))})
    mp=root/'outputs/STAGE04A_UNIVERSE_MANIFEST.tsv'
    fields=['analysis_id','sequence_id','arm','source','release_date','cluster_id','notes','length_bp','fasta_sha256','fasta_header','frozen_fasta']
    with open(mp,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(out)
    gate={'gate':'PASS_STAGE04A_UNIVERSE_FROZEN','gate_pass':True,'n_sequences':len(out),'arms':sorted(set(x['arm'] for x in out)),'known_accessions_excluded':True,'prospective_cutoff_utc':'2026-09-29T00:00:00Z','truth_used':False,'universe_manifest_sha256':sha256_file(mp)}
    write_json(root/'outputs/STAGE04A_UNIVERSE_GATE.json',gate); print(json.dumps(gate,indent=2))
if __name__=='__main__': main()
