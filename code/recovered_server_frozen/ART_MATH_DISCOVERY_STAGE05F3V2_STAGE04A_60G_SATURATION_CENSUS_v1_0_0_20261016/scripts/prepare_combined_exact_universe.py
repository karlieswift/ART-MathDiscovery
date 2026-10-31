#!/usr/bin/env python3
from pathlib import Path
import argparse, gzip, hashlib, json, os, re, sys

EXTS={'.faa','.fa','.fasta','.fas','.fsa'}
KEYWORDS=('unique','exact','dedup','formal')
SKIP_PARTS={'raw','shards','runtime_engine','.git','__pycache__'}

def op(p):
    return gzip.open(p,'rt') if p.suffix=='.gz' else open(p,'rt')

def fasta_iter(p):
    h=None; seq=[]
    with op(p) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if h is not None:
                    yield h,''.join(seq).upper().replace(' ','')
                h=line[1:].strip(); seq=[]
            else:
                seq.append(line)
        if h is not None:
            yield h,''.join(seq).upper().replace(' ','')

def count_records(p):
    n=0
    with op(p) as f:
        for line in f:
            if line.startswith('>'): n+=1
    return n

def candidate_files(root):
    root=Path(root)
    if not root.exists(): return []
    out=[]
    for p in root.rglob('*'):
        if not p.is_file(): continue
        low=p.name.lower()
        parts={x.lower() for x in p.parts}
        if parts & SKIP_PARTS: continue
        suff=p.suffix.lower()
        isfa=suff in EXTS or (suff=='.gz' and p.with_suffix('').suffix.lower() in EXTS)
        if not isfa: continue
        if not any(k in low for k in KEYWORDS): continue
        try:
            size=p.stat().st_size
        except OSError:
            continue
        if size > 2_000_000_000: continue
        out.append(p)
    return sorted(out)

def choose(root, expected, label, fallback_root=None):
    scanned=[]
    for p in candidate_files(root):
        n=count_records(p); scanned.append((str(p),n))
    hits=[Path(p) for p,n in scanned if n==expected]
    if not hits and fallback_root is not None:
        for p in candidate_files(fallback_root):
            if str(p).startswith(str(root)): continue
            low=str(p).lower()
            if 'stage05f3' not in low: continue
            if label=='OLD' and 'f3v2' in low: continue
            if label=='DELTA' and 'f3v2' not in low: continue
            n=count_records(p); scanned.append((str(p),n))
            if n==expected: hits.append(p)
    if not hits:
        print(f'HOLD: no {label} FASTA with expected record count {expected}',file=sys.stderr)
        for p,n in scanned:
            print(f'CANDIDATE\t{n}\t{p}',file=sys.stderr)
        raise SystemExit(2)
    def score(p):
        x=p.name.lower(); s=0
        if 'unique' in x: s+=100
        if 'exact' in x: s+=50
        if 'dedup' in x: s+=20
        if 'formal' in x: s+=5
        if 'all_orfs' in x: s-=500
        s-=len(str(p))/10000
        return (-s, len(str(p)), str(p))
    hits=sorted(set(hits), key=score)
    chosen=hits[0]
    print(f'{label}_UNIQUE_FASTA={chosen}')
    if len(hits)>1:
        print(f'{label}_MATCHING_FASTA_N={len(hits)}')
        for p in hits[:20]: print(f'{label}_MATCH={p}')
    return chosen, scanned

def shaseq(seq): return hashlib.sha256(seq.encode()).hexdigest()
def shafile(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--project-root',required=True)
    ap.add_argument('--old-root',required=True)
    ap.add_argument('--delta-root',required=True)
    ap.add_argument('--out',required=True)
    ap.add_argument('--old-expected',type=int,required=True)
    ap.add_argument('--delta-expected',type=int,required=True)
    a=ap.parse_args()
    P=Path(a.project_root).resolve(); O=Path(a.old_root).resolve(); D=Path(a.delta_root).resolve(); OUT=Path(a.out).resolve(); OUT.mkdir(parents=True,exist_ok=True)
    oldp, oldscan=choose(O,a.old_expected,'OLD',P)
    deltap, deltascan=choose(D,a.delta_expected,'DELTA',None)

    origin={}
    examples={}
    counts={}
    for label,p,exp in [('OLD',oldp,a.old_expected),('DELTA',deltap,a.delta_expected)]:
        n=0; seen=set()
        for h,seq in fasta_iter(p):
            if not seq: raise RuntimeError(f'{label}: empty sequence: {h}')
            n+=1; s=shaseq(seq)
            seen.add(s)
            origin.setdefault(s,set()).add(label)
            examples.setdefault((s,label),h)
        if n!=exp: raise RuntimeError(f'{label}: record count changed {n}!={exp}')
        if len(seen)!=exp: raise RuntimeError(f'{label}: input is not exact-AA unique: unique={len(seen)} records={exp}')
        counts[label]=n

    oldset={s for s,o in origin.items() if 'OLD' in o}
    deltaset={s for s,o in origin.items() if 'DELTA' in o}
    overlap=oldset & deltaset
    union=sorted(origin)

    faa=OUT/'COMBINED_GLOBAL_EXACT_UNIQUE.faa'
    mp=OUT/'GLOBAL_EXACT_AA_MAP.tsv'
    # Need sequence lookup; reread both selected inputs.
    seq_by_sha={}
    for p in [oldp,deltap]:
        for _,seq in fasta_iter(p): seq_by_sha.setdefault(shaseq(seq),seq)

    with open(faa,'w') as fo, open(mp,'w') as fm:
        fm.write('global_id\tsequence_sha256\torigin_class\told_header_example\tdelta_header_example\taa_length\n')
        for i,s in enumerate(union,1):
            gid=f'G{i:08d}'
            oo=origin[s]
            cls='SHARED_EXACT' if oo=={'OLD','DELTA'} else ('OLD_ONLY' if 'OLD' in oo else 'DELTA_ONLY')
            seq=seq_by_sha[s]
            fo.write(f'>{gid}\n')
            for j in range(0,len(seq),80): fo.write(seq[j:j+80]+'\n')
            oh=examples.get((s,'OLD'),'').replace('\t',' ').replace('\n',' ')
            dh=examples.get((s,'DELTA'),'').replace('\t',' ').replace('\n',' ')
            fm.write(f'{gid}\t{s}\t{cls}\t{oh}\t{dh}\t{len(seq)}\n')

    meta={
      'stage':'Stage05F3V2-Stage04A-precluster',
      'old_unique_expected':a.old_expected,
      'delta_unique_expected':a.delta_expected,
      'old_unique_fasta':str(oldp),
      'old_unique_fasta_sha256':shafile(oldp),
      'delta_unique_fasta':str(deltap),
      'delta_unique_fasta_sha256':shafile(deltap),
      'old_exact_unique_n':len(oldset),
      'delta_exact_unique_n':len(deltaset),
      'cross_universe_exact_overlap_n':len(overlap),
      'combined_global_exact_unique_n':len(union),
      'old_only_exact_n':len(oldset-overlap),
      'delta_only_exact_n':len(deltaset-overlap),
      'shared_exact_n':len(overlap),
      'science_change':False,
      'retuning_authorized':False,
      'claim_boundary':'Analysis-only global exact-AA dedupe for saturation census; frozen old and delta discovery outputs are not rewritten.'
    }
    with open(OUT/'PRECLUSTER_EXACT_CENSUS.json','w') as f: json.dump(meta,f,indent=2); f.write('\n')
    print(json.dumps(meta,indent=2))

if __name__=='__main__': main()
