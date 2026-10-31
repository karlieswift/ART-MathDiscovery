#!/usr/bin/env python3
import gzip, hashlib, json, os, re, sys
from collections import defaultdict
from pathlib import Path

P=Path(sys.argv[1])
outdir=Path(sys.argv[2])
outdir.mkdir(parents=True,exist_ok=True)

FA_EXT=(".fa",".faa",".fasta",".fas",".fa.gz",".faa.gz",".fasta.gz",".fas.gz")
role_re=re.compile(r'(rvt|reverse.?transcript|(^|[^a-z])rt([^a-z]|$)|(^|[^a-z])art([^a-z]|$)|swiss|known|positive|calibr)',re.I)

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(16<<20),b''): h.update(b)
    return h.hexdigest()

def opener(path):
    return gzip.open(path,'rt',encoding='utf-8',errors='ignore') if str(path).lower().endswith('.gz') else open(path,'rt',encoding='utf-8',errors='ignore')

def read_fasta(path):
    h=None; seq=[]
    with opener(path) as f:
        for line in f:
            line=line.rstrip('\r\n')
            if line.startswith('>'):
                if h is not None: yield h,''.join(seq).upper()
                h=line[1:]; seq=[]
            elif line.strip(): seq.append(line.strip())
        if h is not None: yield h,''.join(seq).upper()

def looks_protein(path):
    n=0; aa=0; non_dna=0
    try:
        for h,s in read_fasta(path):
            s=re.sub('[^A-Z*]','',s)
            aa += len(s); non_dna += sum(c not in 'ACGTUN' for c in s)
            n += 1
            if n>=20 or aa>=20000: break
    except Exception:
        return False
    if n==0 or aa==0: return False
    return (non_dna/aa)>=0.01

explicit=os.environ.get('F3B3B_REFERENCE_FASTAS','').strip()
candidates=[]
mode='AUTO_FROZEN_PROJECT_REFERENCE_SCAN'
if explicit:
    mode='EXPLICIT_ENV_REFERENCE_FASTAS'
    for x in explicit.split(':'):
        if x.strip(): candidates.append(Path(x.strip()).resolve())
else:
    stage_globs=[
        'ART_MATH_DISCOVERY_STAGE05D1_*',
        'ART_MATH_DISCOVERY_STAGE05D4R5B_*',
        'ART_MATH_DISCOVERY_STAGE05D4R5C_*',
        'ART_MATH_DISCOVERY_STAGE05E3*',
        'ART_MATH_DISCOVERY_STAGE05E4*',
    ]
    roots=[]
    for pat in stage_globs: roots.extend(sorted(P.glob(pat)))
    for root in roots:
        if not root.is_dir(): continue
        for path in root.rglob('*'):
            if not path.is_file(): continue
            low=path.name.lower()
            if not any(low.endswith(x) for x in FA_EXT): continue
            rel=str(path).lower()
            if '/references/' not in rel and '/reference/' not in rel: continue
            if not role_re.search(str(path)): continue
            candidates.append(path.resolve())

# Deterministic unique paths, excluding current F3B stages.
uniq=[]; seen=set()
for p in sorted(candidates,key=lambda x:str(x)):
    if p in seen: continue
    seen.add(p)
    if not p.is_file() or p.stat().st_size==0: continue
    if 'STAGE05F3B' in str(p): continue
    if not looks_protein(p): continue
    uniq.append(p)

# Existing D4R5B RVT-All/SwissProt hit IDs, if available, are used to focus any full SwissProt source.
swiss_hit_ids=set()
for d in sorted(P.glob('ART_MATH_DISCOVERY_STAGE05D4R5B_RT_FOCUSED_STRUCTURAL_PANEL_*')):
    t=d/'raw/RVTALL_SWISSPROT.tblout'
    if t.is_file():
        for line in open(t,encoding='utf-8',errors='ignore'):
            if not line.strip() or line.startswith('#'): continue
            x=line.split()
            if x: swiss_hit_ids.add(x[0])

inventory=[]
seq_sources=defaultdict(list)
header_sources=defaultdict(list)
source_counts={}
for path in uniq:
    src_n=0; kept_n=0
    is_swiss='swiss' in str(path).lower()
    for h,s in read_fasta(path):
        src_n += 1
        s=re.sub('[^A-Z]','',s.replace('*',''))
        if not s: continue
        token=h.split()[0]
        # When a full SwissProt source is present and prior frozen RVT-All hits exist,
        # restrict it to the pre-existing RT-focused positive-calibration IDs.
        if is_swiss and swiss_hit_ids and token not in swiss_hit_ids:
            continue
        kept_n += 1
        sh=hashlib.sha256(s.encode('ascii')).hexdigest()
        seq_sources[sh].append(str(path))
        header_sources[sh].append(h.replace('\t',' '))
        if sh not in source_counts:
            source_counts[sh]=s
    inventory.append({
        'path':str(path),'sha256':sha(path),'source_sequence_n':src_n,
        'retained_reference_sequence_n':kept_n,
        'swissprot_focused_by_prior_rvtall_hits':bool(is_swiss and swiss_hit_ids)
    })

panel=outdir/'F3B3B_KNOWN_SYSTEM_REFERENCE_PANEL.faa'
mapf=outdir/'F3B3B_KNOWN_SYSTEM_REFERENCE_MAP.tsv'
with open(panel,'w') as fo, open(mapf,'w') as fm:
    fm.write('reference_id\tsequence_sha256\taa_length\tsource_path_n\tsource_paths\toriginal_headers\n')
    for i,sh in enumerate(sorted(source_counts),1):
        rid=f'KREF{i:08d}'
        s=source_counts[sh]
        fo.write(f'>{rid}\n')
        for j in range(0,len(s),80): fo.write(s[j:j+80]+'\n')
        paths=sorted(set(seq_sources[sh])); heads=header_sources[sh]
        fm.write(f'{rid}\t{sh}\t{len(s)}\t{len(paths)}\t{";".join(paths)}\t{" || ".join(heads)}\n')

with open(outdir/'F3B3B_REFERENCE_INVENTORY.tsv','w') as f:
    f.write('path\tsha256\tsource_sequence_n\tretained_reference_sequence_n\tswissprot_focused_by_prior_rvtall_hits\n')
    for r in inventory:
        f.write('\t'.join(str(r[k]) for k in ['path','sha256','source_sequence_n','retained_reference_sequence_n','swissprot_focused_by_prior_rvtall_hits'])+'\n')

summary={
    'resolution_mode':mode,
    'reference_source_file_n':len(inventory),
    'reference_unique_aa_n':len(source_counts),
    'prior_rvtall_swissprot_hit_id_n':len(swiss_hit_ids),
    'sequence_reference_panel_available':len(source_counts)>0,
}
with open(outdir/'F3B3B_REFERENCE_PANEL_SUMMARY.json','w') as f: json.dump(summary,f,indent=2); f.write('\n')
print(json.dumps(summary,indent=2))
