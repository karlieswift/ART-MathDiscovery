#!/usr/bin/env python3
import json, os, sys
from pathlib import Path

src=Path(sys.argv[1])
outdir=Path(sys.argv[2])
nshard=int(sys.argv[3])
expected=int(sys.argv[4])
counts_tsv=Path(sys.argv[5])
marker=Path(sys.argv[6])

outdir.mkdir(parents=True,exist_ok=True)

tmps=[outdir/f".shard_{i:03d}.faa.part" for i in range(nshard)]
finals=[outdir/f"shard_{i:03d}.faa" for i in range(nshard)]

for p in tmps+finals:
    try:p.unlink()
    except FileNotFoundError:pass

outs=[open(p,"wb",buffering=1<<20) for p in tmps]
counts=[0]*nshard
sizes=[0]*nshard

record=-1
out=None
sid=None

with open(src,"rb",buffering=16<<20) as f:
    for line in f:
        if line.startswith(b">"):
            record+=1
            sid=record % nshard
            out=outs[sid]
            counts[sid]+=1

        if out is None:
            raise RuntimeError("sequence before first FASTA header")

        out.write(line)
        sizes[sid]+=len(line)

for f in outs:
    f.close()

n=record+1
if n != expected:
    raise SystemExit(f"FAIL record_n {n} != {expected}")

for a,b in zip(tmps,finals):
    os.replace(a,b)

with open(counts_tsv,"w") as f:
    f.write("shard\trecord_n\tbytes\n")
    for i,(n,b) in enumerate(zip(counts,sizes)):
        f.write(f"shard_{i:03d}\t{n}\t{b}\n")

tmp=marker.with_suffix(".tmp")
tmp.write_text(json.dumps({
    "shard_n":nshard,
    "record_n":n,
    "expected_record_n":expected,
    "counts":counts
},indent=2)+"\n")
os.replace(tmp,marker)

print("record_n =",n)
print("shard_n =",nshard)
print("min_shard_n =",min(counts))
print("max_shard_n =",max(counts))
print("PASS_F3B1_SHARDING")
