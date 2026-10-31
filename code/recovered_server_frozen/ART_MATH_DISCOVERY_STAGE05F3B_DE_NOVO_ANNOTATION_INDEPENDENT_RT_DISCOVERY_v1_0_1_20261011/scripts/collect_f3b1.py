#!/usr/bin/env python3
import sys, gzip, json, hashlib
from collections import defaultdict
from pathlib import Path

tbl, dom, lockf, outdir, orf_n = sys.argv[1:6]
outdir = Path(outdir)
outdir.mkdir(parents=True, exist_ok=True)

FORMAL = 1e-4
EXPLORE = 1e-2

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):
            h.update(b)
    return h.hexdigest()

full = {}
with open(tbl) as f:
    for line in f:
        if not line.strip() or line.startswith('#'):
            continue
        x=line.split()
        profile=x[0]
        query=x[2]
        e=float(x[4])
        k=(query,profile)
        if k not in full or e < full[k]:
            full[k]=e

ie = {}
with open(dom) as f:
    for line in f:
        if not line.strip() or line.startswith('#'):
            continue
        x=line.split()
        profile=x[0]
        query=x[3]
        v=float(x[12])
        k=(query,profile)
        if k not in ie or v < ie[k]:
            ie[k]=v

rows=[]
qrows=defaultdict(list)

for k,e in full.items():
    if k not in ie:
        continue
    d=ie[k]

    if e <= FORMAL and d <= FORMAL:
        cls="FORMAL"
    elif e <= EXPLORE and d <= EXPLORE:
        cls="EXPLORATORY"
    else:
        continue

    q,p=k
    rows.append((q,p,e,d,cls))
    qrows[q].append((p,e,d,cls))

def qkey(q):
    try: return (0,int(q))
    except: return (1,q)

rows.sort(key=lambda z:(qkey(z[0]),z[1]))

with gzip.open(
    outdir/"STAGE05F3B1_RVT_PAIR_HITS.tsv.gz",
    "wt"
) as f:
    f.write("query_id\tprofile\tfull_evalue\tbest_domain_i_evalue\tclass\n")
    for r in rows:
        f.write("%s\t%s\t%.8g\t%.8g\t%s\n"%r)

formal=[]
explore=[]

with open(
    outdir/"STAGE05F3B1_CANDIDATE_SUMMARY.tsv",
    "w"
) as f:
    f.write(
        "query_id\tclass\tmodel_n\tbest_full_evalue\t"
        "best_domain_i_evalue\n"
    )

    for q in sorted(qrows,key=qkey):
        r=qrows[q]
        cls="FORMAL" if any(x[3]=="FORMAL" for x in r) else "EXPLORATORY"

        if cls=="FORMAL":
            formal.append(q)
        else:
            explore.append(q)

        f.write(
            f"{q}\t{cls}\t{len(set(x[0] for x in r))}\t"
            f"{min(x[1] for x in r):.8g}\t"
            f"{min(x[2] for x in r):.8g}\n"
        )

(outdir/"STAGE05F3B1_FORMAL_IDS.txt").write_text(
    "".join(x+"\n" for x in formal)
)

(outdir/"STAGE05F3B1_EXPLORATORY_ONLY_IDS.txt").write_text(
    "".join(x+"\n" for x in explore)
)

lock=json.load(open(lockf))

gate={
    "stage":"Stage05F3B1",
    "gate":"PASS_F3B1_FROZEN_ANNOTATION_INDEPENDENT_RVT_SCREEN_COMPLETE",
    "gate_pass":True,
    "source_fasta_sha256":lock["source_fasta_sha256"],
    "protocol_sha256":lock["protocol_sha256"],
    "rvt_hmm_sha256":lock["rvt_hmm_sha256"],
    "rvt_hmm_profile_n":lock["rvt_hmm_profile_n"],
    "orf_peptide_instance_n":int(orf_n),
    "formal_candidate_n":len(formal),
    "exploratory_only_candidate_n":len(explore),
    "qualifying_pair_n":len(rows),
    "formal_full_evalue_threshold":FORMAL,
    "formal_domain_i_evalue_threshold":FORMAL,
    "exploratory_full_evalue_threshold":EXPLORE,
    "exploratory_domain_i_evalue_threshold":EXPLORE,
    "annotation_used_for_candidate_generation":False,
    "claim_boundary":
      "Candidates are annotation-independent RVT-profile-supported RT-like "
      "sequences; biochemical RT/ART identity is not established.",
    "tblout_sha256":sha(tbl),
    "domtblout_sha256":sha(dom)
}

with open(outdir/"STAGE05F3B1_GATE.json","w") as f:
    json.dump(gate,f,indent=2)
    f.write("\n")

print(json.dumps(gate,indent=2))
