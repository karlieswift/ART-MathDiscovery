#!/usr/bin/env python3

import csv
import gzip
import hashlib
import json
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path


def sha256_file(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(16<<20),b""):
            h.update(b)
    return h.hexdigest()


def read_fasta(path):
    h=None
    seq=[]
    with open(path,encoding="ascii") as f:
        for line in f:
            line=line.rstrip("\r\n")
            if line.startswith(">"):
                if h is not None:
                    yield h,"".join(seq).upper()
                h=line[1:]
                seq=[]
            elif line:
                seq.append(line.strip())
        if h is not None:
            yield h,"".join(seq).upper()


def ikey(x):
    try:
        return (0,int(x))
    except:
        return (1,x)


def prepare(args):
    formal_file,numfa,hdrfa,pairgz,outdir=args
    out=Path(outdir)
    out.mkdir(parents=True,exist_ok=True)

    formal=[x.strip() for x in open(formal_file) if x.strip()]
    fset=set(formal)

    assert len(formal)==10464
    assert len(fset)==10464

    num=list(read_fasta(numfa))
    hdr=list(read_fasta(hdrfa))

    assert len(num)==10464, len(num)
    assert len(hdr)==10464, len(hdr)

    records={}
    unknown_acc=0

    for (nh,ns),(hh,hs) in zip(num,hdr):
        q=nh.split()[0]

        assert ns==hs, f"paired FASTA sequence mismatch: {q}"
        assert q in fset, q
        assert q not in records, q

        m=re.search(r"(?:SRR|ERR|DRR)\d+",hh)

        if m:
            acc=m.group(0)
        else:
            acc="UNKNOWN"
            unknown_acc+=1

        records[q]={
            "candidate_id":q,
            "sequence":ns,
            "aa_length":len(ns),
            "source_header":hh.replace("\t"," "),
            "accession":acc
        }

    assert set(records)==fset

    # --------------------------------------------------------
    # Frozen HMM pair evidence
    # --------------------------------------------------------

    pairs=defaultdict(list)

    with gzip.open(pairgz,"rt") as f:
        r=csv.DictReader(f,delimiter="\t")
        for row in r:
            q=row["query_id"]
            if q not in fset:
                continue

            pairs[q].append({
                "profile":row["profile"],
                "full_evalue":float(row["full_evalue"]),
                "domain_ie":float(row["best_domain_i_evalue"]),
                "class":row["class"]
            })

    for q in formal:
        assert q in pairs, f"formal candidate has no HMM pair: {q}"
        assert any(x["class"]=="FORMAL" for x in pairs[q]), q

    # --------------------------------------------------------
    # Exact-AA dedupe
    # --------------------------------------------------------

    seq_groups=defaultdict(list)

    for q,d in records.items():
        seq_groups[d["sequence"]].append(q)

    candidate_to_exact={}
    exact_groups={}

    for seq,members in seq_groups.items():
        members=sorted(members,key=ikey)
        rep=members[0]

        exact_groups[rep]=members

        for q in members:
            candidate_to_exact[q]=rep

    exact_reps=sorted(exact_groups,key=ikey)

    with open(out/"F3B2_EXACT_AA_MAP.tsv","w") as f:
        f.write(
            "candidate_id\texact_rep_id\texact_group_size\t"
            "sequence_sha256\n"
        )

        for rep in exact_reps:
            seq=records[rep]["sequence"]
            sh=hashlib.sha256(seq.encode("ascii")).hexdigest()

            for q in exact_groups[rep]:
                f.write(
                    f"{q}\t{rep}\t{len(exact_groups[rep])}\t{sh}\n"
                )

    with open(out/"F3B2_EXACT_AA_REPRESENTATIVES.faa","w") as f:
        for rep in exact_reps:
            seq=records[rep]["sequence"]
            f.write(f">{rep}\n")
            for i in range(0,len(seq),80):
                f.write(seq[i:i+80]+"\n")

    # --------------------------------------------------------
    # Candidate descriptive table
    # --------------------------------------------------------

    with open(out/"F3B2_CANDIDATE_BASE.tsv","w") as f:
        f.write(
            "candidate_id\taa_length\taccession\texact_rep_id\t"
            "exact_group_size\tformal_model_n\tqualifying_model_n\t"
            "best_full_evalue\tbest_domain_i_evalue\tsource_header\n"
        )

        for q in sorted(formal,key=ikey):
            pp=pairs[q]

            formal_models={
                x["profile"] for x in pp
                if x["class"]=="FORMAL"
            }

            qualifying_models={
                x["profile"] for x in pp
            }

            rep=candidate_to_exact[q]

            f.write(
                "\t".join(map(str,[
                    q,
                    records[q]["aa_length"],
                    records[q]["accession"],
                    rep,
                    len(exact_groups[rep]),
                    len(formal_models),
                    len(qualifying_models),
                    min(x["full_evalue"] for x in pp),
                    min(x["domain_ie"] for x in pp),
                    records[q]["source_header"]
                ]))+"\n"
            )

    summary={
        "formal_candidate_n":len(formal),
        "exact_aa_unique_n":len(exact_reps),
        "exact_duplicate_candidate_n":len(formal)-len(exact_reps),
        "unknown_accession_n":unknown_acc
    }

    with open(out/"F3B2_PREPARE_SUMMARY.json","w") as f:
        json.dump(summary,f,indent=2)
        f.write("\n")

    print(json.dumps(summary,indent=2))
    print("PASS_F3B2_EXACT_AA_COMPRESSION")


def finalize(args):
    base_file,exact_file,family_file,lock_file,outdir=args

    out=Path(outdir)

    lock=json.load(open(lock_file))

    base={}
    with open(base_file) as f:
        for row in csv.DictReader(f,delimiter="\t"):
            base[row["candidate_id"]]=row

    exact={}
    exact_size={}

    with open(exact_file) as f:
        for row in csv.DictReader(f,delimiter="\t"):
            q=row["candidate_id"]
            exact[q]=row["exact_rep_id"]
            exact_size[row["exact_rep_id"]]=int(row["exact_group_size"])

    exact_reps=set(exact.values())

    # MMseqs cluster TSV: representative -> member
    exact_to_family={}

    with open(family_file) as f:
        for line in f:
            if not line.strip():
                continue

            x=line.rstrip("\n").split("\t")

            rep=x[0].split()[0]
            member=x[1].split()[0]

            exact_to_family[member]=rep

    missing=exact_reps-set(exact_to_family)

    if missing:
        raise SystemExit(
            f"family clustering missing exact representatives: {len(missing)}"
        )

    fam_candidates=defaultdict(list)
    fam_exact=defaultdict(set)

    for q in base:
        er=exact[q]
        fr=exact_to_family[er]

        fam_candidates[fr].append(q)
        fam_exact[fr].add(er)

    # --------------------------------------------------------
    # Candidate final map
    # --------------------------------------------------------

    final_file=out/"F3B2_CANDIDATE_FAMILY_MAP.tsv"

    with open(final_file,"w") as f:
        f.write(
            "candidate_id\tfamily_rep_id\texact_rep_id\t"
            "aa_length\taccession\tformal_model_n\t"
            "qualifying_model_n\tbest_full_evalue\t"
            "best_domain_i_evalue\n"
        )

        for q in sorted(base,key=ikey):
            r=base[q]
            er=exact[q]
            fr=exact_to_family[er]

            f.write(
                "\t".join([
                    q,
                    fr,
                    er,
                    r["aa_length"],
                    r["accession"],
                    r["formal_model_n"],
                    r["qualifying_model_n"],
                    r["best_full_evalue"],
                    r["best_domain_i_evalue"]
                ])+"\n"
            )

    # --------------------------------------------------------
    # Family summary
    # --------------------------------------------------------

    family_rows=[]

    for fr,members in fam_candidates.items():
        lengths=[int(base[q]["aa_length"]) for q in members]

        accessions={
            base[q]["accession"]
            for q in members
            if base[q]["accession"]!="UNKNOWN"
        }

        longest=max(
            members,
            key=lambda q:(int(base[q]["aa_length"]),-ikey(q)[1])
            if ikey(q)[0]==0 else (int(base[q]["aa_length"]),0)
        )

        family_rows.append({
            "family_rep_id":fr,
            "candidate_n":len(members),
            "exact_aa_rep_n":len(fam_exact[fr]),
            "accession_n":len(accessions),
            "multirun_replication":len(accessions)>=2,
            "min_aa_length":min(lengths),
            "median_aa_length":statistics.median(lengths),
            "max_aa_length":max(lengths),
            "longest_candidate_id":longest,
            "max_formal_model_n":
                max(int(base[q]["formal_model_n"]) for q in members),
            "max_qualifying_model_n":
                max(int(base[q]["qualifying_model_n"]) for q in members)
        })

    family_rows.sort(
        key=lambda x:(
            -x["accession_n"],
            -x["candidate_n"],
            -x["max_aa_length"],
            ikey(x["family_rep_id"])
        )
    )

    family_file_out=out/"F3B2_FAMILY_SUMMARY.tsv"

    cols=[
        "family_rep_id",
        "candidate_n",
        "exact_aa_rep_n",
        "accession_n",
        "multirun_replication",
        "min_aa_length",
        "median_aa_length",
        "max_aa_length",
        "longest_candidate_id",
        "max_formal_model_n",
        "max_qualifying_model_n"
    ]

    with open(family_file_out,"w") as f:
        f.write("\t".join(cols)+"\n")
        for r in family_rows:
            f.write("\t".join(str(r[c]) for c in cols)+"\n")

    with open(out/"F3B2_FAMILY_TOP100.tsv","w") as f:
        f.write("\t".join(cols)+"\n")
        for r in family_rows[:100]:
            f.write("\t".join(str(r[c]) for c in cols)+"\n")

    gate={
        "stage":"Stage05F3B2",
        "gate":
          "PASS_F3B2_FORMAL_CANDIDATE_FREEZE_AND_FAMILY_COMPRESSION",
        "gate_pass":True,

        "formal_candidate_n":len(base),
        "exact_aa_unique_n":len(exact_reps),
        "sequence_family_n":len(family_rows),

        "multirun_family_n":
          sum(x["multirun_replication"] for x in family_rows),

        "max_accession_replication_n":
          max(x["accession_n"] for x in family_rows),

        "max_family_candidate_n":
          max(x["candidate_n"] for x in family_rows),

        "family_min_seq_id":
          lock["family_min_seq_id"],

        "family_coverage":
          lock["family_coverage"],

        "family_cov_mode":
          lock["family_cov_mode"],

        "family_clustering_role":
          "operational sequence-family compression only; "
          "not proof of biological family identity",

        "known_rt_annotation_used":False,
        "novelty_screen_performed":False,

        "formal_candidate_membership_changed":False,

        "candidate_family_map_sha256":
          sha256_file(final_file),

        "family_summary_sha256":
          sha256_file(family_file_out)
    }

    with open(out/"STAGE05F3B2_GATE.json","w") as f:
        json.dump(gate,f,indent=2)
        f.write("\n")

    print(json.dumps(gate,indent=2))


if __name__=="__main__":
    mode=sys.argv[1]

    if mode=="prepare":
        prepare(sys.argv[2:])
    elif mode=="finalize":
        finalize(sys.argv[2:])
    else:
        raise SystemExit("mode must be prepare or finalize")
