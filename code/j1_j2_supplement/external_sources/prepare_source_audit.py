#!/usr/bin/env python3
"""Literal author-coordinate extraction and source audit; never runs a detector.

All coordinates are verified by uniquely matching the author's entire released
GenBank context to the independently retrieved complete native chromosome.
This is secondary structured-RNA source preparation, not an ART-positive cohort.
"""
from pathlib import Path
import argparse
import csv
import datetime
import hashlib
import json
import re


def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()


def sha_file(p):
    return sha_bytes(p.read_bytes())


def write_tsv(p, rows, fields):
    with p.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exclusion-ledger", required=True)
    a = ap.parse_args()
    root = Path(__file__).resolve().parent
    cache = root / "cache"
    incoming = root / "incoming"
    incoming.mkdir(exist_ok=True)
    now = datetime.datetime.now(datetime.timezone.utc).isoformat()
    fa = cache / "NC_012660.1.fna"
    gb = cache / "Minerva_TwoAYGGAY.commit900fb32.gb"
    whole = "".join(s.strip() for s in fa.read_text().splitlines() if not s.startswith(">"))
    whole = whole.upper()
    assert set(whole) == set("ACGT")
    rc = whole.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    author_url = "https://raw.githubusercontent.com/garykbrixi/minerva/900fb32a427f2ebed601569257877167ea45c9ed/minerva/example_data/TwoAYGGAY_Pseudomonas_fluorescens_SBW25.gb"
    sequence_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012660.1&rettype=fasta&retmode=text"
    sample_id = "MINERVA_TWOAYGGAY_SBW25"
    windows, features, units = [], [], {}
    for rec in gb.read_text().split("//"):
        if "ORIGIN\n" not in rec:
            continue
        acc = re.search(r"^ACCESSION +(.+)$", rec, re.M).group(1)
        definition = re.search(r"^DEFINITION +(.+)$", rec, re.M).group(1)
        declared_length = int(re.search(r"^LOCUS\s+\S+\s+(\d+) bp", rec, re.M).group(1))
        seq = "".join(re.findall(r"[a-z]+", rec.split("ORIGIN\n", 1)[1])).upper()
        assert len(seq) == declared_length
        matches = [m.start() for m in re.finditer(re.escape(seq), whole)]
        assert len(matches) == 1, (acc, "author context is not a unique exact native match", matches)
        begin = matches[0]
        central = re.search(r"\|NC_012660\.1\|(\d+)-(\d+)\|([+-])\|ctx_3000$", acc)
        assert central
        assert begin == int(central.group(1)) - 1 - 3000
        windows.append({"author_record": acc, "native_parent_accession": "NC_012660.1", "native_start0": begin,
                        "native_end0": begin + len(seq), "exact_native_matches": 1, "context_sequence_sha256": sha_bytes(seq.encode()),
                        "author_definition": definition, "use_for_detector_input": "NO_TRUTH_CENTERED_CROP"})
        for f in re.finditer(r"^     ([A-Za-z_]+) +([^\n]+)\n((?:                     [^\n]*\n)*)", rec, re.M):
            if f.group(1) != "ncRNA":
                continue
            loc = re.fullmatch(r"(complement\()?([0-9]+)\.\.([0-9]+)\)?", f.group(2))
            assert loc, ("unsupported author feature location", f.group(2))
            local_start0, local_end0 = int(loc.group(2)) - 1, int(loc.group(3))
            assert 0 <= local_start0 < local_end0 <= len(seq)
            strand = "-" if loc.group(1) else "+"
            start0, end0 = begin + local_start0, begin + local_end0
            q = dict(re.findall(r'/([A-Za-z_]+)="([^\"]*)"', f.group(3)))
            assert q.get("rfam_id") == "RF01731"
            dna = whole[start0:end0]
            if strand == "-":
                dna = dna.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            row = {"author_record": acc, "native_parent_accession": "NC_012660.1", "start0": start0, "end0": end0,
                   "strand": strand, "rfam_id": "RF01731", "author_cmsearch_score": q.get("score", ""),
                   "author_cmsearch_evalue": q.get("evalue", ""), "unit_sequence_sha256": sha_bytes(dna.encode()),
                   "source_feature_location": f.group(2), "source_inference": q.get("inference", "")}
            features.append(row)
            key = (start0, end0, strand)
            if key in units:
                assert units[key]["unit_sequence_sha256"] == row["unit_sequence_sha256"]
                units[key]["n_author_context_occurrences"] += 1
            else:
                units[key] = dict(row, n_author_context_occurrences=1, sequence=dna)
    assert len(windows) == 40
    assert len(features) == 56
    assert len(units) == 40
    write_tsv(root / "AUTHOR_CONTEXT_NATIVE_MAPPING.tsv", windows, list(windows[0]))
    write_tsv(root / "AUTHOR_NCRNA_FEATURES_RAW.tsv", features, list(features[0]))
    unit_rows = []
    truth = []
    unit_fasta = []
    for i, key in enumerate(sorted(units), 1):
        u = dict(units[key])
        dna = u.pop("sequence")
        target = f"RF01731_UNIT_{i:03d}"
        u["target_id"] = target
        unit_rows.append(u)
        unit_fasta.append(f">{target}|NC_012660.1|start0={u['start0']}|end0={u['end0']}|strand={u['strand']}\n{dna}\n")
        truth.append({"sample_id": sample_id, "target_id": target, "label_class": "other_structured_rna",
                      "label_scope": "author_defined_rna_interval", "truth_start0": u["start0"], "truth_end0": u["end0"],
                      "array_start0": "", "array_end0": "", "rt_start0": "", "rt_end0": "",
                      "author_truth_url": author_url, "truth_evidence_path": "../cache/" + gb.name,
                      "truth_evidence_sha256": sha_file(gb), "truth_basis": "INDEPENDENT_AUTHOR_ANNOTATION",
                      "curator": "literal_author_GenBank_extraction;author_review_pending",
                      "truth_fixed_utc": now})
    write_tsv(root / "AUTHOR_NCRNA_UNITS_DEDUPLICATED.tsv", unit_rows, list(unit_rows[0]))
    (root / "AUTHOR_NCRNA_UNITS.fna").write_text("".join(unit_fasta))
    write_tsv(incoming / "TRUTH.secondary_structured_rna.tsv", truth, list(truth[0]))
    ledger = Path(a.exclusion_ledger)
    with ledger.open() as h:
        excluded = list(csv.DictReader(h, delimiter="\t"))
    aliases = ["NC_012660.1", "NC_012660", "AM181176.4", "AM181176", "GCF_000009225.2", "GCF_000009225", "PRJNA224116"]
    match_values = set(aliases + [sha_bytes(whole.encode()), sha_bytes(rc.encode())])
    matched = [x for x in excluded if x["value"] in match_values]
    review = {
        "sample_id": sample_id,
        "audit_utc": now,
        "native_accession": "NC_012660.1",
        "accession_aliases": aliases,
        "alias_evidence": "cache/NC_012660.1.record_header.gb; COMMENT states RefSeq identical to AM181176.4",
        "length_bp": len(whole), "sequence_sha256": sha_bytes(whole.encode()), "reverse_complement_sha256": sha_bytes(rc.encode()),
        "exclusion_ledger_filename": ledger.name, "exclusion_ledger_sha256": sha_file(ledger),
        "exclusion_ledger_entries": len(excluded), "exact_accession_alias_or_complete_sequence_matches": matched,
        "prior_exposure_conclusion": "NO_EXACT_LEDGER_MATCH_FOUND" if not matched else "EXPOSED_MATCH_FOUND",
        "strict_never_seen_status": "UNKNOWN",
        "homology_exclusion_status": "NOT_ESTABLISHED_AS_INDEPENDENT_NEW_FAMILY",
        "source_study_overlap": "Minerva paper and author repository had already informed the historical UG27 evaluation. The project audit contains references to the TwoAYGGAY package example, although no scoring record for this native accession was found. Absence of exact accession/hash matches is not evidence of complete non-exposure.",
        "classification": "SECONDARY_ONLY_EXPOSURE_UNCERTAIN",
        "claims_permitted": ["whole-genome frozen-engine response at externally annotated structured-RNA intervals", "one genome with 40 nested target intervals"],
        "claims_not_established": ["new independent ART-positive benchmark", "40 independent genomes", "ART-negative whole genome", "author ncRNA intervals are tandem-repeat array bounds", "new-family held-out generalization"],
        "human_review": "PENDING; no author approval or biological validation is implied",
    }
    (root / "LEAKAGE_REVIEW_NC_012660.1.json").write_text(json.dumps(review, indent=2) + "\n")
    candidate = {
        "sample_id": sample_id, "accession": "NC_012660.1", "native_parent_accession": "NC_012660.1",
        "genome_path": "../cache/NC_012660.1.fna", "fasta_sha256": sha_file(fa), "sequence_sha256": sha_bytes(whole.encode()),
        "taxon": "Pseudomonas [fluorescens] SBW25", "source_url": sequence_url,
        "source_version": "NC_012660.1; author_annotations_commit_900fb32a427f2ebed601569257877167ea45c9ed",
        "source_study": "Li_Brixi_2026_Minerva", "source_publication_date": "2026-09-22",
        "family_cluster": "Minerva_TwoAYGGAY", "leakage_unit": "NC_012660.1_single_native_genome",
        "match_group": "NC_012660.1_single_native_genome", "exposure_status": "EXPOSURE_UNCERTAIN",
        "eligibility_status": "SECONDARY_ONLY",
        "eligibility_reason": "Author-defined RF01731 ncRNA intervals; one complete genome; original ART/UG27 source-study already exposed; not ART-positive or certified negative.",
        "leakage_review_path": "../LEAKAGE_REVIEW_NC_012660.1.json",
        "leakage_review_sha256": sha_file(root / "LEAKAGE_REVIEW_NC_012660.1.json"),
        "leakage_review_status": "REVIEWED_SOURCE_AND_SEQUENCE_GROUPS", "accession_aliases": ";".join(aliases),
    }
    write_tsv(incoming / "CANDIDATES.secondary_structured_rna.tsv", [candidate], list(candidate))
    # Preserve the author table as an unlabeled source pool, never as negative truth.
    cells = json.loads((cache / "sheet1.cells.json").read_text())
    cells = dict(cells)
    pool = []
    for i in range(2, 152):
        pool.append({"author_row": i, "species": cells.get(f"A{i}", ""), "accession_or_source": cells.get(f"B{i}", ""),
                     "phylum": cells.get(f"C{i}", ""), "author_category": cells.get(f"D{i}", ""),
                     "benchmark_label": "UNLABELED_BACKGROUND_SOURCE_ONLY", "locus_truth_status": "NOT_PROVIDED_IN_TABLE_S1",
                     "new_external_eligibility": "PENDING_FULL_EXPOSURE_NATIVE_SOURCE_AND_LABEL_REVIEW"})
    write_tsv(root / "MINERVA_S1_BACKGROUND_SOURCE_POOL.tsv", pool, list(pool[0]))
    coverage = {
        "status": "SOURCE_AUDIT_COMPLETE_PRIMARY_EXTERNAL_ART_DATA_GAP_SECONDARY_NATIVE_DATA_READY",
        "generated_utc": now, "science_scores_read_or_recomputed": False,
        "art_primary_newly_admitted_genomes": 0, "certified_new_negative_genomes": 0,
        "secondary_native_genomes_prepared": 1, "native_nucleotide_bp": len(whole),
        "author_context_records": len(windows), "unique_exact_context_matches": len(windows),
        "raw_author_ncRNA_features": len(features), "unique_author_ncRNA_intervals": len(units),
        "independent_genome_units_in_secondary_data": 1,
        "strict_new_external_status": "NOT_ASSERTED",
        "exact_available_history_matches": len(matched),
        "minerva_author_repository_commit": "900fb32a427f2ebed601569257877167ea45c9ed",
        "minerva_genbank_examples": ["UG27_systems.gb: 5 previously exposed systems", "TwoAYGGAY_Pseudomonas_fluorescens_SBW25.gb: 40 overlapping author contexts from one native genome"],
        "minerva_zenodo_22903368_files": ["Minerva_manuscript_joined.pdf", "Minerva_SupplementaryTables_preprint.xlsx"],
        "minerva_supplementary_tables": {"Table_S1": "150 native input species/source list; no locus coordinates", "Table_S2": "model hyperparameters", "Table_S3": "oligos and experimental plasmids, including previously used five UG27 systems; not additional native array truth"},
        "art_official_pdf": {"sha256": "d43c2fa66aa1d255aba37a0499390bcd86e5b2cafa1f8d249d728e998a4b3b27", "bytes": 21011275, "pages": 40, "embedded_file_count": 0,
                             "data_code_availability_heading_found": False, "repository_or_catalog_download_pointer_found": False,
                             "audit_scope": "Downloaded official Anthropic PDF text plus PDF hyperlinks/embedded files, official announcement, indexed bioRxiv landing page. This does not prove nonavailability elsewhere."},
        "missing_primary_source_materials": ["ART_01-ART_95 to accession/original contig FASTA/version map", "author array presence/absence/not-assessed labels and repeat-copy coordinates", "RT and partner native CDS coordinates and sequences", "NART assessed loci source map and author locus classification", "full unexposed source/near-homology audit after actual source files are available"],
        "readiness": "Secondary manifests are prepared for separate frozen whole-genome evaluation. No detector was run by this source-audit script. Primary ART benchmark remains unready."
    }
    (root / "SOURCE_COVERAGE_AUDIT.json").write_text(json.dumps(coverage, indent=2) + "\n")
    print(json.dumps(coverage, indent=2))


if __name__ == "__main__":
    main()
