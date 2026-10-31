#!/usr/bin/env python3
import csv
import hashlib
import json
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(16 << 20), b''):
            h.update(b)
    return h.hexdigest()


def clean_db_entry(b):
    return b.rstrip(b'\x00\r\n')


def load_selected_offsets(index_path, wanted_keys):
    out = {}
    with open(index_path, encoding='ascii') as f:
        for line in f:
            if not line.strip():
                continue
            x = line.split()
            if len(x) < 3:
                raise SystemExit(f'FAIL malformed MMseqs index line: {line[:200]!r}')
            key = x[0]
            if key in wanted_keys:
                if key in out:
                    raise SystemExit(f'FAIL duplicate MMseqs index key: {key}')
                out[key] = (int(x[1]), int(x[2]))
    missing = wanted_keys - set(out)
    if missing:
        raise SystemExit(
            f'FAIL parent contig keys absent from genomeDB_h.index n={len(missing)} '
            f'examples={sorted(missing)[:20]}'
        )
    return out


def fetch_one(fh, offset, length):
    fh.seek(offset)
    b = fh.read(length)
    if len(b) != length:
        raise SystemExit(
            f'FAIL short genomeDB_h read offset={offset} got={len(b)} expected={length}'
        )
    return clean_db_entry(b).decode('ascii', errors='replace').replace('\t', ' ').strip()


def parent_key_from_orf_header(header):
    # MMseqs extractorfs calls Orf::writeOrfHeader(buffer, input_db_key, ...).
    # The input contig key is serialized as the leading integer in the ORF header.
    m = re.match(r'^\s*(\d+)(?:\s|\||\[|$)', header)
    if m:
        return m.group(1)
    # Conservative fallback: only accept a leading decimal token before any text.
    m = re.match(r'^\s*(\d+)', header)
    if m:
        return m.group(1)
    raise ValueError(header)


def accession_from_parent_header(header, allowed):
    token = header.split()[0]
    prefix = token.split('|', 1)[0]
    if prefix in allowed:
        return prefix
    m = re.search(r'(?:SRR|ERR|DRR)\d+', header)
    if m and m.group(0) in allowed:
        return m.group(0)
    return None


def main():
    if len(sys.argv) != 12:
        raise SystemExit(
            'usage: rescue_source_accession.py '
            'R3_INSTANCE_MAP R3_FORMAL_IDS R3_BASE R3_FAMILY_MAP '
            'R3_F3B1_GATE R3_F3B2_GATE GENOMEDB_H GENOMEDB_H_INDEX '
            'F3A_MANIFEST OUTDIR LOCK_JSON'
        )

    (instance_map, formal_ids_file, base_file, family_map_file,
     f3b1_gate_file, f3b2_gate_file, genome_h, genome_h_index,
     manifest_file, outdir, lock_file) = sys.argv[1:12]

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    lock = json.load(open(lock_file))
    g1 = json.load(open(f3b1_gate_file))
    g2 = json.load(open(f3b2_gate_file))

    assert g1['gate_pass'] is True
    assert g2['gate_pass'] is True
    assert g1['corrected_formal_instance_n'] == 10598, g1['corrected_formal_instance_n']
    assert g2['formal_candidate_instance_n'] == 10598, g2['formal_candidate_instance_n']
    assert g2['exact_aa_unique_n'] == 10329, g2['exact_aa_unique_n']
    assert g2['sequence_family_n'] == 4912, g2['sequence_family_n']

    # Canonical 131-source universe.
    with open(manifest_file, encoding='utf-8') as f:
        mr = list(csv.DictReader(f, delimiter='\t'))
    allowed = {r['accession'] for r in mr}
    if len(allowed) != 131:
        raise SystemExit(f'FAIL canonical manifest accession_n={len(allowed)} != 131')

    formal_ids = [x.strip() for x in open(formal_ids_file, encoding='ascii') if x.strip()]
    if len(formal_ids) != 10598 or len(set(formal_ids)) != 10598:
        raise SystemExit(
            f'FAIL corrected formal instance IDs n={len(formal_ids)} unique={len(set(formal_ids))}'
        )
    formal_set = set(formal_ids)

    # Read the R3 namespace map and recover parent contig keys from ORF headers.
    rows = []
    wanted_parent_keys = set()
    parse_fail = []
    with open(instance_map, encoding='utf-8') as f:
        r = csv.DictReader(f, delimiter='\t')
        required = {
            'instance_id', 'original_query_name', 'serialized_record_ordinal',
            'orfdb_key', 'aa_length', 'sequence_sha256', 'accession', 'source_header'
        }
        if not required.issubset(set(r.fieldnames or [])):
            raise SystemExit(f'FAIL R3 instance map missing columns: {sorted(required-set(r.fieldnames or []))}')
        for row in r:
            try:
                parent = parent_key_from_orf_header(row['source_header'])
            except ValueError:
                if len(parse_fail) < 20:
                    parse_fail.append((row['instance_id'], row['source_header']))
                continue
            row['parent_contig_key'] = parent
            rows.append(row)
            wanted_parent_keys.add(parent)

    if parse_fail:
        raise SystemExit(f'FAIL could not parse parent contig key from ORF header examples={parse_fail}')
    if len(rows) != g1['selected_namespace_instance_n']:
        raise SystemExit(
            f"FAIL namespace map row_n={len(rows)} != selected_namespace_instance_n={g1['selected_namespace_instance_n']}"
        )

    # Fetch only needed original contig headers from the canonical genomeDB header DB.
    offsets = load_selected_offsets(genome_h_index, wanted_parent_keys)
    parent_header = {}
    with open(genome_h, 'rb', buffering=0) as fh:
        for key, (off, ln) in offsets.items():
            parent_header[key] = fetch_one(fh, off, ln)

    rescued = {}
    unknown = []
    for row in rows:
        ph = parent_header[row['parent_contig_key']]
        acc = accession_from_parent_header(ph, allowed)
        if acc is None:
            if len(unknown) < 20:
                unknown.append((row['instance_id'], row['parent_contig_key'], ph))
            continue
        row['rescued_accession'] = acc
        row['parent_contig_header'] = ph
        rescued[row['instance_id']] = row

    if unknown:
        raise SystemExit(f'FAIL parent contig headers could not be assigned to canonical 131 accessions examples={unknown}')
    if len(rescued) != len(rows):
        raise SystemExit(f'FAIL rescued instance n={len(rescued)} != namespace row n={len(rows)}')

    formal_missing = formal_set - set(rescued)
    if formal_missing:
        raise SystemExit(f'FAIL formal instances absent from rescued map n={len(formal_missing)} examples={sorted(formal_missing)[:20]}')

    # Write complete namespace provenance rescue.
    rescued_map = out / 'F3B1R3_INSTANCE_MAP_WITH_CANONICAL_SOURCE.tsv'
    fields = [
        'instance_id', 'original_query_name', 'serialized_record_ordinal',
        'orfdb_key', 'parent_contig_key', 'aa_length', 'sequence_sha256',
        'rescued_accession', 'source_header', 'parent_contig_header'
    ]
    with open(rescued_map, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for row in rows:
            w.writerow({k: row[k] for k in fields})

    # Load R3 candidate base table; patch accession only.
    base = {}
    with open(base_file, encoding='utf-8') as f:
        r = csv.DictReader(f, delimiter='\t')
        base_fields = list(r.fieldnames or [])
        for row in r:
            base[row['candidate_id']] = row
    if set(base) != formal_set:
        raise SystemExit(
            f'FAIL R3 base candidate set differs from corrected formal set: '
            f'base={len(base)} formal={len(formal_set)}'
        )

    patched_base = out / 'F3B2_R4_CANDIDATE_BASE_WITH_SOURCE.tsv'
    with open(patched_base, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=base_fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for q in sorted(base, key=lambda q: int(base[q]['serialized_record_ordinal'])):
            row = dict(base[q])
            row['accession'] = rescued[q]['rescued_accession']
            w.writerow(row)

    # Patch existing family membership table without reclustering.
    fmap = {}
    with open(family_map_file, encoding='utf-8') as f:
        r = csv.DictReader(f, delimiter='\t')
        fmap_fields = list(r.fieldnames or [])
        for row in r:
            fmap[row['candidate_id']] = row
    if set(fmap) != formal_set:
        raise SystemExit(
            f'FAIL R3 family-map candidate set differs from corrected formal set: '
            f'family_map={len(fmap)} formal={len(formal_set)}'
        )

    patched_fmap = out / 'F3B2_R4_CANDIDATE_FAMILY_MAP_WITH_SOURCE.tsv'
    patched_rows = {}
    with open(patched_fmap, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fmap_fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for q in sorted(fmap, key=lambda q: int(base[q]['serialized_record_ordinal'])):
            row = dict(fmap[q])
            row['accession'] = rescued[q]['rescued_accession']
            patched_rows[q] = row
            w.writerow(row)

    # Recompute family descriptive metadata only; family membership is frozen from R3.
    fam = defaultdict(list)
    for q, row in patched_rows.items():
        fam[row['family_rep_id']].append(q)

    if len(fam) != 4912:
        raise SystemExit(f'FAIL family_n={len(fam)} != frozen R3 family_n=4912')

    family_rows = []
    for fr, members in fam.items():
        lengths = [int(patched_rows[q]['aa_length']) for q in members]
        accessions = {patched_rows[q]['accession'] for q in members}
        qnames = {patched_rows[q]['original_query_name'] for q in members}
        exact_reps = {patched_rows[q]['exact_rep_id'] for q in members}
        longest = max(
            members,
            key=lambda q: (
                int(patched_rows[q]['aa_length']),
                -int(base[q]['serialized_record_ordinal'])
            )
        )
        family_rows.append({
            'family_rep_id': fr,
            'candidate_instance_n': len(members),
            'original_query_name_n': len(qnames),
            'exact_aa_rep_n': len(exact_reps),
            'accession_n': len(accessions),
            'multirun_replication': len(accessions) >= 2,
            'min_aa_length': min(lengths),
            'median_aa_length': statistics.median(lengths),
            'max_aa_length': max(lengths),
            'longest_candidate_id': longest,
            'max_formal_model_n': max(int(patched_rows[q]['formal_model_n']) for q in members),
            'max_qualifying_model_n': max(int(patched_rows[q]['qualifying_model_n']) for q in members),
        })

    family_rows.sort(
        key=lambda x: (
            -x['accession_n'],
            -x['candidate_instance_n'],
            -x['max_aa_length'],
            x['family_rep_id']
        )
    )

    cols = [
        'family_rep_id', 'candidate_instance_n', 'original_query_name_n',
        'exact_aa_rep_n', 'accession_n', 'multirun_replication',
        'min_aa_length', 'median_aa_length', 'max_aa_length',
        'longest_candidate_id', 'max_formal_model_n', 'max_qualifying_model_n'
    ]

    family_summary = out / 'F3B2_R4_FAMILY_SUMMARY_WITH_SOURCE.tsv'
    family_top100 = out / 'F3B2_R4_FAMILY_TOP100_WITH_SOURCE.tsv'
    for path, subset in ((family_summary, family_rows), (family_top100, family_rows[:100])):
        with open(path, 'w', encoding='utf-8', newline='') as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', lineterminator='\n')
            w.writeheader()
            w.writerows(subset)

    formal_accessions = {rescued[q]['rescued_accession'] for q in formal_set}
    namespace_accessions = {row['rescued_accession'] for row in rescued.values()}

    gate = {
        'stage': 'Stage05F3B2-R4',
        'gate': 'PASS_F3B2_R4_CANONICAL_SOURCE_ACCESSION_METADATA_RESCUE',
        'gate_pass': True,
        'rescue_role': 'metadata/provenance repair only; HMM classification, exact-AA grouping, and sequence-family membership are unchanged',
        'formal_candidate_instance_n': len(formal_set),
        'exact_aa_unique_n': g2['exact_aa_unique_n'],
        'sequence_family_n': len(family_rows),
        'selected_namespace_instance_n': len(rows),
        'selected_namespace_unknown_accession_n': 0,
        'formal_unknown_accession_n': 0,
        'canonical_source_universe_accession_n': len(allowed),
        'namespace_accession_n': len(namespace_accessions),
        'formal_accession_n': len(formal_accessions),
        'multirun_family_n': sum(bool(r['multirun_replication']) for r in family_rows),
        'max_accession_replication_n': max(r['accession_n'] for r in family_rows),
        'max_family_candidate_instance_n': max(r['candidate_instance_n'] for r in family_rows),
        'family_membership_changed': False,
        'exact_aa_membership_changed': False,
        'hmm_results_changed': False,
        'scientific_threshold_changed': False,
        'known_rt_annotation_used': False,
        'novelty_screen_performed': False,
        'accession_source': 'original canonical F3A contig header in F3B1 v1.0.0 genomeDB_h, reached via parent contig key encoded by MMseqs extractorfs in each ORF header',
        'family_clustering_role': g2['family_clustering_role'],
        'f3b1r3_gate_sha256': sha256_file(f3b1_gate_file),
        'f3b2_r3_gate_sha256': sha256_file(f3b2_gate_file),
        'r3_instance_map_sha256': sha256_file(instance_map),
        'r3_family_map_sha256': sha256_file(family_map_file),
        'canonical_manifest_sha256': sha256_file(manifest_file),
        'rescued_instance_map_sha256': sha256_file(rescued_map),
        'rescued_candidate_family_map_sha256': sha256_file(patched_fmap),
        'rescued_family_summary_sha256': sha256_file(family_summary),
        'implementation_lock_sha256': sha256_file(lock_file),
        'claim_boundary': 'R4 repairs source-run/accession metadata only. It does not alter candidate identity, HMM support, exact-AA identity, operational family membership, or novelty status.'
    }

    with open(out / 'STAGE05F3B2_R4_GATE.json', 'w') as f:
        json.dump(gate, f, indent=2)
        f.write('\n')

    print(json.dumps(gate, indent=2))
    print('PASS_STAGE05F3B2_R4_SOURCE_ACCESSION_METADATA_RESCUE_COMPLETE')


if __name__ == '__main__':
    main()
