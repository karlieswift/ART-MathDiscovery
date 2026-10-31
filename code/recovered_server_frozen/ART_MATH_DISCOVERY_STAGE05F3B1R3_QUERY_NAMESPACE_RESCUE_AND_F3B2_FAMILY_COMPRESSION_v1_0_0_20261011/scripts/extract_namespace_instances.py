#!/usr/bin/env python3
import csv
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path

EXPECTED_ORF_N = 65315646
EXPECTED_FORMAL_QNAME_N = 10464


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(16 << 20), b''):
            h.update(b)
    return h.hexdigest()


def iter_fasta(path):
    header = None
    seq = []
    with open(path, 'rb', buffering=16 << 20) as f:
        for line in f:
            if line.startswith(b'>'):
                if header is not None:
                    yield header, b''.join(seq).upper()
                header = line[1:].strip().decode('ascii')
                seq = []
            elif line.strip():
                seq.append(line.strip())
        if header is not None:
            yield header, b''.join(seq).upper()


def clean_db_entry(b):
    return b.rstrip(b'\x00\r\n')


def load_index_offsets(index_path, wanted_keys):
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
        raise SystemExit(f'FAIL selected DB keys absent from index n={len(missing)} examples={sorted(missing)[:20]}')
    return out


def fetch_one(fh, offset, length):
    fh.seek(offset)
    b = fh.read(length)
    if len(b) != length:
        raise SystemExit(f'FAIL short MMseqs DB read offset={offset} got={len(b)} expected={length}')
    return clean_db_entry(b)


def write_fasta_record(fh, header, seq):
    fh.write(('>' + header + '\n').encode('ascii'))
    for i in range(0, len(seq), 80):
        fh.write(seq[i:i+80] + b'\n')


def accession_from_header(h):
    m = re.search(r'(?:SRR|ERR|DRR)\d+', h)
    return m.group(0) if m else 'UNKNOWN'


def main():
    if len(sys.argv) != 11:
        raise SystemExit(
            'usage: extract_namespace_instances.py FORMAL_QNAMES HMM_INPUT_FASTA '
            'ORFDB ORFDB_INDEX ORFDB_H ORFDB_H_INDEX OUT_FASTA OUT_MAP OUT_AUDIT OUT_QNAME_COUNTS'
        )
    formal_file, hmm_fa, orfdb, orfidx, orfdbh, orfhdx, out_fa, out_map, out_audit, out_counts = sys.argv[1:11]

    formal = [x.strip() for x in open(formal_file, encoding='ascii') if x.strip()]
    if len(formal) != EXPECTED_FORMAL_QNAME_N or len(set(formal)) != EXPECTED_FORMAL_QNAME_N:
        raise SystemExit(f'FAIL frozen formal query-name set n={len(formal)} unique={len(set(formal))}')
    wanted = set(formal)

    Path(out_fa).parent.mkdir(parents=True, exist_ok=True)
    Path(out_map).parent.mkdir(parents=True, exist_ok=True)

    tmp_fa = Path(str(out_fa) + '.part')
    tmp_map0 = Path(str(out_map) + '.preheader.part')
    for p in (tmp_fa, tmp_map0):
        try: p.unlink()
        except FileNotFoundError: pass

    qcounts = Counter()
    selected_n = 0
    record_n = 0
    selected_keys = set()
    selected_meta = []

    with open(orfidx, encoding='ascii') as idxf, open(tmp_fa, 'wb') as fout:
        for ordinal, (header, seq) in enumerate(iter_fasta(hmm_fa)):
            idxline = idxf.readline()
            if not idxline:
                raise SystemExit(f'FAIL ORFDB.index ended early at ordinal {ordinal}')
            x = idxline.split()
            if len(x) < 3:
                raise SystemExit(f'FAIL malformed ORFDB.index at ordinal {ordinal}')
            dbkey = x[0]
            qname = header.split()[0]
            if qname in wanted:
                qcounts[qname] += 1
                instance_id = f'F3BI{ordinal:08d}'
                seq_sha = hashlib.sha256(seq).hexdigest()
                write_fasta_record(fout, instance_id, seq)
                if dbkey in selected_keys:
                    raise SystemExit(f'FAIL non-unique mapped ORFDB key {dbkey}')
                selected_keys.add(dbkey)
                selected_meta.append((instance_id, qname, ordinal, dbkey, len(seq), seq_sha))
                selected_n += 1
            record_n += 1
        if idxf.readline():
            raise SystemExit('FAIL ORFDB.index has more rows than HMM-input FASTA records')

    if record_n != EXPECTED_ORF_N:
        raise SystemExit(f'FAIL serialized_orf_record_n {record_n} != {EXPECTED_ORF_N}')
    missing = wanted - set(qcounts)
    if missing:
        raise SystemExit(f'FAIL frozen formal query names absent from HMM input n={len(missing)} examples={sorted(missing)[:20]}')
    if selected_n < EXPECTED_FORMAL_QNAME_N:
        raise SystemExit(f'FAIL selected namespace instances {selected_n} < {EXPECTED_FORMAL_QNAME_N}')

    seq_offsets = load_index_offsets(orfidx, selected_keys)
    hdr_offsets = load_index_offsets(orfhdx, selected_keys)

    mismatch = []
    unknown_acc = 0
    with open(orfdb, 'rb', buffering=0) as fs, open(orfdbh, 'rb', buffering=0) as fh, \
         open(out_map, 'w', encoding='utf-8') as fm:
        fm.write('instance_id\toriginal_query_name\tserialized_record_ordinal\torfdb_key\taa_length\tsequence_sha256\taccession\tsource_header\n')
        for instance_id, qname, ordinal, dbkey, aa_len, seq_sha in selected_meta:
            soff, slen = seq_offsets[dbkey]
            dbseq = fetch_one(fs, soff, slen).upper()
            got_sha = hashlib.sha256(dbseq).hexdigest()
            if len(dbseq) != aa_len or got_sha != seq_sha:
                mismatch.append((instance_id, qname, ordinal, dbkey, aa_len, len(dbseq), seq_sha, got_sha))
                if len(mismatch) >= 20:
                    break
            hoff, hlen = hdr_offsets[dbkey]
            hdr = fetch_one(fh, hoff, hlen).decode('ascii', errors='replace').replace('\t', ' ').strip()
            acc = accession_from_header(hdr)
            if acc == 'UNKNOWN':
                unknown_acc += 1
            fm.write(f'{instance_id}\t{qname}\t{ordinal}\t{dbkey}\t{aa_len}\t{seq_sha}\t{acc}\t{hdr}\n')
    if mismatch:
        raise SystemExit(f'FAIL HMM-input ordinal to ORFDB-key sequence validation mismatch: {mismatch}')

    tmp_fa.replace(out_fa)

    with open(out_counts, 'w', encoding='utf-8') as f:
        f.write('original_query_name\thmm_input_instance_n\n')
        for q in sorted(formal, key=lambda z: int(z) if z.isdigit() else z):
            f.write(f'{q}\t{qcounts[q]}\n')

    duplicated = {q:n for q,n in qcounts.items() if n > 1}
    audit = {
        'stage': 'Stage05F3B1R3-namespace-extraction',
        'gate': 'PASS_F3B1R3_FORMAL_QUERY_NAMESPACE_EXPANDED_TO_UNIQUE_INSTANCES',
        'gate_pass': True,
        'serialized_orf_record_n': record_n,
        'frozen_formal_query_name_n': EXPECTED_FORMAL_QNAME_N,
        'selected_namespace_instance_n': selected_n,
        'formal_query_name_with_duplicate_instances_n': len(duplicated),
        'max_query_name_instance_multiplicity': max(qcounts.values()),
        'query_name_instance_multiplicity_sum': sum(qcounts.values()),
        'mapped_unique_orfdb_key_n': len(selected_keys),
        'sequence_validation_mismatch_n': 0,
        'unknown_accession_n': unknown_acc,
        'hmm_input_fasta_sha256': sha256_file(hmm_fa),
        'selection_rule': 'all HMM-input FASTA records whose query-name token belongs to the frozen 10,464 formal query-name set',
        'instance_id_rule': 'F3BI + zero-padded global serialized FASTA record ordinal',
        'scientific_threshold_changed': False,
    }
    with open(out_audit, 'w') as f:
        json.dump(audit, f, indent=2)
        f.write('\n')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
