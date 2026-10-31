#!/usr/bin/env python3
import hashlib
import json
import sys
from pathlib import Path


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def read_formal_ids(path):
    ids = [x.strip() for x in open(path, encoding='ascii') if x.strip()]
    if len(ids) != 10464:
        raise SystemExit(f'FAIL formal_candidate_n {len(ids)} != 10464')
    if len(set(ids)) != len(ids):
        raise SystemExit('FAIL duplicate formal candidate ids')
    return ids


def iter_fasta(path):
    header = None
    seq = []
    with open(path, 'rb', buffering=16 << 20) as f:
        for line in f:
            if line.startswith(b'>'):
                if header is not None:
                    yield header, b''.join(seq).replace(b'\n', b'').replace(b'\r', b'').upper()
                header = line[1:].strip().decode('ascii')
                seq = []
            elif line.strip():
                seq.append(line.strip())
        if header is not None:
            yield header, b''.join(seq).replace(b'\n', b'').replace(b'\r', b'').upper()


def clean_db_entry(b):
    # MMseqs entries normally end in newline + NUL. Remove only terminators.
    return b.rstrip(b'\x00\r\n')


def load_index_offset_for_keys(index_path, wanted_keys):
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
                out[key] = (int(x[1]), int(x[2]))
                if len(out) == len(wanted_keys):
                    break
    missing = wanted_keys - set(out)
    if missing:
        raise SystemExit(f'FAIL selected DB keys absent from index: n={len(missing)} examples={sorted(missing)[:20]}')
    return out


def fetch_entries(data_path, key_to_offset):
    out = {}
    # Offset-sorted reads reduce random I/O on network storage.
    items = sorted(key_to_offset.items(), key=lambda kv: kv[1][0])
    with open(data_path, 'rb', buffering=0) as f:
        for key, (off, length) in items:
            f.seek(off)
            b = f.read(length)
            if len(b) != length:
                raise SystemExit(f'FAIL short MMseqs DB read key={key} got={len(b)} expected={length}')
            out[key] = clean_db_entry(b)
    return out


def write_fasta_record(f, header, seq):
    f.write(('>' + header + '\n').encode('ascii'))
    for i in range(0, len(seq), 80):
        f.write(seq[i:i+80] + b'\n')


def main():
    if len(sys.argv) != 10:
        raise SystemExit(
            'usage: extract_formal_from_serialized_fasta.py '
            'FORMAL_IDS HMM_INPUT_FASTA ORFDB ORFDB_INDEX ORFDB_H ORFDB_H_INDEX '
            'OUT_NUMERIC_FASTA OUT_HEADER_FASTA OUT_MAP_TSV'
        )

    formal_file, hmm_fa, orfdb, orfidx, orfdbh, orfhdx, out_num, out_hdr, out_map = sys.argv[1:10]
    formal = read_formal_ids(formal_file)
    wanted = set(formal)

    selected = {}
    record_n = 0

    # Critical mapping: HMM query IDs came from this serialized FASTA. Pair each
    # FASTA record position with the MMseqs DB index row at the same position.
    with open(orfidx, encoding='ascii') as idxf:
        for ordinal, (header, seq) in enumerate(iter_fasta(hmm_fa)):
            idxline = idxf.readline()
            if not idxline:
                raise SystemExit(f'FAIL ORFDB.index ended early at serialized record ordinal {ordinal}')
            x = idxline.split()
            if len(x) < 3:
                raise SystemExit(f'FAIL malformed ORFDB.index at ordinal {ordinal}')
            dbkey = x[0]
            qid = header.split()[0]
            if qid in wanted:
                if qid in selected:
                    raise SystemExit(f'FAIL duplicate HMM-input query id: {qid}')
                selected[qid] = {
                    'ordinal': ordinal,
                    'dbkey': dbkey,
                    'seq': seq,
                    'seq_sha256': sha256_bytes(seq),
                }
            record_n += 1

        # Ensure no additional MMseqs sequence entries exist beyond serialized FASTA.
        extra = idxf.readline()
        if extra:
            raise SystemExit('FAIL ORFDB.index has more entries than serialized HMM-input FASTA')

    if record_n != 65315646:
        raise SystemExit(f'FAIL serialized_orf_record_n {record_n} != 65315646')
    missing = wanted - set(selected)
    if missing:
        raise SystemExit(f'FAIL formal IDs absent from exact HMM-input FASTA: n={len(missing)} examples={sorted(missing)[:20]}')
    if len(selected) != 10464:
        raise SystemExit(f'FAIL selected formal records {len(selected)} != 10464')

    keyset = {v['dbkey'] for v in selected.values()}
    if len(keyset) != 10464:
        raise SystemExit(f'FAIL formal query IDs mapped to non-unique ORFDB keys: {len(keyset)}')

    seq_offsets = load_index_offset_for_keys(orfidx, keyset)
    hdr_offsets = load_index_offset_for_keys(orfhdx, keyset)
    db_seqs = fetch_entries(orfdb, seq_offsets)
    db_hdrs = fetch_entries(orfdbh, hdr_offsets)

    # Strong mapping audit: each selected HMM-input peptide must equal the peptide
    # stored at the mapped ORFDB key byte-for-byte after terminator removal.
    mismatch = []
    for qid in formal:
        d = selected[qid]
        dbseq = db_seqs[d['dbkey']].upper()
        if dbseq != d['seq']:
            mismatch.append((qid, d['ordinal'], d['dbkey'], len(d['seq']), len(dbseq)))
            if len(mismatch) >= 20:
                break
    if mismatch:
        raise SystemExit(f'FAIL serialized-record-to-ORFDB-key mapping sequence mismatch: {mismatch}')

    Path(out_num).parent.mkdir(parents=True, exist_ok=True)
    Path(out_hdr).parent.mkdir(parents=True, exist_ok=True)
    Path(out_map).parent.mkdir(parents=True, exist_ok=True)

    with open(out_num, 'wb') as fn, open(out_hdr, 'wb') as fh, open(out_map, 'w', encoding='utf-8') as fm:
        fm.write('candidate_id\tserialized_record_ordinal\torfdb_key\taa_length\tsequence_sha256\tsource_header\n')
        for qid in formal:
            d = selected[qid]
            key = d['dbkey']
            seq = d['seq']
            hdr = db_hdrs[key].decode('ascii', errors='replace').replace('\t', ' ').strip()
            # Keep candidate IDs exactly as used by hmmscan in the numeric FASTA.
            write_fasta_record(fn, qid, seq)
            write_fasta_record(fh, hdr, seq)
            fm.write(
                f"{qid}\t{d['ordinal']}\t{key}\t{len(seq)}\t{d['seq_sha256']}\t{hdr}\n"
            )

    audit = {
        'gate': 'PASS_F3B2_R2_SERIALIZED_QUERYID_TO_ORFDB_KEY_MAPPING',
        'gate_pass': True,
        'serialized_orf_record_n': record_n,
        'formal_candidate_n': len(formal),
        'mapped_unique_orfdb_key_n': len(keyset),
        'sequence_validation_mismatch_n': 0,
        'selection_basis': 'exact hmmscan input FASTA query identifiers',
        'mapping_basis': 'record ordinal paired with original MMseqs ORFDB.index row',
        'scientific_candidate_membership_changed': False,
    }
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
