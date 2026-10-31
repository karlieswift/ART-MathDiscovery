#!/usr/bin/env python3
"""Read-only R5C recovery, validated hit aggregation and guarded replay (stdlib)."""
import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

R5C_STAGE = 'ART_MATH_DISCOVERY_STAGE05D4R5C_CLEAN_RT_CALIBRATION_v1_0_0_20261007'
FORMAT = 'query,target,evalue,bits,qcov,tcov,alnlen,qlen,tlen,fident,alntmscore'
HIT_FIELDS = FORMAT.split(',')
ORDER = {'NO_HIT': 0, 'WEAK': 1, 'SUGGESTIVE': 2, 'STRONG': 3}
POLICY_SHA256 = '71c43302fa3ff859b9575597a3760ebf4bcf0e5bdeb4694f2ac5a0faf8d38902'
ARCHIVE_IDENTITY_SHA256 = '355dd00a319a8ecaba0952823fb892b1ef3ac52bb89698d8c703d5769f739275'
HERE = Path(__file__).resolve().parent


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def fingerprint(path, cache=None):
    """Hash file targets, including symlinks; reject in-flight modification."""
    path = Path(path)
    rec = {'path': str(path), 'is_symlink': path.is_symlink()}
    if path.is_symlink():
        rec['link_text'] = os.readlink(path)
    try:
        rec['resolved_path'] = str(path.resolve(strict=True))
        before = path.stat()
        if not path.is_file():
            rec['state'] = 'NOT_REGULAR_FILE'
            return rec
        rec.update({'bytes': before.st_size, 'mtime_ns': before.st_mtime_ns, 'ctime_ns': before.st_ctime_ns, 'inode': before.st_ino, 'device': before.st_dev})
        cache_key = (str(path.resolve()), before.st_size, before.st_mtime_ns, before.st_ctime_ns, before.st_ino, before.st_dev)
        if cache is not None and cache_key in cache:
            digest = cache[cache_key]
            rec['hash_reused_for_identical_target'] = True
        else:
            digest = sha256(path)
            rec['hash_reused_for_identical_target'] = False
        rec['sha256'] = digest
        after = path.stat()
        stable = (before.st_size, before.st_mtime_ns, before.st_ctime_ns, before.st_ino) == (after.st_size, after.st_mtime_ns, after.st_ctime_ns, after.st_ino)
        rec['state'] = 'HASHED_STABLE' if stable else 'CHANGED_DURING_HASH'
        if stable and cache is not None:
            cache[cache_key] = digest
    except (OSError, RuntimeError) as exc:
        rec.update({'state': 'UNREADABLE_OR_MISSING', 'error': str(exc)})
    return rec


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def write_tsv(path, rows, fields):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def within(path, root):
    try:
        Path(path).resolve().relative_to(Path(root).resolve())
        return True
    except ValueError:
        return False


def forbid_source_output_overlap(output, *sources):
    for source in sources:
        if source and within(output, source):
            raise ValueError('Output would modify the read-only source tree: ' + str(source))


def make_output(project, output, label, protected_sources=()):
    project = Path(project).resolve()
    output = Path(output).resolve() if output else project / 'exports' / ('J1_' + label + '_' + datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%f'))
    if not within(output, project) or not any(within(output, project / name) for name in ('outputs', 'exports')):
        raise ValueError('Output must be under PROJECT/outputs or PROJECT/exports: ' + str(output))
    forbid_source_output_overlap(output, *protected_sources)
    if output.exists() and any(output.iterdir()):
        raise ValueError('Refusing to overwrite nonempty output: ' + str(output))
    output.mkdir(parents=True, exist_ok=True)
    return output


def resolve_stage(project, explicit=None):
    if explicit:
        path = Path(explicit).resolve()
        return ([path] if path.is_dir() else []), 'EXPLICIT_SOURCE'
    project = Path(project)
    candidates = []
    direct = project / R5C_STAGE
    if direct.is_dir():
        candidates.append(direct.resolve())
    # Exact basename only. A tar/unzip wrapper can create one extra level.
    for parent in sorted(project.iterdir()):
        if not parent.is_dir() or parent.name in ('outputs', 'exports', 'tools', '.git'):
            continue
        candidate = parent / R5C_STAGE
        if candidate.is_dir():
            candidates.append(candidate.resolve())
    candidates = sorted(set(candidates))
    # An outer directory consisting solely of the exact same inner stage is packaging.
    candidates = [p for p in candidates if not ((p / R5C_STAGE).is_dir() and len(list(p.iterdir())) == 1)]
    return candidates, 'EXACT_BASENAME_ONLY'


def db_snapshot(prefix, cache=None):
    prefix = Path(prefix)
    files = sorted(prefix.parent.glob(prefix.name + '*')) if prefix.parent.is_dir() else []
    records = []
    for path in files:
        print('[J1 DB hash] ' + str(path), flush=True)
        rec = fingerprint(path, cache)
        rec['component'] = path.name
        records.append(rec)
    stable = [r for r in records if r['state'] == 'HASHED_STABLE']
    names = {r['component'] for r in stable if r['bytes'] > 0}
    families = {}
    for stem in (prefix.name, prefix.name + '_ss', prefix.name + '_ca', prefix.name + '_h'):
        payload = [n for n in names if n == stem or re.fullmatch(re.escape(stem) + r'\.\d+', n)]
        families[stem] = {'index': stem + '.index' in names, 'dbtype': stem + '.dbtype' in names, 'nonempty_payload_components': sorted(payload)}
    ready = bool(records) and len(stable) == len(records) and all(v['index'] and v['dbtype'] and v['nonempty_payload_components'] for v in families.values())
    return {'prefix': str(prefix), 'created_utc': utc(), 'components': records,
            'required_core_families': families,
            'current_prefix_snapshot_state': 'HASHED_WITH_CORE_PAYLOAD_PRESENT' if ready else 'NOT_COMPLETE',
            'historical_full_database_identity': 'NOT_ESTABLISHED_BY_PARTIAL_LEGACY_MANIFEST',
            'boundary': 'All currently visible prefix entries are inventoried, including symlink file targets. A complete historic payload fingerprint or rebuild provenance is still needed to establish original full database identity.'}


def legacy_comparison(snapshot, path):
    current = {r['component']: r for r in snapshot['components']}
    results = []
    for line in Path(path).read_text().splitlines():
        match = re.match(r'^([0-9a-fA-F]{64})\s+\*?(.+)$', line)
        if not match:
            raise ValueError('Invalid sha256 manifest line in ' + str(path))
        expected, original_path = match.groups()
        name = Path(original_path).name
        got = current.get(name)
        state = 'MISSING' if not got or 'sha256' not in got else ('MATCH' if got['sha256'] == expected.lower() else 'MISMATCH')
        results.append({'component': name, 'expected_sha256': expected.lower(), 'observed_sha256': (got or {}).get('sha256', ''), 'state': state})
    return {'manifest': str(path), 'manifest_sha256': sha256(path), 'records': results,
            'all_listed_components_match': bool(results) and all(r['state'] == 'MATCH' for r in results),
            'whole_database_identity_proven': False}


def collect(project, output, explicit=None, max_copy_bytes=64 << 20):
    project, output = Path(project).resolve(), Path(output)
    candidates, discovery = resolve_stage(project, explicit)
    report = {'created_utc': utc(), 'project_root': str(project), 'stage_candidates': [str(p) for p in candidates],
              'discovery_rule': discovery, 'source_mutated': False, 'database_payload_copied': False}
    if len(candidates) != 1:
        report.update({'status': 'NOT_COMPLETE', 'reason': 'STAGE_MISSING' if not candidates else 'AMBIGUOUS_EXACT_STAGE',
                       'scientific_completeness': 'NOT_ASSERTED'})
        write_json(output / 'J1_COLLECTION_REPORT.json', report)
        return report
    stage = candidates[0]
    forbid_source_output_overlap(output, stage)
    records, table_index = [], []
    skipped_dirs = {'tmp', 'temp', 'cache', '.cache', '__pycache__', '.git'}
    for root, dirs, names in os.walk(stage, followlinks=False):
        root = Path(root)
        dirs[:] = sorted(d for d in dirs if d not in skipped_dirs)
        for name in sorted(names):
            path = root / name
            relative = path.relative_to(stage)
            if path.parent == stage / 'references' and path.name.startswith(('clean_rt_afsp', 'matched_nonrt_afsp')):
                continue  # Hashed separately in full, never copied into the return archive.
            rec = {'relative_path': str(relative), 'source_path': str(path), 'is_symlink': path.is_symlink()}
            try:
                rec['bytes'] = path.stat().st_size
                if not path.is_file():
                    rec['state'] = 'NOT_REGULAR_FILE'
                elif rec['bytes'] > max_copy_bytes:
                    rec['state'] = 'DEFERRED_SIZE_LIMIT'
                else:
                    hashed = fingerprint(path)
                    rec.update(hashed)
                    if hashed['state'] == 'HASHED_STABLE':
                        destination = output / 'source' / R5C_STAGE / relative
                        destination.parent.mkdir(parents=True, exist_ok=True)
                        shutil.copyfile(path, destination)
                        rec['copied_sha256'] = sha256(destination)
                        rec['state'] = 'COPIED_VERIFIED' if rec['copied_sha256'] == rec['sha256'] else 'COPY_HASH_MISMATCH'
                        if path.suffix.lower() in ('.tsv', '.csv'):
                            try:
                                with path.open(encoding='utf-8-sig', newline='') as stream:
                                    reader = csv.reader(stream, delimiter='\t' if path.suffix.lower() == '.tsv' else ',')
                                    header = next(reader, [])
                                    n = sum(1 for row in reader if any(row))
                                table_index.append({'relative_path': str(relative), 'data_rows': n, 'header_json': json.dumps(header)})
                            except (UnicodeError, csv.Error) as exc:
                                table_index.append({'relative_path': str(relative), 'data_rows': '', 'header_json': 'UNPARSED: ' + str(exc)})
            except (OSError, RuntimeError) as exc:
                rec.update({'state': 'UNREADABLE_OR_MISSING', 'error': str(exc)})
            records.append(rec)
    cache = {}
    snapshot = db_snapshot(stage / 'references' / 'clean_rt_afsp', cache)
    negative_snapshot = db_snapshot(stage / 'references' / 'matched_nonrt_afsp', cache)
    write_json(output / 'DB_FULL_COMPONENT_SNAPSHOT.json', snapshot)
    write_json(output / 'NONRT_DB_FULL_COMPONENT_SNAPSHOT.json', negative_snapshot)
    write_tsv(output / 'DB_FULL_COMPONENT_SHA256.tsv', snapshot['components'], ['component', 'bytes', 'sha256', 'state', 'path', 'is_symlink', 'link_text', 'resolved_path'])
    write_json(output / 'FILE_COLLECTION_INVENTORY.json', records)
    write_tsv(output / 'TABLE_SCHEMA_INDEX.tsv', table_index, ['relative_path', 'data_rows', 'header_json'])
    policy = stage / 'config' / 'POLICY.json'
    report['policy_historical_fingerprint'] = {'expected_sha256': POLICY_SHA256, 'observed_sha256': sha256(policy) if policy.is_file() else '', 'state': 'MATCH' if policy.is_file() and sha256(policy) == POLICY_SHA256 else 'MISSING_OR_MISMATCH'}
    legacy = HERE / 'reference' / 'R5C_CLEAN_RT_AFSP_SHA256.tsv'
    if legacy.is_file():
        report['legacy_component_comparison'] = legacy_comparison(snapshot, legacy)
    report.update({'status': 'COLLECTED_REVIEW_REQUIRED', 'scientific_completeness': 'NOT_ASSERTED',
                   'copied_files': sum(r['state'] == 'COPIED_VERIFIED' for r in records),
                   'copied_zero_byte_files': sum(r['state'] == 'COPIED_VERIFIED' and r['bytes'] == 0 for r in records),
                   'deferred_or_failed_files': [r for r in records if r['state'] != 'COPIED_VERIFIED'],
                   'database_snapshot_state': snapshot['current_prefix_snapshot_state'],
                   'nonrt_database_snapshot_state': negative_snapshot['current_prefix_snapshot_state'],
                   'unique_database_targets_fully_hashed': len(cache),
                   'unique_database_bytes_read_for_hashing': sum(key[1] for key in cache),
                   'next_step': 'Audit the original control roster and all raw control hits. Current inventory and six matching index/type hashes alone are not a 128-control or full historical database replay.'})
    write_json(output / 'J1_COLLECTION_REPORT.json', report)
    return report


def audit_targets(source_project, manifest, output):
    """Repeat numerical grading of the prior 114-member target audit; do not relabel censored queries."""
    source_project, manifest, output = Path(source_project), Path(manifest), Path(output)
    members = read_tsv(manifest)
    rows, failures = [], []
    for member in members:
        row = {'route': member['route'], 'query_id': member['query_id'], 'raw_path': member['raw_m8_relative_path'],
               'expected_tier': member['recomputed_tier'], 'old_expected_tier': member['old_highest_TM_only_tier']}
        if not member['raw_m8_relative_path']:
            row.update({'state': 'RETAINED_TECHNICAL_UNRESOLVED', 'all_hit_tier': '', 'hit_n': 0})
            rows.append(row)
            continue
        path = source_project / member['raw_m8_relative_path']
        parsed = parse_hits(path)
        valid_mapping = all(h['query'] in (member['query_id'], member['query_id'] + '_model_0') for h in parsed['hits'])
        valid_hash = parsed.get('sha256') == member['raw_m8_sha256']
        cif = source_project / member['cif_relative_path']
        valid_cif = cif.is_file() and sha256(cif) == member['cif_sha256']
        agg = aggregate_hits(parsed['hits'])
        complete = parsed['state'] in ('VALID', 'EMPTY', 'EMPTY_TEXT') and valid_mapping and valid_hash and valid_cif
        matches = agg['max_qualifying_tier'] == member['recomputed_tier'] and len(parsed['hits']) == int(member['raw_hit_count'])
        state = 'MATCH_EXISTING_AUDIT' if complete and matches else 'NOT_COMPLETE_OR_MISMATCH'
        if state != 'MATCH_EXISTING_AUDIT':
            failures.append(member['route'] + ':' + member['query_id'])
        row.update({'state': state, 'raw_file_state': parsed['state'], 'raw_sha256_matches': valid_hash,
                    'cif_sha256_matches': valid_cif, 'explicit_query_mapping_valid': valid_mapping,
                    'all_hit_tier': agg['max_qualifying_tier'], 'highest_tm_only_tier': agg['highest_tm_only_tier'], 'hit_n': len(parsed['hits'])})
        rows.append(row)
    write_tsv(output / 'TARGET_114_RULE_CONSISTENCY.tsv', rows,
              ['route', 'query_id', 'raw_path', 'expected_tier', 'old_expected_tier', 'state', 'raw_file_state',
               'raw_sha256_matches', 'cif_sha256_matches', 'explicit_query_mapping_valid', 'all_hit_tier', 'highest_tm_only_tier', 'hit_n'])
    report = {'status': 'PASS_TARGET_RULE_REAGGREGATION' if not failures else 'NOT_COMPLETE', 'failures': failures,
              'manifest_sha256': sha256(manifest), 'target_members': len(members),
              'raw_results_regraded': sum(r['state'] == 'MATCH_EXISTING_AUDIT' for r in rows),
              'technical_unresolved_retained': sum(r['state'] == 'RETAINED_TECHNICAL_UNRESOLVED' for r in rows),
              'hit_rows_regraded': sum(r['hit_n'] for r in rows),
              'route_tiers': {route: dict(Counter(r['all_hit_tier'] for r in rows if r['route'] == route and r['state'] == 'MATCH_EXISTING_AUDIT')) for route in sorted({r['route'] for r in rows})},
              'new_search_performed': False,
              'scope': 'The same strict parser and all-hit grader as the R5C control reaggregation, applied to hash-identical existing target outputs. Empty raw outputs retain the adjudication of the prior archived audit; this step is not a new search-completion adjudication or structure prediction.'}
    write_json(output / 'TARGET_114_RULE_CONSISTENCY.json', report)
    return report


def hit_tier(tm, qcov, tcov):
    if tm >= 0.50 and qcov >= 0.50 and tcov >= 0.50:
        return 'STRONG'
    if tm >= 0.40 and qcov >= 0.30 and tcov >= 0.30:
        return 'SUGGESTIVE'
    return 'WEAK'


def parse_hits(path, expected_query=None, permit_excluded_self_tm_overshoot=False):
    """Strict 11-column raw Foldseek format; missing, empty and invalid remain distinct."""
    path = Path(path)
    if not path.is_file():
        return {'state': 'MISSING', 'hits': [], 'errors': []}
    digest = sha256(path)
    if path.stat().st_size == 0:
        return {'state': 'EMPTY', 'hits': [], 'errors': [], 'sha256': digest}
    hits, errors, warnings = [], [], []
    with path.open(encoding='utf-8-sig') as stream:
        for index, line in enumerate(stream, 1):
            if not line.strip():
                continue
            values = line.rstrip('\r\n').split('\t')
            if len(values) != 11:
                errors.append({'line': index, 'reason': 'EXPECTED_EXACTLY_11_COLUMNS', 'observed_columns': len(values)})
                continue
            row = dict(zip(HIT_FIELDS, values))
            try:
                if not row['query'] or not row['target']:
                    raise ValueError('empty query/target')
                if expected_query is not None and row['query'] != expected_query:
                    raise ValueError('query ID differs from explicit mapping: ' + row['query'])
                for key in HIT_FIELDS[2:]:
                    row[key] = float(row[key])
                    if not math.isfinite(row[key]):
                        raise ValueError('nonfinite ' + key)
                for key in ('qcov', 'tcov', 'fident'):
                    if not 0 <= row[key] <= 1:
                        raise ValueError('out of range ' + key)
                if not 0 <= row['alntmscore'] <= 1:
                    if permit_excluded_self_tm_overshoot and row['query'] == row['target'] and row['alntmscore'] > 1:
                        warnings.append({'line': index, 'reason': 'EXCLUDED_SELF_HIT_TM_ABOVE_ONE', 'raw_value': row['alntmscore']})
                    else:
                        raise ValueError('out of range alntmscore')
                for key in ('evalue', 'bits'):
                    if row[key] < 0:
                        raise ValueError('negative ' + key)
                for key in ('alnlen', 'qlen', 'tlen'):
                    if row[key] <= 0 or not row[key].is_integer():
                        raise ValueError('invalid positive integer ' + key)
                row['tier'] = hit_tier(row['alntmscore'], row['qcov'], row['tcov'])
                row['source_line'] = index
                hits.append(row)
            except ValueError as exc:
                errors.append({'line': index, 'reason': str(exc)})
    state = 'INVALID' if errors else ('VALID' if hits else 'EMPTY_TEXT')
    return {'state': state, 'hits': hits, 'errors': errors, 'warnings': warnings, 'sha256': digest}


def aggregate_hits(hits):
    if not hits:
        return {'max_qualifying_tier': 'NO_HIT', 'highest_tm_only_tier': 'NO_HIT', 'representative': None}
    # Highest-TM-only is a comparison algorithm, not an assertion about the original R5C implementation.
    old_best = max(hits, key=lambda row: row['alntmscore'])
    best_tier = max((row['tier'] for row in hits), key=ORDER.get)
    pool = [row for row in hits if row['tier'] == best_tier]
    best = min(pool, key=lambda row: (-row['alntmscore'], -row['qcov'], -row['tcov'], row['evalue'], row['target']))
    return {'max_qualifying_tier': best_tier, 'highest_tm_only_tier': old_best['tier'], 'representative': best}


def original_manifest_check(root, manifest):
    root, manifest = Path(root), Path(manifest)
    rows = []
    for line in manifest.read_text().splitlines():
        m = re.fullmatch(r'([0-9a-f]{64})\s+\*?(.+)', line)
        if not m:
            raise ValueError('Invalid archived checksum line: ' + str(manifest))
        expected, rel = m.groups()
        path = root / rel
        if not within(path, root):
            raise ValueError('Archived manifest path escapes source')
        observed = sha256(path) if path.is_file() else ''
        rows.append({'relative_path': rel, 'expected_sha256': expected, 'observed_sha256': observed,
                     'state': 'MATCH' if observed == expected else ('MISSING' if not observed else 'MISMATCH')})
    return rows


def result_db_evidence(prefix, roster):
    """Independently join frozen AF keys to every binary-result record, including NUL-only results."""
    prefix = Path(prefix)
    index_path = Path(str(prefix) + '.index')
    result = {'index_path': str(index_path), 'state': 'NOT_COMPLETE', 'records': {}, 'errors': []}
    if not index_path.is_file():
        result['errors'].append('MISSING_RESULT_INDEX')
        return result
    shards = [p for p in prefix.parent.glob(prefix.name + '.*') if p.suffix[1:].isdigit()]
    shards.sort(key=lambda p: int(p.suffix[1:]))
    if prefix.is_file():
        if shards:
            result['errors'].append('AMBIGUOUS_SINGLE_AND_SHARDED_PAYLOAD')
            return result
        shards = [prefix]
    elif not shards or [int(p.suffix[1:]) for p in shards] != list(range(len(shards))):
        result['errors'].append('MISSING_OR_NONCONTIGUOUS_RESULT_SHARDS')
        return result
    if sum(p.stat().st_size for p in shards) > 1 << 30:
        result['errors'].append('RESULT_PAYLOAD_OVER_1_GIB_REQUIRES_STREAMING_REVIEW')
        return result
    data = b''.join(p.read_bytes() for p in shards)
    result.update({'index_sha256': sha256(index_path), 'payload_bytes': len(data),
                   'shards': [{'name': p.name, 'bytes': p.stat().st_size, 'sha256': sha256(p)} for p in shards]})
    spans = []
    for line_no, line in enumerate(index_path.read_text().splitlines(), 1):
        try:
            key, offset, size = line.split('\t')
            offset, size = int(offset), int(size)
            if key in result['records'] or offset < 0 or size < 1 or offset + size > len(data):
                raise ValueError('duplicate key or invalid offset/size')
            blob = data[offset:offset + size]
            if blob[-1:] != b'\0' or b'\0' in blob[:-1]:
                raise ValueError('invalid NUL-delimited record')
            record_lines = blob[:-1].splitlines()
            target_keys = [s.split(b'\t', 1)[0].decode('ascii') for s in record_lines]
            result['records'][key] = {'index_line': line_no, 'offset': offset, 'bytes': size,
                                      'raw_alignment_rows': len(record_lines), 'target_keys': target_keys,
                                      'result_state': 'COMPLETED_EMPTY' if size == 1 else 'COMPLETED_NONEMPTY',
                                      'record_sha256': hashlib.sha256(blob).hexdigest()}
            spans.append((offset, offset + size))
        except (ValueError, UnicodeError) as exc:
            result['errors'].append('INDEX_LINE_%s: %s' % (line_no, exc))
    expected = {r['afsp_key'] for r in roster}
    if set(result['records']) != expected or len(expected) != len(roster):
        result['errors'].append('INDEX_KEYS_DIFFER_FROM_FROZEN_ROSTER')
    spans.sort()
    if not spans or spans[0][0] != 0 or spans[-1][1] != len(data) or any(a[1] != b[0] for a, b in zip(spans, spans[1:])):
        result['errors'].append('INDEX_OFFSETS_NOT_EXACT_PAYLOAD_PARTITION')
    result['state'] = 'COMPLETE_RESULT_RECORD_COVERAGE' if not result['errors'] else 'NOT_COMPLETE'
    return result


def audit_calibration(source, output):
    """Audit the recovered original R5C schema without inventing a generic roster mapping."""
    source, output = Path(source).resolve(), Path(output)
    forbid_source_output_overlap(output, source)
    files = {
        'positive_roster': 'frozen/STAGE05D4R5C_CLEAN_RT_REFERENCES.tsv',
        'negative_roster': 'frozen/STAGE05D4R5C_MATCHED_NONRT_CONTROLS.tsv',
        'positive_hits': 'raw/POS_VS_CLEAN.m8', 'negative_hits': 'raw/NEG_VS_CLEAN.m8',
        'positive_original_summary': 'outputs/STAGE05D4R5C_KNOWN_RT_CALIBRATION.tsv',
        'negative_original_summary': 'outputs/STAGE05D4R5C_MATCHED_NONRT_CALIBRATION.tsv',
        'original_summary': 'outputs/STAGE05D4R5C_SUMMARY.json',
        'policy': 'config/POLICY.json', 'classifier': 'scripts/summarize_r5c.py',
        'freeze_hashes': 'frozen/STAGE05D4R5C_FREEZE_SHA256.tsv',
        'result_hashes': 'outputs/STAGE05D4R5C_SHA256.tsv'}
    missing = [v for v in files.values() if not (source / v).is_file()]
    report = {'created_utc': utc(), 'source_directory': str(source), 'status': 'NOT_COMPLETE',
              'new_foldseek_search_performed': False, 'new_structure_prediction_performed': False,
              'full_database_replay_complete': False, 'claim_scope': 'Reaggregation of original frozen calibration results, not an independent holdout experiment.'}
    if missing:
        report['missing'] = missing
        write_json(output / 'CALIBRATION_AUDIT_REPORT.json', report)
        return report
    checks = original_manifest_check(source, source / files['freeze_hashes']) + original_manifest_check(source, source / files['result_hashes'])
    write_tsv(output / 'ARCHIVED_SOURCE_CHECKSUM_AUDIT.tsv', checks, ['relative_path', 'expected_sha256', 'observed_sha256', 'state'])
    failures = ['ARCHIVED_CHECKSUM_' + r['state'] + ':' + r['relative_path'] for r in checks if r['state'] != 'MATCH']
    pinned_identity = HERE / 'reference' / 'R5C_ARCHIVE_IDENTITY.json'
    fixed_checks = []
    if not pinned_identity.is_file() or sha256(pinned_identity) != ARCHIVE_IDENTITY_SHA256:
        failures.append('PINNED_ARCHIVE_IDENTITY_MANIFEST_MISSING_OR_CHANGED')
    else:
        for item in json.loads(pinned_identity.read_text())['files']:
            path = source / item['path']
            actual = sha256(path) if path.is_file() else ''
            state = 'MATCH' if actual == item['sha256'] else ('MISSING' if not actual else 'MISMATCH')
            fixed_checks.append({'path': item['path'], 'expected_sha256': item['sha256'], 'observed_sha256': actual, 'state': state})
            if state != 'MATCH':
                failures.append('PINNED_ORIGINAL_ARCHIVE_' + state + ':' + item['path'])
    write_tsv(output / 'PINNED_ORIGINAL_SOURCE_IDENTITY.tsv', fixed_checks, ['path', 'expected_sha256', 'observed_sha256', 'state'])
    policy = json.loads((source / files['policy']).read_text())
    if sha256(source / files['policy']) != POLICY_SHA256:
        failures.append('POLICY_DOES_NOT_MATCH_FROZEN_HISTORICAL_HASH')
    pos, neg = read_tsv(source / files['positive_roster']), read_tsv(source / files['negative_roster'])
    pos_ids, neg_ids = {r['accession'] for r in pos}, {r['accession'] for r in neg}
    if len(pos) != 32 or len(pos_ids) != 32 or len(neg) != 96 or len(neg_ids) != 96 or pos_ids & neg_ids:
        failures.append('CONTROL_MEMBERSHIP_NOT_UNIQUE_32_PLUS_96_DISJOINT')
    expected_target_map = {r['afsp_target']: r['accession'] for r in pos}
    target_key_map = {r['afsp_target']: r['afsp_key'] for r in pos}
    match_counts = Counter(r['matched_rt_accession'] for r in neg)
    if set(match_counts) != pos_ids or any(n != 3 for n in match_counts.values()):
        failures.append('NEGATIVE_MATCHING_NOT_EXACTLY_THREE_PER_RT')
    pos_lengths = {r['accession']: int(r['aa_length']) for r in pos}
    if any(int(r['length_delta']) != int(r['aa_length']) - pos_lengths.get(r['matched_rt_accession'], -1) for r in neg):
        failures.append('LENGTH_MATCHING_DELTA_MISMATCH')
    all_rows, hit_rows, changed = [], [], []
    mode_reports = {}
    for mode, roster in [('POS', pos), ('NEG', neg)]:
        word = 'positive' if mode == 'POS' else 'negative'
        hit_path = source / files[word + '_hits']
        parsed = parse_hits(hit_path, permit_excluded_self_tm_overshoot=mode == 'POS')
        db_evidence = result_db_evidence(source / 'raw' / (mode + '_VS_CLEAN.db'), roster)
        write_json(output / (mode + '_RESULT_RECORD_AUDIT.json'), db_evidence)
        if parsed['state'] not in ('VALID', 'EMPTY', 'EMPTY_TEXT'):
            failures.append(mode + '_RAW_HITS_' + parsed['state'])
        if db_evidence['state'] != 'COMPLETE_RESULT_RECORD_COVERAGE':
            failures.append(mode + '_RESULT_RECORD_COVERAGE_NOT_COMPLETE')
        expected_queries = {r['afsp_target']: r['accession'] for r in roster}
        if len(expected_queries) != len(roster) or len(expected_target_map) != len(pos):
            failures.append(mode + '_NONUNIQUE_FROZEN_AF_TARGET_NAMES')
        by_query = {r['accession']: [] for r in roster}
        for hit in parsed['hits']:
            if hit['query'] not in expected_queries or hit['target'] not in expected_target_map:
                failures.append(mode + '_HIT_MAPPING_OUTSIDE_EXACT_FROZEN_AF_TARGETS_LINE_' + str(hit['source_line']))
                continue
            acc, target_acc = expected_queries[hit['query']], expected_target_map[hit['target']]
            excluded = mode == 'POS' and acc == target_acc
            by_query[acc].append(hit)
            hit_rows.append(dict(hit, cohort=mode, accession=acc, target_accession=target_acc,
                                 excluded_self_hit=excluded, raw_file=files[word + '_hits'], raw_sha256=parsed['sha256']))
        original_rows = read_tsv(source / files[word + '_original_summary'])
        originals = {r['accession']: r for r in original_rows}
        if len(originals) != len(roster) or len(original_rows) != len(roster) or set(originals) != set(by_query):
            failures.append(mode + '_ORIGINAL_SUMMARY_ROSTER_MISMATCH')
        orig_tier_key = 'best_nonself_tier' if mode == 'POS' else 'best_tier'
        orig_tm_key = 'best_nonself_alntmscore' if mode == 'POS' else 'best_alntmscore'
        orig_target_key = 'best_nonself_target' if mode == 'POS' else 'best_rt_target'
        for roster_line, member in enumerate(roster, 2):
            acc = member['accession']
            raw_hits = by_query[acc]
            hits = [h for h in raw_hits if not (mode == 'POS' and expected_target_map[h['target']] == acc)]
            rec = db_evidence['records'].get(member['afsp_key'], {})
            result_matches = rec.get('raw_alignment_rows') == len(raw_hits) and Counter(rec.get('target_keys', [])) == Counter(target_key_map[h['target']] for h in raw_hits)
            state = 'COMPLETE' if result_matches and parsed['state'] in ('VALID', 'EMPTY', 'EMPTY_TEXT') and db_evidence['state'] == 'COMPLETE_RESULT_RECORD_COVERAGE' else 'NOT_COMPLETE'
            if state != 'COMPLETE':
                failures.append(mode + '_QUERY_RESULT_MISMATCH:' + acc)
            agg = aggregate_hits(hits)
            old = max(hits, key=lambda h: h['alntmscore']) if hits else None
            rep = agg['representative']
            original = originals.get(acc, {})
            original_match = original.get(orig_tier_key) == agg['highest_tm_only_tier']
            if old:
                original_match = original_match and original.get(orig_target_key) == expected_target_map[old['target']] and float(original[orig_tm_key]) == old['alntmscore']
                original_match = original_match and float(original['best_nonself_qcov' if mode == 'POS' else 'best_qcov']) == old['qcov'] and float(original['best_nonself_tcov' if mode == 'POS' else 'best_tcov']) == old['tcov']
            original_match = original_match and int(original.get('n_nonself_hits' if mode == 'POS' else 'n_rt_target_hits', -1)) == len(hits)
            if not original_match:
                failures.append(mode + '_ORIGINAL_HIGHEST_TM_SUMMARY_MISMATCH:' + acc)
            row = {'cohort': mode, 'accession': acc, 'afsp_key': member['afsp_key'], 'afsp_target': member['afsp_target'],
                   'aa_length': member['aa_length'], 'annotation_class': member.get('annotation_class', ''),
                   'matched_rt_accession': member.get('matched_rt_accession', ''), 'length_delta': member.get('length_delta', ''),
                   'membership_source': files[word + '_roster'], 'membership_source_line': roster_line,
                   'raw_hits_file': files[word + '_hits'], 'raw_hits_sha256': parsed.get('sha256', ''),
                   'query_result_state': state, 'result_index_key': member['afsp_key'],
                   'result_record_bytes': rec.get('bytes', ''), 'result_record_sha256': rec.get('record_sha256', ''),
                   'raw_hit_count': len(raw_hits), 'self_hits_excluded': len(raw_hits) - len(hits), 'eligible_hit_count': len(hits),
                   'original_reported_tier': original.get(orig_tier_key, ''),
                   'highest_tm_only_tier': agg['highest_tm_only_tier'], 'highest_tm_matches_original_summary': original_match,
                   'max_qualifying_tier': agg['max_qualifying_tier'],
                   'old_best_target': old['target'] if old else '', 'old_best_tm': old['alntmscore'] if old else '',
                   'old_best_qcov': old['qcov'] if old else '', 'old_best_tcov': old['tcov'] if old else '',
                   'old_best_source_line': old['source_line'] if old else '',
                   'representative_target': rep['target'] if rep else '', 'representative_tm': rep['alntmscore'] if rep else '',
                   'representative_qcov': rep['qcov'] if rep else '', 'representative_tcov': rep['tcov'] if rep else '',
                   'representative_evalue': rep['evalue'] if rep else '', 'representative_bits': rep['bits'] if rep else '',
                   'representative_source_line': rep['source_line'] if rep else '',
                   'aggregation_changed': agg['highest_tm_only_tier'] != agg['max_qualifying_tier']}
            all_rows.append(row)
            if row['aggregation_changed']:
                changed.append(row)
        mode_rows = [r for r in all_rows if r['cohort'] == mode]
        mode_reports[mode] = {'members': len(roster), 'raw_hit_rows': len(parsed['hits']), 'raw_query_count': len({r['query'] for r in parsed['hits']}),
                              'excluded_self_hits': sum(r['self_hits_excluded'] for r in mode_rows),
                              'original_highest_tm_tiers': dict(Counter(r['highest_tm_only_tier'] for r in mode_rows)),
                              'all_hit_tiers': dict(Counter(r['max_qualifying_tier'] for r in mode_rows)),
                              'completed_nul_only_result_records': sum(r['bytes'] == 1 for r in db_evidence['records'].values()),
                              'parser_state': parsed['state'], 'parser_errors': parsed['errors'], 'parser_warnings': parsed.get('warnings', [])}
    write_tsv(output / 'CALIBRATION_128_AUDIT.tsv', all_rows, list(all_rows[0]) if all_rows else ['accession'])
    write_tsv(output / 'CALIBRATION_AGGREGATION_CHANGES.tsv', changed, list(all_rows[0]) if all_rows else ['accession'])
    write_tsv(output / 'CALIBRATION_ALL_HITS_AUDIT.tsv', hit_rows,
              ['cohort', 'accession', 'target_accession', 'excluded_self_hit', 'raw_file', 'raw_sha256', 'source_line'] + HIT_FIELDS + ['tier'])
    write_tsv(output / 'NEGATIVE_NO_HIT_EVIDENCE.tsv', [r for r in all_rows if r['cohort'] == 'NEG' and r['max_qualifying_tier'] == 'NO_HIT'], list(all_rows[0]) if all_rows else ['accession'])
    report.update({'status': 'PASS_ARCHIVED_CONTROL_REAGGREGATION' if not failures else 'NOT_COMPLETE',
                   'failures': failures, 'cohorts': mode_reports, 'original_summary': json.loads((source / files['original_summary']).read_text()),
                   'controls_with_changed_aggregation': [r['accession'] for r in changed],
                   'policy_sha256': sha256(source / files['policy']), 'classifier_sha256': sha256(source / files['classifier']),
                   'pinned_original_identity_manifest_sha256': ARCHIVE_IDENTITY_SHA256,
                   'pinned_original_scientific_files_checked': len(fixed_checks),
                   'query_mapping': 'Exact frozen afsp_target string; self exclusion by frozen accession, as in the original source.',
                   'grading': {'strong': 'alntmscore>=0.50 AND qcov>=0.50 AND tcov>=0.50',
                               'suggestive': 'alntmscore>=0.40 AND qcov>=0.30 AND tcov>=0.30',
                               'query_aggregation': 'Maximum qualifying tier across all eligible returned hits',
                               'extra_evalue_or_bits_filter': False, 'coverage_units': 'Raw fractional qcov/tcov; no percent rescaling',
                               'representative': 'Maximum TM within assigned tier; ties by qcov, tcov, E-value, target, as frozen Stage04B9'},
                   'self_hit_tm_note': 'The 32 excluded self hits have raw TM scores 1.001–1.007. They are recorded without clamping and never used in nonself calibration. Nonself out-of-range or nonfinite scores fail validation.',
                   'negative_selection': {'three_per_rt': match_counts == Counter({p: 3 for p in pos_ids}),
                                          'all_96_exact_length_matches': all(int(r['length_delta']) == 0 for r in neg),
                                          'independent_family_holdout': False, 'experimental_nonrt_ground_truth': False},
                   'historical_full_database_identity': 'NOT_ESTABLISHED; archived result replay does not restore search database payloads',
                   'source_files': [{'role': role, 'path': path, 'sha256': sha256(source / path)} for role, path in files.items()]})
    write_json(output / 'CALIBRATION_AUDIT_REPORT.json', report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest='command', required=True)
    col = subs.add_parser('collect', help='Read-only exact-stage collection and full database component hashing.')
    col.add_argument('--project-root', required=True)
    col.add_argument('--r5c-dir')
    col.add_argument('--output')
    col.add_argument('--max-copy-mib', type=int, default=64)
    aud = subs.add_parser('audit-calibration', help='Reaggregate original R5C 32+96 raw hits, with result-record completeness evidence.')
    aud.add_argument('--project-root', required=True)
    aud.add_argument('--source', required=True)
    aud.add_argument('--output')
    targets = subs.add_parser('audit-targets', help='Apply the identical grader to the prior 114-member target audit.')
    targets.add_argument('--project-root', required=True)
    targets.add_argument('--source-project-root', required=True)
    targets.add_argument('--manifest', required=True)
    targets.add_argument('--output')
    plan = subs.add_parser('plan-replay', help='Freeze guarded original DB-vs-DB 128-control commands without running searches.')
    plan.add_argument('--project-root', required=True)
    plan.add_argument('--collection', required=True)
    plan.add_argument('--source', default=str(HERE / 'reference' / 'r5c_archived_results'))
    plan.add_argument('--threads', type=int, default=32)
    plan.add_argument('--output')
    for name in ('run-replay', 'export-db'):
        execution = subs.add_parser(name, help='Execute only the verified fixed protocol in an isolated project export directory.')
        execution.add_argument('--project-root', required=True)
        execution.add_argument('--plan', required=True)
        execution.add_argument('--output')
    portable = subs.add_parser('replay-portable', help='Re-run the exported portable database after moving it to another machine.')
    portable.add_argument('--project-root', required=True)
    portable.add_argument('--portable-root', required=True)
    portable.add_argument('--source', default=str(HERE / 'reference' / 'r5c_archived_results'))
    portable.add_argument('--threads', type=int, default=32)
    portable.add_argument('--output')
    args = parser.parse_args(argv)
    try:
        protected = [getattr(args, name, None) for name in ('source', 'collection', 'portable_root')]
        if args.command == 'collect':
            protected.extend(resolve_stage(Path(args.project_root).resolve(), args.r5c_dir)[0])
        if getattr(args, 'plan', None) and Path(args.plan).is_file():
            plan_data = json.loads(Path(args.plan).read_text())
            protected.extend([plan_data.get('archived_source'), plan_data.get('collection')])
        output = make_output(args.project_root, args.output, args.command, protected)
        if args.command == 'collect':
            report = collect(args.project_root, output, args.r5c_dir, args.max_copy_mib << 20)
        elif args.command == 'audit-calibration':
            report = audit_calibration(args.source, output)
        elif args.command == 'audit-targets':
            report = audit_targets(args.source_project_root, args.manifest, output)
        else:
            import j1_execution
            if args.command == 'plan-replay':
                report = j1_execution.plan_replay(args.collection, args.source, output, args.threads)
            elif args.command == 'run-replay':
                report = j1_execution.run_replay(args.plan, output)
            elif args.command == 'export-db':
                report = j1_execution.export_db(args.plan, output)
            elif args.command == 'replay-portable':
                report = j1_execution.replay_portable(args.portable_root, args.source, output, args.threads)
        print(json.dumps({'output': str(output), 'status': report['status']}, indent=2))
        completed_statuses = {'COLLECTED_REVIEW_REQUIRED', 'READY_FOR_GUARDED_REPLAY',
                              'PASS_ARCHIVED_CONTROL_REAGGREGATION', 'PASS_TARGET_RULE_REAGGREGATION',
                              'PASS_EXACT_SEMANTIC_SEARCH_REPLAY', 'PASS_PORTABLE_DB_AND_EXACT_REPLAY'}
        return 0 if report['status'] in completed_statuses else 2
    except (OSError, ValueError, KeyError) as exc:
        print('J1 ERROR: ' + str(exc), file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
