"""Execution backend: fixed original R5C protocol; no database downloads or rebuilding."""
import json
import os
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path

import j1_reproducibility as core

RUNNER_SHA = '554d6ed8cc3d2a8477c097e79eadf1c5659ad26023ae9d7b95e7c2b8b996a782'
VERSIONS = {'10.941cd33', '941cd33ff0771cd2e3f144e3293e22a2b87e9fda'}


def unchanged(snapshot):
    errors = []
    prefix = Path(snapshot['prefix'])
    current_names = {p.name for p in prefix.parent.glob(prefix.name + '*')}
    recorded_names = {item['component'] for item in snapshot.get('components', [])}
    if current_names != recorded_names:
        errors.append('DATABASE_COMPONENT_SET_CHANGED:' + str(prefix))
    for item in snapshot.get('components', []):
        path = Path(item['path'])
        try:
            stat = path.stat()
            if item.get('state') != 'HASHED_STABLE' or not path.is_file() or str(path.resolve()) != item['resolved_path'] or (stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns, stat.st_ino, stat.st_dev) != (item['bytes'], item['mtime_ns'], item.get('ctime_ns'), item['inode'], item['device']):
                errors.append('INPUT_CHANGED:' + str(path))
        except OSError:
            errors.append('INPUT_UNAVAILABLE:' + str(path))
    return errors


def index_keys(path):
    keys = []
    for line in Path(path).read_text().splitlines():
        values = line.split('\t')
        if len(values) != 3 or not values[0].isdigit() or not values[1].isdigit() or not values[2].isdigit():
            raise ValueError('Malformed database index: ' + str(path))
        keys.append(values[0])
    if len(keys) != len(set(keys)):
        raise ValueError('Duplicate keys in ' + str(path))
    return set(keys)


def commands(executable, positive_db, negative_db, threads, run='{RUN}'):
    output = str(run)
    result = []
    for mode, query in [('POS', positive_db), ('NEG', negative_db)]:
        result_db = output + '/raw/' + mode + '_VS_CLEAN.db'
        result.append({'name': mode.lower() + '_search', 'argv': [executable, 'search', query, positive_db, result_db,
                      output + '/work/' + mode.lower() + '_search', '-e', '1.0', '--max-seqs', '100',
                      '--threads', str(threads), '--exact-tmscore', '1', '-a']})
        result.append({'name': mode.lower() + '_convertalis', 'argv': [executable, 'convertalis', query, positive_db,
                      result_db, output + '/raw/' + mode + '_VS_CLEAN.m8', '--format-output', core.FORMAT, '--exact-tmscore', '1']})
    return result


def runtime_info():
    executable = shutil.which('foldseek')
    info = {'path': executable or '', 'state': 'MISSING'}
    if executable:
        proc = subprocess.run([executable, 'version'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=30)
        version = proc.stdout.strip()
        info.update({'version': version, 'sha256': core.sha256(executable),
                     'state': 'ORIGINAL_VERSION' if proc.returncode == 0 and version in VERSIONS else 'DIFFERENT_OR_UNRECOGNIZED_VERSION'})
    return info


def plan_replay(collection, source, output, threads=32):
    collection, source, output = Path(collection).resolve(), Path(source).resolve(), Path(output)
    core.forbid_source_output_overlap(output, source, collection)
    report = {'created_utc': core.utc(), 'status': 'NOT_COMPLETE', 'threads': threads,
              'archived_source': str(source), 'collection': str(collection), 'blocking_reasons': []}
    reasons = report['blocking_reasons']
    if threads < 1:
        reasons.append('THREADS_MUST_BE_POSITIVE')
    required = [collection / 'J1_COLLECTION_REPORT.json', collection / 'DB_FULL_COMPONENT_SNAPSHOT.json', collection / 'NONRT_DB_FULL_COMPONENT_SNAPSHOT.json']
    if any(not p.is_file() for p in required):
        reasons.append('COLLECTION_OR_BOTH_DATABASE_SNAPSHOTS_MISSING')
        core.write_json(output / 'REPLAY_PLAN.json', report)
        return report
    collection_report, pos_snapshot, neg_snapshot = [json.loads(p.read_text()) for p in required]
    report['database_snapshots'] = {'POS': pos_snapshot, 'NEG': neg_snapshot}
    for mode, snapshot in report['database_snapshots'].items():
        if snapshot['current_prefix_snapshot_state'] != 'HASHED_WITH_CORE_PAYLOAD_PRESENT':
            reasons.append(mode + '_DATABASE_PAYLOAD_NOT_COMPLETE')
        reasons.extend(unchanged(snapshot))
    legacy = collection_report.get('legacy_component_comparison', {})
    if not legacy.get('all_listed_components_match'):
        reasons.append('ORIGINAL_SIX_CLEAN_REFERENCE_COMPONENTS_NOT_MATCHED')
    if collection_report.get('policy_historical_fingerprint', {}).get('state') != 'MATCH':
        reasons.append('SERVER_POLICY_NOT_MATCHED_TO_FROZEN_HASH')
    source_audit = core.audit_calibration(source, output / 'archived_source_audit')
    if source_audit['status'] != 'PASS_ARCHIVED_CONTROL_REAGGREGATION':
        reasons.append('ORIGINAL_CALIBRATION_SOURCE_AUDIT_INCOMPLETE')
    report['original_expected_tiers'] = {mode: d['all_hit_tiers'] for mode, d in source_audit.get('cohorts', {}).items()}
    report['source_files'] = source_audit.get('source_files', [])
    for mode, name in [('POS', 'CLEAN_RT_REFERENCES'), ('NEG', 'MATCHED_NONRT_CONTROLS')]:
        roster_path = source / 'frozen' / ('STAGE05D4R5C_' + name + '.tsv')
        if roster_path.is_file():
            roster = core.read_tsv(roster_path)
            prefix = report['database_snapshots'][mode]['prefix']
            try:
                if index_keys(prefix + '.index') != {r['afsp_key'] for r in roster}:
                    reasons.append(mode + '_DB_PRIMARY_KEYS_DIFFER_FROM_ORIGINAL_ROSTER')
            except (OSError, ValueError) as exc:
                reasons.append(str(exc))
    protocol = core.HERE / 'reference' / 'original_r5c_protocol' / 'run_stage05d4r5c.sh'
    report['original_runner'] = {'path': str(protocol), 'expected_sha256': RUNNER_SHA,
                                 'observed_sha256': core.sha256(protocol) if protocol.is_file() else ''}
    if report['original_runner']['observed_sha256'] != RUNNER_SHA:
        reasons.append('ORIGINAL_PROTOCOL_SOURCE_HASH_MISMATCH')
    report['runtime'] = runtime_info()
    if report['runtime']['state'] != 'ORIGINAL_VERSION':
        reasons.append('FOLDSEEK_10_941cd33_REQUIRED_IN_CURRENT_ENVIRONMENT')
    report['commands'] = commands(report['runtime'].get('path') or 'foldseek', pos_snapshot['prefix'], neg_snapshot['prefix'], threads)
    report['status'] = 'READY_FOR_GUARDED_REPLAY' if not reasons else 'NOT_COMPLETE'
    report['scientific_parameters_changed'] = False
    report['database_guard'] = 'Full SHA-256 once per unique resolved file during collection; component-set, resolved target, size, mtime, ctime, inode and device are checked at planning, execution preflight and after execution. No repeated full scan of shared AlphaFold payloads.'
    report['boundary'] = 'This plan replays the original 32 known RT and 96 matched controls using the original DB-vs-DB search and convertalis protocol. Historical full payload identity remains distinct from a newly fingerprinted current database and an output-matching replay.'
    core.write_json(output / 'REPLAY_PLAN.json', report)
    return report


def validate_plan(plan):
    reasons = []
    if plan.get('status') != 'READY_FOR_GUARDED_REPLAY':
        return ['PLAN_NOT_READY'] + plan.get('blocking_reasons', [])
    for snapshot in plan['database_snapshots'].values():
        reasons.extend(unchanged(snapshot))
    runtime = runtime_info()
    if runtime != plan.get('runtime'):
        reasons.append('FOLDSEEK_EXECUTABLE_OR_VERSION_CHANGED_SINCE_PLAN')
    protocol = Path(plan['original_runner']['path'])
    if not protocol.is_file() or core.sha256(protocol) != RUNNER_SHA:
        reasons.append('ORIGINAL_PROTOCOL_SOURCE_CHANGED')
    for item in plan['source_files']:
        path = Path(plan['archived_source']) / item['path']
        if not path.is_file() or core.sha256(path) != item['sha256']:
            reasons.append('ARCHIVED_SOURCE_CHANGED:' + item['path'])
    expected = commands(plan['runtime']['path'], plan['database_snapshots']['POS']['prefix'], plan['database_snapshots']['NEG']['prefix'], plan['threads'])
    if plan['commands'] != expected:
        reasons.append('COMMAND_PLAN_DIFFERS_FROM_FIXED_ORIGINAL_PROTOCOL')
    return reasons


def execute(argv, log_path):
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    print('[J1 execute] ' + json.dumps(argv), flush=True)
    started = core.utc()
    with log_path.open('w') as stream:
        proc = subprocess.run(argv, stdout=stream, stderr=subprocess.STDOUT)
    return {'argv': argv, 'exit_code': proc.returncode, 'started_utc': started, 'finished_utc': core.utc(),
            'log_path': str(log_path), 'log_sha256': core.sha256(log_path)}


def compare_replay(source, replay, output):
    source, replay, output = Path(source), Path(replay), Path(output)
    comparisons, errors, details = [], [], {}
    for mode, name in [('POS', 'CLEAN_RT_REFERENCES'), ('NEG', 'MATCHED_NONRT_CONTROLS')]:
        roster = core.read_tsv(source / 'frozen' / ('STAGE05D4R5C_' + name + '.tsv'))
        original = core.parse_hits(source / 'raw' / (mode + '_VS_CLEAN.m8'), permit_excluded_self_tm_overshoot=mode == 'POS')
        new = core.parse_hits(replay / 'raw' / (mode + '_VS_CLEAN.m8'), permit_excluded_self_tm_overshoot=mode == 'POS')
        evidence = core.result_db_evidence(replay / 'raw' / (mode + '_VS_CLEAN.db'), roster)
        core.write_json(output / (mode + '_REPLAY_RESULT_RECORD_AUDIT.json'), evidence)
        if original['state'] != 'VALID' or new['state'] not in ('VALID', 'EMPTY', 'EMPTY_TEXT') or evidence['state'] != 'COMPLETE_RESULT_RECORD_COVERAGE':
            errors.append(mode + '_INVALID_OR_INCOMPLETE_SEARCH_RESULT')
        expected_names = {r['afsp_target'] for r in roster}
        target_map = {r['afsp_target']: r['afsp_key'] for r in core.read_tsv(source / 'frozen' / 'STAGE05D4R5C_CLEAN_RT_REFERENCES.tsv')}
        if any(h['query'] not in expected_names or h['target'] not in target_map for h in new['hits']):
            errors.append(mode + '_UNMAPPED_QUERY_OR_TARGET')
        for member in roster:
            name = member['afsp_target']
            old = [h for h in original['hits'] if h['query'] == name and not (mode == 'POS' and h['target'] == name)]
            raw_new = [h for h in new['hits'] if h['query'] == name]
            new_hits = [h for h in raw_new if not (mode == 'POS' and h['target'] == name)]
            old_agg, new_agg = core.aggregate_hits(old), core.aggregate_hits(new_hits)
            old_sem = Counter(tuple(h[key] for key in core.HIT_FIELDS) for h in old)
            new_sem = Counter(tuple(h[key] for key in core.HIT_FIELDS) for h in new_hits)
            rec = evidence['records'].get(member['afsp_key'], {})
            complete = rec.get('raw_alignment_rows') == len(raw_new) and Counter(rec.get('target_keys', [])) == Counter(target_map.get(h['target'], 'UNMAPPED') for h in raw_new)
            row = {'cohort': mode, 'accession': member['accession'], 'original_all_hit_tier': old_agg['max_qualifying_tier'],
                   'replayed_all_hit_tier': new_agg['max_qualifying_tier'], 'original_eligible_hits': len(old),
                   'replayed_eligible_hits': len(new_hits), 'exact_semantic_hit_multiset_match': old_sem == new_sem,
                   'result_record_coverage_complete': complete,
                   'tier_matches': old_agg['max_qualifying_tier'] == new_agg['max_qualifying_tier']}
            comparisons.append(row)
            if not complete or not row['tier_matches']:
                errors.append(mode + '_QUERY_REPLAY_MISMATCH:' + member['accession'])
        details[mode] = {'raw_output_sha256': new.get('sha256'), 'original_raw_output_sha256': original.get('sha256'),
                         'raw_file_byte_identical': new.get('sha256') == original.get('sha256'), 'parse_errors': new['errors']}
    core.write_tsv(output / 'CONTROL_128_SEARCH_REPLAY_COMPARISON.tsv', comparisons, list(comparisons[0]))
    exact = all(r['exact_semantic_hit_multiset_match'] for r in comparisons)
    return {'status': 'PASS_EXACT_SEMANTIC_SEARCH_REPLAY' if not errors and exact else ('TIERS_MATCH_HIT_VALUES_DIFFER_REVIEW_REQUIRED' if not errors else 'NOT_COMPLETE'),
            'failures': errors, 'all_128_tiers_match': not errors, 'all_eligible_hit_values_match': exact,
            'queries': len(comparisons), 'raw_outputs': details}


def run_searches(plan, output, prefixes=None):
    output = Path(output)
    for name in ('raw', 'work', 'logs'):
        (output / name).mkdir(parents=True, exist_ok=True)
    prefixes = prefixes or {mode: plan['database_snapshots'][mode]['prefix'] for mode in ('POS', 'NEG')}
    receipts = []
    for command in commands(plan['runtime']['path'], prefixes['POS'], prefixes['NEG'], plan['threads'], output):
        receipt = execute(command['argv'], output / 'logs' / (command['name'] + '.log'))
        receipts.append(receipt)
        core.write_json(output / 'COMMAND_RECEIPTS.json', receipts)
        if receipt['exit_code'] != 0:
            return {'status': 'NOT_COMPLETE', 'failures': ['COMMAND_FAILED:' + command['name']], 'command_receipts': receipts}
    report = compare_replay(plan['archived_source'], output, output)
    report['command_receipts'] = receipts
    return report


def run_replay(plan_path, output):
    plan_path, output = Path(plan_path), Path(output)
    plan = json.loads(plan_path.read_text())
    core.forbid_source_output_overlap(output, plan.get('archived_source'), plan.get('collection'))
    errors = validate_plan(plan)
    if errors:
        report = {'status': 'NOT_COMPLETE', 'blocking_reasons': errors, 'search_performed': False}
    else:
        report = run_searches(plan, output)
        after = [e for snapshot in plan['database_snapshots'].values() for e in unchanged(snapshot)]
        if after:
            report.update({'status': 'NOT_COMPLETE', 'input_changed_after_search': after})
        report['search_performed'] = True
    report['plan_sha256'] = core.sha256(plan_path)
    report['source_inputs_modified'] = False
    core.write_json(output / 'SEARCH_REPLAY_REPORT.json', report)
    return report


def fasta(path):
    result, name, sequence = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if name is not None:
                if name in result:
                    raise ValueError('Duplicate FASTA header: ' + name)
                result[name] = ''.join(sequence)
            name, sequence = line[1:].split()[0], []
        elif line.strip():
            sequence.append(line.strip())
    if name is not None:
        if name in result:
            raise ValueError('Duplicate FASTA header: ' + name)
        result[name] = ''.join(sequence)
    return result


def materialize_core_database(executable, source_prefix, key_file, target_prefix, work, log_dir):
    """Official base:createsubdb per component; only regular core outputs are assembled."""
    source_prefix, target_prefix, work, log_dir = str(source_prefix), Path(target_prefix), Path(work), Path(log_dir)
    target_prefix.parent.mkdir(parents=True, exist_ok=True)
    keyset = set(Path(key_file).read_text().splitlines())
    receipts = []
    for suffix, label in [('', 'aa'), ('_ss', '3di'), ('_ca', 'ca'), ('_h', 'headers')]:
        built = work / label / 'db'
        built.parent.mkdir(parents=True, exist_ok=False)
        receipt = execute([executable, 'base:createsubdb', str(key_file), source_prefix + suffix, str(built), '--subdb-mode', '0', '--id-mode', '0'], log_dir / (label + '.log'))
        receipts.append(receipt)
        if receipt['exit_code'] != 0:
            raise ValueError('Official createsubdb failed for ' + label)
        if index_keys(str(built) + '.index') != keyset:
            raise ValueError('Exported component keys differ: ' + label)
        selected = [p for p in built.parent.iterdir() if p.name in ('db', 'db.index', 'db.dbtype') or re.fullmatch(r'db\.\d+', p.name)]
        if not any(p.name == 'db' or re.fullmatch(r'db\.\d+', p.name) for p in selected):
            raise ValueError('Materialized component payload missing: ' + label)
        for path in selected:
            if path.is_symlink() or not path.is_file():
                raise ValueError('Official mode 0 left a nonregular core component: ' + str(path))
            destination = Path(str(target_prefix) + suffix + path.name[2:])
            if destination.exists() or destination.is_symlink():
                raise ValueError('Refusing export overwrite: ' + str(destination))
            shutil.copyfile(path, destination)
    return receipts


def export_db(plan_path, output):
    plan_path, output = Path(plan_path), Path(output)
    plan = json.loads(plan_path.read_text())
    core.forbid_source_output_overlap(output, plan.get('archived_source'), plan.get('collection'))
    reasons = validate_plan(plan)
    report = {'created_utc': core.utc(), 'status': 'NOT_COMPLETE', 'blocking_reasons': reasons,
              'source_inputs_modified': False, 'official_materialization_only': True, 'whole_afsp_payload_copied': False}
    if reasons:
        core.write_json(output / 'PORTABLE_DB_REPORT.json', report)
        return report
    executable = plan['runtime']['path']
    help_proc = subprocess.run([executable, 'base:createsubdb', '-h'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=30)
    (output / 'CREATESUBDB_HELP.txt').write_text(help_proc.stdout)
    if not ('--subdb-mode' in help_proc.stdout and '0: copy data' in help_proc.stdout and '--id-mode' in help_proc.stdout):
        reasons.append('INSTALLED_OFFICIAL_SUBDB_COPY_MODE_SEMANTICS_NOT_CONFIRMED')
        core.write_json(output / 'PORTABLE_DB_REPORT.json', report)
        return report
    receipts, identities, snapshots, prefixes = [], [], {}, {}
    try:
        for mode, name, basename in [('POS', 'CLEAN_RT_REFERENCES', 'clean_rt_afsp'), ('NEG', 'MATCHED_NONRT_CONTROLS', 'matched_nonrt_afsp')]:
            roster_path = Path(plan['archived_source']) / 'frozen' / ('STAGE05D4R5C_' + name + '.tsv')
            roster = core.read_tsv(roster_path)
            key_file = output / (mode + '_FROZEN_KEYS.txt')
            key_file.write_text(''.join(k + '\n' for k in sorted({r['afsp_key'] for r in roster}, key=int)))
            original_prefix = plan['database_snapshots'][mode]['prefix']
            new_prefix = output / 'reference_portable' / basename
            receipts.extend(materialize_core_database(executable, original_prefix, key_file, new_prefix, output / 'work' / mode, output / 'logs' / mode))
            prefixes[mode] = str(new_prefix)
            for which, prefix in [('source', original_prefix), ('portable', str(new_prefix))]:
                receipt = execute([executable, 'convert2fasta', prefix, str(output / (mode + '_' + which + '.faa'))], output / 'logs' / (mode + '_' + which + '_convert2fasta.log'))
                receipts.append(receipt)
                if receipt['exit_code'] != 0:
                    raise ValueError('FASTA conversion failed: ' + mode + '/' + which)
            a, b = fasta(output / (mode + '_source.faa')), fasta(output / (mode + '_portable.faa'))
            expected = {r['afsp_target']: r for r in roster}
            if set(a) != set(expected) or set(b) != set(expected):
                raise ValueError('Source or portable FASTA names do not exactly match frozen AF names: ' + mode)
            for key, member in expected.items():
                match = a[key] == b[key] and len(a[key]) == int(member['aa_length'])
                identities.append({'cohort': mode, 'accession': member['accession'], 'afsp_target': key,
                                   'original_sequence_sha256': core.hashlib.sha256(a[key].encode()).hexdigest(),
                                   'portable_sequence_sha256': core.hashlib.sha256(b[key].encode()).hexdigest(),
                                   'length': len(a[key]), 'frozen_length': member['aa_length'], 'exact_match': match})
                if not match:
                    raise ValueError('Sequence or length mismatch: ' + member['accession'])
            snapshots[mode] = core.db_snapshot(new_prefix)
            if snapshots[mode]['current_prefix_snapshot_state'] != 'HASHED_WITH_CORE_PAYLOAD_PRESENT' or any(r['is_symlink'] for r in snapshots[mode]['components']):
                raise ValueError('Portable component inventory is incomplete or has external links')
        core.write_tsv(output / 'CONTROL_128_SEQUENCE_IDENTITY.tsv', identities, list(identities[0]))
        core.write_json(output / 'PORTABLE_DB_COMPONENTS.json', snapshots)
        # Round-trip search uses only the materialized four-component databases.
        replay = run_searches(plan, output / 'portable_replay', prefixes)
        core.write_json(output / 'portable_replay' / 'SEARCH_REPLAY_REPORT.json', replay)
        report['portable_search_replay'] = replay
        if replay['status'] != 'PASS_EXACT_SEMANTIC_SEARCH_REPLAY':
            reasons.append('PORTABLE_DATABASE_SEARCH_REPLAY_REQUIRES_REVIEW')
        reasons.extend(e for snapshot in plan['database_snapshots'].values() for e in unchanged(snapshot))
    except (OSError, ValueError, KeyError) as exc:
        reasons.append(str(exc))
    report.update({'status': 'PASS_PORTABLE_DB_AND_EXACT_REPLAY' if not reasons else 'NOT_COMPLETE',
                   'command_receipts': receipts, 'portable_prefixes': prefixes, 'sequence_identity_records': len(identities),
                   'database_component_snapshots': snapshots,
                   'plan_sha256': core.sha256(plan_path),
                   'source_database_snapshot_sha256': {mode: core.hashlib.sha256(json.dumps(s, sort_keys=True).encode()).hexdigest() for mode, s in plan['database_snapshots'].items()},
                   'distribution_boundary': 'The portable files contain the selected original AA, 3Di, C-alpha and header records generated by official mode-0 subsetting. Redistribution follows the original data licenses; this program assigns no new data license.'})
    core.write_json(output / 'COMMAND_RECEIPTS.json', receipts)
    core.write_json(output / 'PORTABLE_DB_REPORT.json', report)
    return report


def replay_portable(portable_root, source, output, threads=32):
    """Relocate by recorded component basenames and verify every component before searching."""
    portable_root, source, output = Path(portable_root).resolve(), Path(source).resolve(), Path(output)
    core.forbid_source_output_overlap(output, portable_root, source)
    errors = []
    report = {'status': 'NOT_COMPLETE', 'failures': errors, 'portable_root': str(portable_root), 'new_search_performed': False}
    metadata = portable_root / 'PORTABLE_DB_COMPONENTS.json'
    if not metadata.is_file():
        errors.append('PORTABLE_COMPONENT_HASH_MANIFEST_MISSING')
        core.write_json(output / 'SEARCH_REPLAY_REPORT.json', report)
        return report
    snapshots = json.loads(metadata.read_text())
    runtime = runtime_info()
    if runtime['state'] != 'ORIGINAL_VERSION':
        errors.append('ORIGINAL_FOLDSEEK_VERSION_REQUIRED')
    if threads < 1:
        errors.append('THREADS_MUST_BE_POSITIVE')
    prefixes = {}
    for mode, basename in [('POS', 'clean_rt_afsp'), ('NEG', 'matched_nonrt_afsp')]:
        if mode not in snapshots:
            errors.append(mode + '_SNAPSHOT_MISSING')
            continue
        prefix = portable_root / 'reference_portable' / basename
        prefixes[mode] = str(prefix)
        expected_names = set()
        for component in snapshots[mode]['components']:
            name = component['component']
            if Path(name).name != name or not name.startswith(basename):
                errors.append('INVALID_COMPONENT_NAME:' + name)
                continue
            expected_names.add(name)
            path = prefix.parent / name
            if path.is_symlink() or not path.is_file() or core.sha256(path) != component.get('sha256'):
                errors.append('PORTABLE_COMPONENT_MISSING_OR_HASH_MISMATCH:' + name)
        if {p.name for p in prefix.parent.glob(basename + '*')} != expected_names:
            errors.append(mode + '_PORTABLE_COMPONENT_SET_CHANGED')
    audit = core.audit_calibration(source, output / 'archived_source_audit')
    if audit['status'] != 'PASS_ARCHIVED_CONTROL_REAGGREGATION':
        errors.append('ORIGINAL_SOURCE_REAGGREGATION_FAILED')
    if not errors:
        plan = {'runtime': runtime, 'threads': threads, 'archived_source': str(source)}
        report = run_searches(plan, output, prefixes)
        report['new_search_performed'] = True
        report['portable_root'] = str(portable_root)
    core.write_json(output / 'SEARCH_REPLAY_REPORT.json', report)
    return report
