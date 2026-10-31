#!/usr/bin/env python3
"""Numerical, missingness, corruption and real-tool integration checks. Fixture data are not scientific results."""
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import j1_reproducibility as j
import j1_execution as execution

RUNTIME = HERE / 'runtime'
RUNTIME.mkdir(exist_ok=True)
ARCHIVE = HERE.parent / 'reference' / 'r5c_archived_results'


class J1Tests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=RUNTIME)
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_all_hit_counterexample_and_boundaries(self):
        rows = []
        for tm, q, t, name in [(0.90, 0.10, 0.80, 'highest_tm'), (0.50, 0.50, 0.50, 'qualifying_strong')]:
            rows.append({'alntmscore': tm, 'qcov': q, 'tcov': t, 'tier': j.hit_tier(tm, q, t), 'evalue': 0.01, 'target': name})
        agg = j.aggregate_hits(rows)
        self.assertEqual(agg['highest_tm_only_tier'], 'WEAK')
        self.assertEqual(agg['max_qualifying_tier'], 'STRONG')
        self.assertEqual(agg['representative']['target'], 'qualifying_strong')
        self.assertEqual(j.hit_tier(.49999, .5, .5), 'SUGGESTIVE')
        self.assertEqual(j.hit_tier(.4, .3, .3), 'SUGGESTIVE')
        self.assertEqual(j.hit_tier(.39999, .9, .9), 'WEAK')
        self.assertEqual(j.hit_tier(.8, .29999, .8), 'WEAK')

    def test_missing_empty_and_invalid_are_distinct(self):
        missing, empty, bad = [self.root / name for name in ('missing.m8', 'empty.m8', 'bad.m8')]
        empty.write_bytes(b'')
        bad.write_text('q\tt\t1e-9\t100\t0.8\t0.8\t100\t120\t120\t0.5\tNaN\n')
        self.assertEqual(j.parse_hits(missing)['state'], 'MISSING')
        self.assertEqual(j.parse_hits(empty)['state'], 'EMPTY')
        self.assertEqual(j.parse_hits(bad)['state'], 'INVALID')
        bad.write_text('q\tt\t1e-9\t100\t80\t0.8\t100\t120\t120\t0.5\t0.9\n')
        self.assertEqual(j.parse_hits(bad)['state'], 'INVALID')
        bad.write_text('q\tt\t1e-9\t100\t0.8\t0.8\t100\t120\t120\t0.5\t0.9\n')
        self.assertEqual(j.parse_hits(bad, expected_query='other')['state'], 'INVALID')

    def test_result_index_proves_completed_empty_and_rejects_absence(self):
        prefix = self.root / 'result'
        Path(str(prefix) + '.0').write_bytes(b'\0' + b'7\talignment\n\0')
        Path(str(prefix) + '.index').write_text('1\t0\t1\n2\t1\t13\n')
        roster = [{'afsp_key': '1'}, {'afsp_key': '2'}]
        evidence = j.result_db_evidence(prefix, roster)
        self.assertEqual(evidence['state'], 'COMPLETE_RESULT_RECORD_COVERAGE')
        self.assertEqual(evidence['records']['1']['result_state'], 'COMPLETED_EMPTY')
        self.assertEqual(evidence['records']['2']['raw_alignment_rows'], 1)
        self.assertEqual(j.result_db_evidence(prefix, roster + [{'afsp_key': '3'}])['state'], 'NOT_COMPLETE')
        Path(str(prefix) + '.0').unlink()
        self.assertEqual(j.result_db_evidence(prefix, roster)['state'], 'NOT_COMPLETE')

    def test_exact_stage_ambiguity_and_empty_output_preservation(self):
        stage = self.root / j.R5C_STAGE
        (stage / 'raw').mkdir(parents=True)
        (stage / 'raw' / 'empty.m8').write_bytes(b'')
        output = self.root / 'exports' / 'collection'
        report = j.collect(self.root, output)
        self.assertEqual(report['copied_zero_byte_files'], 1)
        self.assertTrue((output / 'source' / j.R5C_STAGE / 'raw' / 'empty.m8').is_file())
        self.assertEqual(report['database_snapshot_state'], 'NOT_COMPLETE')
        (self.root / 'another' / j.R5C_STAGE).mkdir(parents=True)
        ambiguous = j.collect(self.root, self.root / 'exports' / 'ambiguous')
        self.assertEqual(ambiguous['reason'], 'AMBIGUOUS_EXACT_STAGE')
        self.assertEqual(ambiguous['status'], 'NOT_COMPLETE')

    def test_shared_database_target_is_hashed_once_and_changes_detected(self):
        payload = self.root / 'shared_payload'
        payload.write_bytes(b'fixture payload')
        prefix = self.root / 'db'
        for suffix in ('', '_ss', '_ca', '_h'):
            Path(str(prefix) + suffix).symlink_to(payload)
            Path(str(prefix) + suffix + '.index').write_text('1\t0\t16\n')
            Path(str(prefix) + suffix + '.dbtype').write_bytes(b'\0\0\0\0')
        cache = {}
        snapshot = j.db_snapshot(prefix, cache)
        self.assertEqual(snapshot['current_prefix_snapshot_state'], 'HASHED_WITH_CORE_PAYLOAD_PRESENT')
        self.assertEqual(sum(r['hash_reused_for_identical_target'] for r in snapshot['components']), 3)
        self.assertFalse(execution.unchanged(snapshot))
        payload.write_bytes(b'altered payload')
        self.assertTrue(execution.unchanged(snapshot))

    def test_original_archive_numerical_result_and_corrupt_inputs(self):
        output = self.root / 'audit'
        report = j.audit_calibration(ARCHIVE, output)
        self.assertEqual(report['status'], 'PASS_ARCHIVED_CONTROL_REAGGREGATION')
        self.assertEqual(report['cohorts']['POS']['original_highest_tm_tiers'], {'STRONG': 29, 'WEAK': 2, 'SUGGESTIVE': 1})
        self.assertEqual(report['cohorts']['POS']['all_hit_tiers'], {'STRONG': 31, 'SUGGESTIVE': 1})
        self.assertEqual(report['cohorts']['NEG']['all_hit_tiers'], {'NO_HIT': 53, 'WEAK': 43})
        altered = self.root / 'altered'
        shutil.copytree(ARCHIVE, altered)
        (altered / 'raw' / 'NEG_VS_CLEAN.m8').write_bytes(b'')
        self.assertEqual(j.audit_calibration(altered, self.root / 'bad_empty_audit')['status'], 'NOT_COMPLETE')
        (altered / 'raw' / 'NEG_VS_CLEAN.m8').unlink()
        self.assertEqual(j.audit_calibration(altered, self.root / 'missing_audit')['status'], 'NOT_COMPLETE')

    def test_output_location_and_no_overwrite(self):
        with self.assertRaises(ValueError):
            j.make_output(self.root, self.root / 'unapproved-location', 'test')
        out = j.make_output(self.root, None, 'test')
        self.assertEqual(out.parent, self.root / 'exports')
        (out / 'preserve').write_text('preserve original')
        with self.assertRaises(ValueError):
            j.make_output(self.root, out, 'test')
        other_project = self.root / 'other_project'
        other_project.mkdir()
        outside = self.root / 'outside'
        outside.mkdir()
        (other_project / 'exports').symlink_to(outside, target_is_directory=True)
        with self.assertRaises(ValueError):
            j.make_output(other_project, other_project / 'exports' / 'forbidden', 'test')
        with self.assertRaises(ValueError):
            j.make_output(self.root, self.root / 'exports' / 'source' / 'nested', 'test', [self.root / 'exports' / 'source'])

    def test_component_addition_and_preserved_mtime_rewrite_blocked(self):
        prefix = self.root / 'db'
        prefix.write_bytes(b'abcdef')
        snapshot = j.db_snapshot(prefix)
        before = prefix.stat()
        prefix.write_bytes(b'ABCDEF')
        os.utime(prefix, ns=(before.st_atime_ns, before.st_mtime_ns))
        self.assertTrue(execution.unchanged(snapshot))
        fresh = j.db_snapshot(prefix)
        Path(str(prefix) + '.idx').write_bytes(b'new search index')
        self.assertTrue(execution.unchanged(fresh))

    def test_changed_roster_and_rewritten_local_checksums_cannot_pass(self):
        altered = self.root / 'altered'
        shutil.copytree(ARCHIVE, altered)
        path = altered / 'frozen' / 'STAGE05D4R5C_MATCHED_NONRT_CONTROLS.tsv'
        path.write_text(path.read_text().replace('Q9TY96', 'Q9TY96_RELABELED'))
        for filename in ('frozen/STAGE05D4R5C_FREEZE_SHA256.tsv', 'outputs/STAGE05D4R5C_SHA256.tsv'):
            manifest = altered / filename
            changed = []
            for line in manifest.read_text().splitlines():
                _, rel = line.split('  ', 1)
                changed.append(j.sha256(altered / rel) + '  ' + rel)
            manifest.write_text('\n'.join(changed) + '\n')
        result = j.audit_calibration(altered, self.root / 'audit')
        self.assertEqual(result['status'], 'NOT_COMPLETE')
        self.assertTrue(any('PINNED_ORIGINAL_ARCHIVE_MISMATCH' in s for s in result['failures']))

    def test_incomplete_plan_never_runs_searches(self):
        report = execution.plan_replay(self.root / 'missing_collection', ARCHIVE, self.root / 'plan')
        self.assertEqual(report['status'], 'NOT_COMPLETE')
        replay = execution.run_replay(self.root / 'plan' / 'REPLAY_PLAN.json', self.root / 'run')
        self.assertFalse(replay['search_performed'])
        self.assertEqual(replay['status'], 'NOT_COMPLETE')


class OfficialFoldseekIntegration(unittest.TestCase):
    def test_actual_copy_mode_portability_and_search_equivalence(self):
        executable = os.environ.get('J1_TEST_FOLDSEEK')
        cif_root = os.environ.get('J1_TEST_CIF_ROOT')
        if not executable or not cif_root:
            self.skipTest('Set J1_TEST_FOLDSEEK and J1_TEST_CIF_ROOT for the real official-tool integration check.')
        with tempfile.TemporaryDirectory(dir=RUNTIME) as tmp:
            root = Path(tmp)
            inputs = root / 'input_structures'
            inputs.mkdir()
            for name in ('SMQ0001.cif', 'SMQ0002.cif'):
                shutil.copyfile(Path(cif_root) / name, inputs / name)
            source = root / 'original_db' / 'db'
            source.parent.mkdir()
            self.assertEqual(execution.execute([executable, 'createdb', str(inputs), str(source)], root / 'logs' / 'createdb.log')['exit_code'], 0)
            chosen = sorted(execution.index_keys(str(source) + '.index'), key=int)[0]
            keys = root / 'one_key.txt'
            keys.write_text(chosen + '\n')
            soft = root / 'soft_db' / 'db'
            soft.parent.mkdir()
            self.assertEqual(execution.execute([executable, 'createsubdb', str(keys), str(source), str(soft), '--subdb-mode', '1'], root / 'logs' / 'soft_subset.log')['exit_code'], 0)
            before = j.db_snapshot(soft)
            self.assertTrue(any(r['is_symlink'] for r in before['components']))
            portable = root / 'portable_db' / 'db'
            execution.materialize_core_database(executable, soft, keys, portable, root / 'materialize', root / 'logs' / 'materialize')
            snapshot = j.db_snapshot(portable)
            self.assertEqual(snapshot['current_prefix_snapshot_state'], 'HASHED_WITH_CORE_PAYLOAD_PRESENT')
            self.assertFalse(any(r['is_symlink'] for r in snapshot['components']))
            outputs = []
            for label, prefix in [('original', soft), ('portable', portable)]:
                run = root / label
                (run / 'raw').mkdir(parents=True)
                (run / 'work').mkdir()
                for command in execution.commands(executable, str(prefix), str(prefix), 1, run):
                    self.assertEqual(execution.execute(command['argv'], root / 'logs' / (label + '_' + command['name'] + '.log'))['exit_code'], 0)
                outputs.append((run / 'raw' / 'POS_VS_CLEAN.m8').read_bytes())
                self.assertEqual(execution.execute([executable, 'convert2fasta', str(prefix), str(run / 'sequence.faa')], root / 'logs' / (label + '_fasta.log'))['exit_code'], 0)
            self.assertEqual(outputs[0], outputs[1])
            self.assertTrue(outputs[0])
            self.assertEqual(execution.fasta(root / 'original' / 'sequence.faa'), execution.fasta(root / 'portable' / 'sequence.faa'))
            self.assertFalse(execution.unchanged(before))
            core_paths = [Path(r['path']) for r in snapshot['components']]
            for path in core_paths:
                self.assertTrue(j.within(path.resolve(), portable.parent))


if __name__ == '__main__':
    suite = unittest.defaultTestLoader.loadTestsFromModule(sys.modules[__name__])
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    summary = {'tests_run': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors),
               'skipped': len(result.skipped), 'successful': result.wasSuccessful(),
               'scientific_boundary': 'Synthetic and real-tool interface tests only; archived-control audit tests replay existing data, not a new experiment.',
               'official_foldseek_executable': os.environ.get('J1_TEST_FOLDSEEK', ''),
               'official_integration_executed': bool(os.environ.get('J1_TEST_FOLDSEEK') and os.environ.get('J1_TEST_CIF_ROOT'))}
    (HERE / 'TEST_REPORT.json').write_text(json.dumps(summary, indent=2) + '\n')
    raise SystemExit(0 if result.wasSuccessful() else 1)
