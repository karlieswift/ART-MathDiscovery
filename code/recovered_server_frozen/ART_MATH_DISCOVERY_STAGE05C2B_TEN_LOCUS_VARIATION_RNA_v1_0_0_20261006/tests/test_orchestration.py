"""I/O orchestration only: a TEST DOUBLE is never scientific RNA evidence."""
from pathlib import Path
import unittest,tempfile,os,sys,shutil,contextlib,io
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from common import *
from fold import fold_all
class ExecutionTests(unittest.TestCase):
 def test_mock_execution_resume_and_damage_detection(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)
   for d in ['config','precomputed','bin']:(p/d).mkdir()
   pol=read_json(ROOT/'config/POLICY.json');pol['RNAfold']['null_replicates']=3
   write_json(p/'config/POLICY.json',pol)
   seqs={'RNA_'+sha_text('ACGU'*30):'ACGU'*30,'RNA_'+sha_text('GCAU'*15):'GCAU'*15}
   write_fasta(p/'precomputed/RNA_TARGETS.fasta',seqs)
   write_tsv(p/'precomputed/RNA_OCCURRENCE_REGISTRY.tsv',[{'occurrence_id':'O1','target_id':next(iter(seqs))}])
   write_json(p/'precomputed/INPUT_FREEZE.json',{'RNA_targets_sha256':sha_file(p/'precomputed/RNA_TARGETS.fasta'),'RNA_occurrences_sha256':sha_file(p/'precomputed/RNA_OCCURRENCE_REGISTRY.tsv'),'policy_sha256':sha_file(p/'config/POLICY.json')})
   fake=p/'bin/RNAfold';fake.write_text('#!'+sys.executable+'\nimport sys\nif "--version" in sys.argv:\n print("RNAfold 2.7.2 TEST_DOUBLE_NOT_SCIENTIFIC")\nelse:\n k=None;b=[]\n def emit():\n  if k:\n   s="".join(b);print(k);print(s);print("."*len(s)+" (0.00)")\n for l in sys.stdin:\n  l=l.strip()\n  if l.startswith(">"):\n   emit();k=l;b=[]\n  elif l:b.append(l)\n emit()\n');fake.chmod(0o755)
   with patch.dict(os.environ,{'PATH':str(p/'bin')+os.pathsep+os.environ['PATH']}),contextlib.redirect_stdout(io.StringIO()):
    fold_all(p,2);before=(p/'outputs/RNA_TARGET_FOLD_RESULTS.tsv').read_bytes()
    gate=read_json(p/'outputs/RNA_FOLD_GATE.json');self.assertTrue(gate['gate_pass']);self.assertEqual(gate['n_completed_targets'],2)
    fold_all(p,2);self.assertEqual(before,(p/'outputs/RNA_TARGET_FOLD_RESULTS.tsv').read_bytes())
    first=next((p/'work/fold').glob('*/stdout.txt'));first.write_text('corrupted TEST DOUBLE output')
    with self.assertRaises(RuntimeError):fold_all(p,2)
 def test_missing_binary_does_not_create_success(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td);(p/'config').mkdir();(p/'precomputed').mkdir()
   shutil.copyfile(ROOT/'config/POLICY.json',p/'config/POLICY.json')
   for n in ['INPUT_FREEZE.json','RNA_TARGETS.fasta','RNA_OCCURRENCE_REGISTRY.tsv']:shutil.copyfile(ROOT/'precomputed'/n,p/'precomputed'/n)
   with patch('shutil.which',return_value=None):
    with self.assertRaises(RuntimeError):fold_all(p,1)
   self.assertFalse((p/'outputs/RNA_FOLD_GATE.json').exists())
