from pathlib import Path
import sys,tempfile,csv
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/'scripts'))
from sixframe_orfs import sixframe
from inventory_and_freeze import canonicalize_proposals, canonicalize_union
seq='ATG'+'AAA'*50+'TAA'+'CCC'*10
assert list(sixframe(seq,45))
contract=(HERE/'config/BENCHMARK_CONTRACT.json').read_text(); assert 'truth_must_not_be_used_before_score_freeze' in contract
score=(HERE/'scripts/run_score_freeze.py').read_text(); assert "read_tsv(here/'frozen/BASELINE_TRUTH_SEAL.tsv')" not in score; assert 'PRODIGAL_RVT' in score and 'SIXFRAME_RVT' in score; assert "tool['tools']['prodigal']['path']" in score
assert 'C1_SELECTED_PROPOSALS.tsv' in (HERE/'config/SOURCE_POLICY.tsv').read_text()
assert 'CP001363.1' in (HERE/'config/CASES.tsv').read_text() and 'AP006725.1' in (HERE/'config/CASES.tsv').read_text()
with tempfile.TemporaryDirectory() as td:
    td=Path(td); p=td/'p.tsv'; u=td/'u.tsv'
    p.write_text('family\tstart0\tend0\tqvalue\nACNV6MER\t10\t20\t0.001\nRECURRENCE\t15\t25\t0.002\n')
    u.write_text('rank\tstart0\tend0\n1\t9\t26\n')
    props,ps=canonicalize_proposals(p,td/'cp.tsv','X'); unions,us=canonicalize_union(u,td/'cu.tsv','X',props)
    assert ps['acnv_rows']==1 and ps['recurrence_rows']==1 and unions[0]['dual_arm']==1
print('[PASS] Stage05F4A R1 syntax/policy/canonical-board selftest')
print('[PASS] score-freeze runner cannot read truth seal by construction')
print('[PASS] corrected control-source accessions + C1_SELECTED_PROPOSALS support frozen')
