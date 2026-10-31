#!/usr/bin/env python3
import importlib.util, pathlib, tempfile
p=pathlib.Path(__file__).resolve().parents[1]/'scripts/stage05e3h.py'
s=importlib.util.spec_from_file_location('m',p); m=importlib.util.module_from_spec(s); s.loader.exec_module(m)
# Exact SHA and feature primitives.
seq='MAHAHQRPPEAASEETEQTPQRTNAMKTLIERLRAAHLEAEAERVAQYLSAQALEA'
r=m.seq_features('ARCH250_0010','x witness=True',seq,True)
assert r['aa_length']==56 and r['aa_sha256']=='196925f6919925ec4b8f70c4c54ac4956a656a9a61543838655ee1427cd76268'
assert r['strong_hydrophobic_segment_heuristic'] in ('True','False')
assert m.classify_existing(['x NO_SWISSPROT_HIT_NOT_NEGATIVE y'],'D4_R2_SWISSPROT')=='NO_SWISSPROT_HIT_NOT_NEGATIVE'
assert m.classify_existing(['x TIER_S_SHORT_LT180 y'],'D4_R4_PROSTT5_FOLDSEEK')=='NOT_SEARCHED_OLD_D4_SHORT_LT180'
assert m.classify_existing(['x NO_FULL3D_AFSP_HIT_NOT_NEGATIVE y'],'D4_R5_FULL3D')=='NO_FULL3D_CURATED_OR_RT_REFERENCE_HIT_NOT_NEGATIVE'
print('[PASS] Stage05E3H sequence-feature/evidence-semantics selftest')
