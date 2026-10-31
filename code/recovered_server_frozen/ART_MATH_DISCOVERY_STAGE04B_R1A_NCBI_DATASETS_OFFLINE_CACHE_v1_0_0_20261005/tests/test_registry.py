from pathlib import Path
import hashlib
root=Path(__file__).resolve().parents[1]
p=root/'inputs/STAGE04A_R1_DUAL_ARM_REGIONS.tsv'
h=hashlib.sha256(p.read_bytes()).hexdigest()
assert h=='73007864fdef7dba662c9b90a02e2b915359c1fbd93a6061edeea35afde07349', h
acc=[x.strip() for x in (root/'inputs/STAGE04B_R1A_ACCESSIONS_174.txt').read_text().splitlines() if x.strip()]
assert len(acc)==174 and len(set(acc))==174
assert 'NC_002666.1' in acc
print('[PASS] Stage04B-R1A registry selftest: 383 frozen regions -> 174 exact versioned accessions')
