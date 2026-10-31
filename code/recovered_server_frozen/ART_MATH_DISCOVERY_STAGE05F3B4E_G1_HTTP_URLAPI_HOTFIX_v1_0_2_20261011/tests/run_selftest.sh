#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")/.." && pwd)
python - "$S/scripts/ncbi_blast_urlapi.py" <<'PY'
import importlib.util,sys
p=sys.argv[1]; spec=importlib.util.spec_from_file_location('u',p); u=importlib.util.module_from_spec(spec); spec.loader.exec_module(u)
xml='''<?xml version="1.0"?><BlastXML2 xmlns="http://www.ncbi.nlm.nih.gov"><BlastOutput2><report><Report><results><Results><search><Search><query-id>Query_1</query-id><query-title>Q1 family=x</query-title><query-len>100</query-len><hits><Hit><description><HitDescr><id>ref|XP_123.4|</id><accession>XP_123</accession><title>reverse transcriptase [Testus]</title><sciname>Testus</sciname></HitDescr></description><len>110</len><hsps><Hsp><bit-score>200</bit-score><evalue>1e-40</evalue><identity>72</identity><query-from>1</query-from><query-to>80</query-to><hit-from>3</hit-from><hit-to>82</hit-to><align-len>80</align-len></Hsp></hsps></Hit></hits></Search></search></Results></results></Report></report></BlastOutput2></BlastXML2>'''
r=u.parse_xml2_to_rows(xml)
assert len(r)==1,r
assert r[0][0]=='Q1',r[0]
assert r[0][1]=='XP_123.4',r[0]
assert abs(float(r[0][2])-90.0)<1e-6,r[0]
assert abs(float(r[0][12])-80.0)<1e-6,r[0]
print('[PASS] URL API XML2 -> frozen TSV schema parser selftest')
PY
python -m py_compile "$S/scripts/ncbi_blast_urlapi.py" "$S/scripts/summarize_global.py" "$S/scripts/common.py"
bash -n "$S/run_stage05f3b4e_global_urlapi_r2.sh"
echo '[PASS] F3B4E G1-R2 URLAPI hotfix syntax/policy selftest'
