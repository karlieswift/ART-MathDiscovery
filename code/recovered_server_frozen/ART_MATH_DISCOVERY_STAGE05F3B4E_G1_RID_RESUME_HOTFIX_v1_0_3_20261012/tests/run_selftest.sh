#!/usr/bin/env bash
set -euo pipefail
S=$(cd "$(dirname "$0")/.." && pwd)
python -m py_compile "$S/scripts/ncbi_blast_urlapi_resume.py" "$S/scripts/summarize_global.py"
bash -n "$S/run_stage05f3b4e_global_urlapi_r3_resume.sh" "$S/check_rid_status.sh"
python - "$S/scripts/ncbi_blast_urlapi_resume.py" <<'PY'
import importlib.util,sys
p=sys.argv[1];sp=importlib.util.spec_from_file_location('m',p);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
xml='''<BlastXML2><BlastOutput2><report><Report><results><Results><search><Search><query-id>Q1</query-id><query-title>Q1</query-title><query-len>100</query-len><hits><Hit><len>120</len><description><HitDescr><id>ref|XP_1.1|</id><accession>XP_1</accession><title>x</title><sciname>y</sciname></HitDescr></description><hsps><Hsp><query-from>1</query-from><query-to>90</query-to><hit-from>2</hit-from><hit-to>91</hit-to><align-len>90</align-len><identity>81</identity><evalue>1e-20</evalue><bit-score>80</bit-score></Hsp></hsps></Hit></hits></Search></search></Results></results></Report></report></BlastOutput2></BlastXML2>'''
r=m.parse_xml2_to_rows(xml);assert len(r)==1 and abs(float(r[0][2])-90)<1e-9 and abs(float(r[0][12])-90)<1e-9
print('[PASS] RID-resume URLAPI parser/syntax selftest')
PY
