#!/usr/bin/env python3
import sys,csv,pathlib,hashlib,json
root=pathlib.Path(sys.argv[1]);result=pathlib.Path(sys.argv[2]);checks=[]
for r in csv.DictReader(open(root/'EXPECTED_OUTPUT_SHA256.tsv'),delimiter='\t'):
 name=r['filename']
 sub='N1' if name.startswith('N1') else 'N3' if name.startswith('N3') else 'N4'
 p=result/sub/name
 hash=hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else ''
 checks.append({'table':name,'match':hash==r['sha256']})
summary={'status':'PASS_SEVEN_REFERENCE_TABLES_BITWISE' if all(x['match'] for x in checks) else 'FAIL_REPRODUCIBILITY','n_tables':len(checks),'failures':[x['table'] for x in checks if not x['match']]}
print(json.dumps(summary,indent=2))
(result/'EXPECTED_REPLAY_REPORT.json').write_text(json.dumps(summary,indent=2)+'\n')
if summary['failures']:raise SystemExit(2)
