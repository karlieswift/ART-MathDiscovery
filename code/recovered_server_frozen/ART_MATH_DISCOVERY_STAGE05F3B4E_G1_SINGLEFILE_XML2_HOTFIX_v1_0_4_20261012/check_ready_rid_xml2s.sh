#!/usr/bin/env bash
set -euo pipefail
RID=${1:?usage: $0 RID}
python - "$RID" <<'PY'
import sys,urllib.parse,urllib.request
rid=sys.argv[1]
url='https://blast.ncbi.nlm.nih.gov/Blast.cgi'
p={'CMD':'Get','RID':rid,'FORMAT_TYPE':'XML2_S'}
req=urllib.request.Request(url,data=urllib.parse.urlencode(p).encode(),headers={'User-Agent':'ART-MathDiscovery-F3B4E/1.0.4','Content-Type':'application/x-www-form-urlencoded'})
with urllib.request.urlopen(req,timeout=180) as r:b=r.read()
print('bytes=',len(b))
s=b.decode('utf-8','replace')
print(s[:500].replace('\n',' '))
print('single_file_xml2=',('<BlastXML2' in s or '<BlastOutput2' in s))
PY
