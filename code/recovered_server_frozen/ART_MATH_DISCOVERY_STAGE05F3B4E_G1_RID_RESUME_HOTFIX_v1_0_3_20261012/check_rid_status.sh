#!/usr/bin/env bash
set -euo pipefail
RID=${1:?usage: check_rid_status.sh RID}
python - "$RID" <<'PY'
import sys,urllib.parse,urllib.request,re
rid=sys.argv[1]
data=urllib.parse.urlencode({'CMD':'Get','RID':rid,'FORMAT_OBJECT':'SearchInfo','TOOL':'ART_MathDiscovery_F3B4E'}).encode()
req=urllib.request.Request('https://blast.ncbi.nlm.nih.gov/Blast.cgi',data=data,headers={'User-Agent':'ART-MathDiscovery-F3B4E/1.0.3','Content-Type':'application/x-www-form-urlencoded'})
with urllib.request.urlopen(req,timeout=180) as r: body=r.read().decode('utf-8','replace')
m=re.search(r'Status\s*=\s*([A-Za-z]+)',body,re.I)
print('RID =',rid)
print('STATUS =',m.group(1).upper() if m else 'UNPARSED')
for key in ['ThereAreHits','QBlastInfoBegin','QBlastInfoEnd']:
    if key in body: print(body[body.find(key):body.find(key)+300].replace('\n',' '))
PY
