#!/usr/bin/env python3
import sys,json
n=0; bp=0; cur=False
for raw in sys.stdin:
    if raw.startswith('>'):
        n+=1; cur=True
    else:
        s=''.join(raw.split()).upper()
        if s: bp+=len(s)
print(json.dumps({'records':n,'bp':bp}))
