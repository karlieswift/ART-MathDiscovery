#!/usr/bin/env python3
from pathlib import Path
import hashlib, tarfile
from datetime import datetime
R=Path(__file__).resolve().parents[1]
name=R.name+'_RESULTS_'+datetime.now().strftime('%Y%m%d_%H%M%S')+'.tar.gz'
p=R/name
items=['outputs','logs','GATE_STAGE04B_R1E.json','OUTPUT_MANIFEST.tsv','config','docs','scripts','tests','refs','inputs/controls','inputs/candidates.tsv','inputs/STAGE04A_UNIVERSE_MANIFEST.tsv','inputs/PRIOR_R1D2_GATE.json','PACKAGE_MANIFEST.tsv','PACKAGE_MANIFEST.sha256','README.md']
# Include new 66 predictions only. Earlier 174 are already frozen in the prior returned package.
for d in sorted((R/'work/prodigal').glob('*')):
 if d.is_dir() and (d/'R1E_COMPLETE.json').exists():
  import json
  if not json.loads((d/'R1E_COMPLETE.json').read_text())['fingerprint']['reused_R1D2']:
   items.append(str(d.relative_to(R)))
   items.append('inputs/genomes/'+d.name+'.fna')
if not (R/'GATE_STAGE04B_R1E.json').exists() and not (R/'outputs').exists():raise SystemExit('No execution results to collect')
with tarfile.open(p,'w:gz') as tf:
 for item in items:
  z=R/item
  if z.exists():tf.add(z,arcname=item)
h=hashlib.sha256(p.read_bytes()).hexdigest();p.with_name(p.name+'.sha256').write_text(h+'  '+p.name+'\n')
print(p);print(p.stat().st_size);print(h)
