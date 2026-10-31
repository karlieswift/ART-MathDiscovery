from pathlib import Path
import sys,tarfile
from io_utils import manifest_tree
root=Path(sys.argv[1]);out=root/sys.argv[2]
if (root/'.run_lock').exists():raise SystemExit('Refusing to package an active run')
paths=['outputs','logs','config','docs','evidence','tests','scripts','README.md','PACKAGE_MANIFEST.tsv','PACKAGE_MANIFEST.sha256']
paths += [p.name for p in sorted(root.glob('*.sh'))]
if (root/'inputs').exists():paths.append('inputs')
manifest_tree(root/'outputs',root/'outputs'/'RETURNED_OUTPUTS_MANIFEST.tsv',exclude=['RETURNED_OUTPUTS_MANIFEST.tsv'])
tmp=out.with_name(out.name+'.partial')
try:
    with tarfile.open(tmp,'w:gz') as t:
        for rel in paths:
            p=root/rel
            if not p.exists():raise FileNotFoundError(p)
            t.add(p,arcname=rel,filter=lambda info: None if '__pycache__' in info.name or info.name.endswith('.pyc') else info)
    with tarfile.open(tmp,'r:gz') as t:
        for m in t:
            if m.isfile():
                f=t.extractfile(m)
                while f.read(1<<20):pass
    tmp.rename(out)
except Exception:
    if tmp.exists():tmp.unlink()
    raise
print(out)
