#!/usr/bin/env python3
import subprocess,tempfile,json,sys
from pathlib import Path
R=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory() as td:
    td=Path(td); fa=td/'x.fa'; gb=td/'x.gb'; out=td/'o.json'
    seq='A'*1000001
    fa.write_text('>CPTEST.1 synthetic\n'+seq+'\n')
    header='LOCUS       CPTEST 1000001 bp DNA\nACCESSION   CPTEST\nVERSION     CPTEST.1\n'
    gb.write_text(header + ('COMMENT     synthetic validator padding\n'*40))
    subprocess.run([sys.executable,str(R/'scripts/validate_record.py'),'--control-id','X','--registry-accession','NC_X.1','--fetch-accession','CPTEST.1','--relation','TEST','--expected-length-bp','1000001','--fasta',str(fa),'--flatfile',str(gb),'--out',str(out)],check=True,stdout=subprocess.DEVNULL)
    assert json.loads(out.read_text())['gate_pass']
print('[PASS] Stage04C0-R2 synthetic record validator selftest')
