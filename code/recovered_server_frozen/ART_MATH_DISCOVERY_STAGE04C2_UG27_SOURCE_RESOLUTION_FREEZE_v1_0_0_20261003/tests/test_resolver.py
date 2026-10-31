#!/usr/bin/env python3
from pathlib import Path
import tempfile, subprocess, json
import openpyxl
ROOT=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory() as td:
    td=Path(td); (td/'inputs').mkdir(); (td/'outputs').mkdir()
    wb=openpyxl.Workbook(); ws=wb.active; ws.title='UG27'
    for i,acc in enumerate(['CP100001.1','CP100002.1','CP100003.1','CP100004.1','CP100005.1'],1):
        ws.append([f'UG27 experimental system {i}', acc, 'cDNA hairpin'])
    x=td/'inputs'/'Minerva_SupplementaryTables_preprint.xlsx'; wb.save(x)
    p=subprocess.run(['python',str(ROOT/'scripts'/'resolve_ug27_sources.py'),'--root',str(td)],capture_output=True,text=True)
    assert p.returncode==0, (p.returncode,p.stdout,p.stderr)
    g=json.loads((td/'outputs'/'GATE.json').read_text())
    assert g['gate_pass'] is True and g['explicit_ug27_experimental_identifier_n']==5
print('PASS')
