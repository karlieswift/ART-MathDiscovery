#!/usr/bin/env python3
from pathlib import Path
import tempfile, subprocess, json, zipfile, re
import openpyxl
ROOT=Path(__file__).resolve().parents[1]

def run_case(xlsx: Path, td: Path):
    p=subprocess.run(['python',str(ROOT/'scripts'/'resolve_ug27_sources.py'),'--root',str(td)],capture_output=True,text=True)
    assert p.returncode==0, (p.returncode,p.stdout,p.stderr)
    g=json.loads((td/'outputs'/'GATE.json').read_text())
    assert g['gate_pass'] is True and g['explicit_ug27_experimental_identifier_n']==5, g

def make_workbook(path: Path):
    wb=openpyxl.Workbook(); ws=wb.active; ws.title='UG27'
    for i,acc in enumerate(['CP100001.1','CP100002.1','CP100003.1','CP100004.1','CP100005.1'],1):
        ws.append([f'UG27 experimental system {i}', acc, 'cDNA hairpin'])
    wb.save(path)

def strip_dimension(src: Path, dst: Path):
    # Reproduce the real Minerva workbook condition that can make max_row=None
    # in openpyxl read_only mode by removing worksheet <dimension .../> tags.
    with zipfile.ZipFile(src,'r') as zin, zipfile.ZipFile(dst,'w',zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data=zin.read(item.filename)
            if item.filename.startswith('xl/worksheets/sheet') and item.filename.endswith('.xml'):
                txt=data.decode('utf-8')
                txt=re.sub(r'<dimension\s+ref="[^"]+"\s*/>', '', txt, count=1)
                data=txt.encode('utf-8')
            zout.writestr(item, data)

with tempfile.TemporaryDirectory() as raw:
    td=Path(raw); (td/'inputs').mkdir(); (td/'outputs').mkdir()
    normal=td/'normal.xlsx'; make_workbook(normal)
    x=td/'inputs'/'Minerva_SupplementaryTables_preprint.xlsx'; normal.replace(x)
    run_case(x, td)

with tempfile.TemporaryDirectory() as raw:
    td=Path(raw); (td/'inputs').mkdir(); (td/'outputs').mkdir()
    normal=td/'normal.xlsx'; make_workbook(normal)
    x=td/'inputs'/'Minerva_SupplementaryTables_preprint.xlsx'; strip_dimension(normal, x)
    wb=openpyxl.load_workbook(x, read_only=True, data_only=True)
    assert wb['UG27'].max_row is None, wb['UG27'].max_row
    run_case(x, td)

print('PASS')
