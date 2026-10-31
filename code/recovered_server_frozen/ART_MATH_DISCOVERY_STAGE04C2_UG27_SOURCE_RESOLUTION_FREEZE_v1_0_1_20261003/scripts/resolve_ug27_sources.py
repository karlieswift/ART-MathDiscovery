#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json, re, sys
from pathlib import Path
from typing import Iterable

try:
    import openpyxl
except Exception as e:
    print(f"[ERROR] openpyxl unavailable: {e}", file=sys.stderr)
    sys.exit(2)

ACC_PATTERNS = [
    re.compile(r"\b(?:NC|NZ|CP|AP|AE|AB|AJ|AM|FM|FN|FO|FP|HE|HG|JQ|JR|KC|KF|KJ|KM|KP|KT|KU|KX|KY|LC|LR|LT|MF|MG|MH|MK|MN|MT|MW|MZ|OK|OL|OM|ON|OP|OR|OQ|PQ)_[0-9]+(?:\.[0-9]+)?\b", re.I),
    re.compile(r"\b(?:CP|AP|AE|AB|AJ|AM|FM|FN|FO|FP|HE|HG|JQ|JR|KC|KF|KJ|KM|KP|KT|KU|KX|KY|LC|LR|LT|MF|MG|MH|MK|MN|MT|MW|MZ|OK|OL|OM|ON|OP|OR|OQ|PQ)[0-9]{5,}(?:\.[0-9]+)?\b", re.I),
    re.compile(r"\bGCA_[0-9]+(?:\.[0-9]+)?\b", re.I),
    re.compile(r"\bGCF_[0-9]+(?:\.[0-9]+)?\b", re.I),
]
KEYWORDS = ("ug27", "unknown group 27", "experimental", "experimentally", "cdna", "hairpin", "prophage")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def clean(v) -> str:
    if v is None:
        return ""
    return str(v).strip().replace("\n", " ")


def extract_ids(text: str) -> list[str]:
    out=[]
    for p in ACC_PATTERNS:
        out.extend(m.group(0) for m in p.finditer(text))
    return sorted(set(out), key=str.upper)


def tsv_write(path: Path, rows: Iterable[dict], fields: list[str]) -> None:
    import csv
    with path.open('w', newline='', encoding='utf-8') as f:
        w=csv.DictWriter(f, delimiter='\t', fieldnames=fields, extrasaction='ignore')
        w.writeheader(); w.writerows(rows)


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--root', default='.')
    ap.add_argument('--xlsx', default='inputs/Minerva_SupplementaryTables_preprint.xlsx')
    args=ap.parse_args()
    root=Path(args.root).resolve(); xlsx=(root/args.xlsx).resolve() if not Path(args.xlsx).is_absolute() else Path(args.xlsx)
    out=root/'outputs'; out.mkdir(parents=True, exist_ok=True)
    if not xlsx.exists():
        gate={"gate":"HOLD_STAGE04C2_UG27_EXACT_SOURCE_IDS_UNRESOLVED","gate_pass":False,"reason":"supplementary XLSX missing","xlsx":str(xlsx),"scoring_authorized":False}
        (out/'GATE.json').write_text(json.dumps(gate, indent=2), encoding='utf-8')
        print(json.dumps(gate, indent=2)); return 3
    wb=openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    hits=[]; id_rows=[]
    for ws in wb.worksheets:
        for r_idx,row in enumerate(ws.iter_rows(values_only=True), start=1):
            vals=[clean(v) for v in row]
            rowtext=" | ".join(v for v in vals if v)
            low=rowtext.lower()
            if any(k in low for k in KEYWORDS):
                ids=extract_ids(rowtext)
                hits.append({"sheet":ws.title,"row":r_idx,"text":rowtext,"accession_like_ids":",".join(ids)})
                for ident in ids:
                    id_rows.append({"identifier":ident,"sheet":ws.title,"row":r_idx,"context":rowtext,"ug27_context":("ug27" in low or "unknown group 27" in low),"experimental_context":("experiment" in low)})
    # Include neighboring rows around direct UG27 hits, because identifiers may be in adjacent cells/rows.
    direct=[h for h in hits if ('ug27' in h['text'].lower() or 'unknown group 27' in h['text'].lower())]
    for h in direct:
        ws=wb[h['sheet']]
        # In openpyxl read_only mode, worksheets whose XML omits a <dimension>
        # element can legitimately report ws.max_row=None. Do not use min(None, ...).
        # Cap by max_row only when it is a valid integer; iter_rows safely stops
        # at the physical end of the sheet otherwise.
        ws_max_row = ws.max_row if isinstance(ws.max_row, int) and ws.max_row >= 1 else None
        lo = max(1, h['row'] - 3)
        hi = h['row'] + 3 if ws_max_row is None else min(ws_max_row, h['row'] + 3)
        for rr in range(lo, hi + 1):
            row_iter = ws.iter_rows(min_row=rr, max_row=rr, values_only=True)
            row_vals = next(row_iter, None)
            if row_vals is None:
                continue
            vals=[clean(v) for v in row_vals]
            text=" | ".join(v for v in vals if v)
            for ident in extract_ids(text):
                id_rows.append({"identifier":ident,"sheet":ws.title,"row":rr,"context":text,"ug27_context":True,"experimental_context":("experiment" in text.lower())})
    # Deduplicate identifier+context location.
    dedup={}
    for r in id_rows:
        dedup[(r['identifier'].upper(),r['sheet'],r['row'])]=r
    id_rows=list(dedup.values())
    tsv_write(out/'ug27_keyword_hits.tsv', hits, ['sheet','row','text','accession_like_ids'])
    tsv_write(out/'ug27_candidate_identifiers.tsv', id_rows, ['identifier','sheet','row','ug27_context','experimental_context','context'])

    # Conservative auto-resolution: only identifiers explicitly on rows containing both UG27 and experimental wording.
    explicit=sorted({r['identifier'] for r in id_rows if r['ug27_context'] and r['experimental_context']}, key=str.upper)
    # It is safer to HOLD than to assume one identifier == one system when the supplement schema is unknown.
    pass_exact = len(explicit)==5
    gate_name = 'PASS_STAGE04C2_UG27_EXPERIMENTAL_FIVE_EXACT_SOURCE_IDS_FROZEN' if pass_exact else 'HOLD_STAGE04C2_UG27_EXACT_SOURCE_IDS_UNRESOLVED'
    gate={
        "gate":gate_name,
        "gate_pass":pass_exact,
        "expected_experimental_systems_n":5,
        "explicit_ug27_experimental_identifier_n":len(explicit),
        "explicit_identifiers":explicit,
        "keyword_hit_rows_n":len(hits),
        "candidate_identifier_rows_n":len(id_rows),
        "source_xlsx":str(xlsx),
        "source_xlsx_sha256":sha256(xlsx),
        "scoring_authorized":False if not pass_exact else False,
        "note":"Even a PASS here only freezes source identities. A separate sequence-fetch/QC stage is required before detector scoring."
    }
    (out/'GATE.json').write_text(json.dumps(gate, indent=2), encoding='utf-8')
    with (out/'UG27_EXACT_ID_FREEZE.tsv').open('w',encoding='utf-8') as f:
        f.write('system_index\texact_source_identifier\tstatus\n')
        for i,ident in enumerate(explicit,1): f.write(f'{i}\t{ident}\tAUTO_EXPLICIT\n')
    print(json.dumps(gate, indent=2))
    return 0 if pass_exact else 4

if __name__=='__main__': raise SystemExit(main())
