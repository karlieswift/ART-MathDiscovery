#!/usr/bin/env python3
import argparse, csv, hashlib, json, os, subprocess, tempfile, time
from pathlib import Path
from urllib.parse import quote
from Bio import SeqIO

VERSION = "v1.0.4"
TARGET = "OK905446.1"
UA = "ART-MathDiscovery-Stage05E3A-DDBJ-v1.0.4"


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def first_nonempty(text):
    for line in text.splitlines():
        if line.strip():
            return line.strip()[:300]
    return ""


def curl_fetch(url, out_path, hdr_path, timeout=180):
    cmd = [
        "curl", "-fsSL", "--retry", "5", "--retry-all-errors", "--retry-delay", "3",
        "--connect-timeout", "20", "--max-time", str(timeout),
        "-A", UA, "-D", str(hdr_path), "-o", str(out_path), url,
    ]
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return p.returncode, p.stderr.strip()


def parse_status_and_type(hdr):
    status = ctype = ""
    if not hdr.exists():
        return status, ctype
    for ln in hdr.read_text(errors="replace").splitlines():
        if ln.startswith("HTTP/"):
            sp = ln.split()
            if len(sp) >= 2:
                status = sp[1]
        elif ln.lower().startswith("content-type:"):
            ctype = ln.split(":", 1)[1].strip()
    return status, ctype


def accession_ok(rec, requested):
    base = requested.split(".")[0]
    vals = [str(rec.id), str(rec.name)]
    vals += [str(x) for x in (rec.annotations.get("accessions", []) or [])]
    return any(v.split(".")[0] == base for v in vals if v)


def validate_genbank(path, requested):
    data = path.read_bytes() if path.exists() else b""
    text = data.decode("utf-8", "replace")
    if not text.lstrip().startswith("LOCUS"):
        return None, "does_not_start_LOCUS", data
    try:
        recs = list(SeqIO.parse(str(path), "genbank"))
    except Exception as e:
        return None, f"biopython_genbank_parse_error:{type(e).__name__}:{e}", data
    if len(recs) != 1:
        return None, f"genbank_record_count_{len(recs)}", data
    rec = recs[0]
    if len(rec.seq) <= 0:
        return None, "empty_sequence", data
    if not accession_ok(rec, requested):
        return None, f"accession_mismatch_requested_{requested}_record_{rec.id}", data
    return rec, "OK", data


def write_atomic(src, out):
    tmp = out.with_suffix(out.suffix + ".tmp")
    tmp.write_bytes(src.read_bytes())
    rs = list(SeqIO.parse(str(tmp), "genbank"))
    if len(rs) != 1 or len(rs[0].seq) <= 0:
        raise RuntimeError("generated/copied GenBank validation failed")
    os.replace(tmp, out)


def routes(sid):
    base = sid.split(".")[0]
    ids = [sid] + ([] if base == sid else [base])
    out = []
    for q in ids:
        qenc = quote(q, safe="._-")
        out.extend([
            (f"DDBJ_GETENTRY_QUERY:{q}", f"https://getentry.ddbj.nig.ac.jp/getentry?database=na&accession_number={qenc}&format=flatfile&filetype=text"),
            (f"DDBJ_GETENTRY_SMART:{q}", f"https://getentry.ddbj.nig.ac.jp/getentry/na/{qenc}?format=flatfile&filetype=text"),
            (f"DDBJ_GETENTRY_QUERY_TRACE:{q}", f"https://getentry.ddbj.nig.ac.jp/getentry?database=na&accession_number={qenc}&format=flatfile&filetype=text&trace=true"),
        ])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--diagnostics", required=True)
    ap.add_argument("--sleep", type=float, default=1.0)
    a = ap.parse_args()

    sid = TARGET
    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    manifest = Path(a.manifest); manifest.parent.mkdir(parents=True, exist_ok=True)
    diagnostics = Path(a.diagnostics); diagnostics.parent.mkdir(parents=True, exist_ok=True)

    dfields = ["sequence_id", "route", "http_status", "content_type", "bytes", "sha256", "valid_genbank", "validation_or_error", "first_nonempty_line", "curl_error", "elapsed_s"]
    mfields = ["sequence_id", "status", "bytes", "sha256", "selected_route", "provider", "source_format", "output_format"]
    drows = []
    selected = None
    last = ""

    for route, url in routes(sid):
        with tempfile.TemporaryDirectory(prefix="e3a_v104_") as td0:
            td = Path(td0); body = td / "body.txt"; hdr = td / "headers.txt"
            t0 = time.time(); rc, cerr = curl_fetch(url, body, hdr); elapsed = time.time() - t0
            status, ctype = parse_status_and_type(hdr)
            data = body.read_bytes() if body.exists() else b""
            rec = None
            if rc == 0:
                rec, verr, data = validate_genbank(body, sid)
            else:
                verr = f"curl_rc_{rc}"
            drows.append({
                "sequence_id": sid, "route": route, "http_status": status, "content_type": ctype,
                "bytes": len(data), "sha256": sha256_bytes(data) if data else "",
                "valid_genbank": int(rec is not None), "validation_or_error": verr,
                "first_nonempty_line": first_nonempty(data.decode("utf-8", "replace")) if data else "",
                "curl_error": cerr[:500], "elapsed_s": round(elapsed, 3),
            })
            if rec is not None:
                out = outdir / f"{sid}.gbff"
                write_atomic(body, out)
                ob = out.read_bytes()
                selected = (route, len(ob), sha256_bytes(ob))
                break
            last = f"{route}:{verr}:{cerr[:180]}"
        time.sleep(max(0.0, a.sleep))

    with diagnostics.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=dfields, delimiter="\t"); w.writeheader(); w.writerows(drows)

    if selected:
        route, n, hsh = selected
        mrow = {
            "sequence_id": sid, "status": "VALID", "bytes": n, "sha256": hsh,
            "selected_route": route, "provider": "DDBJ_GETENTRY_INSDC",
            "source_format": "DDBJ_GenBank_compatible_flatfile", "output_format": "GenBank_direct_validated",
        }
        print(f"[PASS] {sid} FETCHED route={route}")
        summary = {"requested": 1, "valid": 1, "failed": 0, "failed_ids": [], "provider": "DDBJ getentry / INSDC", "scientific_settings_changed": False}
    else:
        mrow = {
            "sequence_id": sid, "status": "FAILED", "bytes": 0, "sha256": "", "selected_route": "NONE",
            "provider": "DDBJ_GETENTRY_INSDC", "source_format": "AUTO", "output_format": "GenBank",
        }
        summary = {"requested": 1, "valid": 0, "failed": 1, "failed_ids": [sid], "provider": "DDBJ getentry / INSDC", "scientific_settings_changed": False}

    with manifest.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=mfields, delimiter="\t"); w.writeheader(); w.writerow(mrow)

    print(json.dumps(summary, indent=2))
    if not selected:
        raise SystemExit(f"v1.0.4 DDBJ recovery failed: {last}")


if __name__ == "__main__":
    main()
