#!/usr/bin/env python3
import argparse, csv, hashlib, json, os, re, shutil, subprocess, sys, tempfile, time
from pathlib import Path
from urllib.parse import quote
from Bio import SeqIO

VERSION = "v1.0.2"
BASE_URL = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ena_sequence&id={id}&format=embl&style=raw"


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def first_nonempty(text):
    for line in text.splitlines():
        if line.strip():
            return line.strip()[:300]
    return ""


def curl_fetch(url, out_path, hdr_path, timeout=90):
    cmd = [
        "curl", "-fL", "--retry", "2", "--retry-delay", "2",
        "--connect-timeout", "20", "--max-time", str(timeout),
        "-A", "ART-MathDiscovery-Stage05E3A-ENA-v1.0.2",
        "-D", str(hdr_path), "-o", str(out_path), url,
    ]
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return p.returncode, p.stderr.strip()


def parse_status_and_type(hdr_path):
    status = ""
    ctype = ""
    if not hdr_path.exists():
        return status, ctype
    for ln in hdr_path.read_text(errors="replace").splitlines():
        if ln.startswith("HTTP/"):
            parts = ln.split()
            if len(parts) >= 2:
                status = parts[1]
        if ln.lower().startswith("content-type:"):
            ctype = ln.split(":", 1)[1].strip()
    return status, ctype


def validate_embl(path, requested):
    try:
        data = path.read_bytes()
    except Exception as e:
        return None, f"read_error:{e}", b""
    text = data.decode("utf-8", "replace")
    if not text.lstrip().startswith("ID"):
        return None, "does_not_start_EMBL_ID", data
    try:
        recs = list(SeqIO.parse(str(path), "embl"))
    except Exception as e:
        return None, f"biopython_embl_parse_error:{type(e).__name__}:{e}", data
    if len(recs) != 1:
        return None, f"embl_record_count_{len(recs)}", data
    rec = recs[0]
    if len(rec.seq) <= 0:
        return None, "empty_sequence", data
    req_base = requested.split(".")[0]
    rec_base = str(rec.id).split(".")[0]
    if req_base != rec_base:
        # Some flatfiles expose accession via annotations rather than rec.id; accept only if requested base is in accessions.
        accs = [str(x).split(".")[0] for x in rec.annotations.get("accessions", [])]
        if req_base not in accs:
            return None, f"accession_mismatch_requested_{requested}_record_{rec.id}", data
    return rec, "OK", data


def write_genbank(rec, out_path):
    if "molecule_type" not in rec.annotations:
        rec.annotations["molecule_type"] = "DNA"
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()
    SeqIO.write(rec, str(tmp), "genbank")
    # Re-parse generated GenBank before atomic move.
    rs = list(SeqIO.parse(str(tmp), "genbank"))
    if len(rs) != 1 or len(rs[0].seq) <= 0:
        raise RuntimeError("generated GenBank validation failed")
    os.replace(tmp, out_path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--diagnostics", required=True)
    ap.add_argument("--sleep", type=float, default=0.35)
    a = ap.parse_args()

    ids = [x.strip() for x in Path(a.ids).read_text().splitlines() if x.strip()]
    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    manifest = Path(a.manifest); manifest.parent.mkdir(parents=True, exist_ok=True)
    diagnostics = Path(a.diagnostics); diagnostics.parent.mkdir(parents=True, exist_ok=True)

    dfields = ["sequence_id","route","http_status","content_type","bytes","sha256","valid_embl","validation_or_error","first_nonempty_line","curl_error","elapsed_s"]
    mfields = ["sequence_id","status","bytes","sha256","selected_route","provider","source_format","output_format"]
    drows=[]; mrows=[]; failed=[]

    for i, sid in enumerate(ids, 1):
        selected = None
        attempts = [(sid, f"ENA_DBFETCH_EMBL:{sid}")]
        base = sid.split(".")[0]
        if base != sid:
            attempts.append((base, f"ENA_DBFETCH_EMBL:{base}"))
        last_err = ""
        for qid, route in attempts:
            with tempfile.TemporaryDirectory(prefix="e3a_ena_") as td:
                td = Path(td)
                body = td / "body.embl"
                hdr = td / "headers.txt"
                url = BASE_URL.format(id=quote(qid, safe="._-"))
                t0=time.time()
                rc, cerr = curl_fetch(url, body, hdr)
                elapsed=time.time()-t0
                status, ctype = parse_status_and_type(hdr)
                data = body.read_bytes() if body.exists() else b""
                rec = None; verr = ""
                if rc == 0:
                    rec, verr, data = validate_embl(body, sid)
                else:
                    verr = f"curl_rc_{rc}"
                drows.append({
                    "sequence_id":sid,"route":route,"http_status":status,"content_type":ctype,
                    "bytes":len(data),"sha256":sha256_bytes(data) if data else "","valid_embl":int(rec is not None),
                    "validation_or_error":verr,"first_nonempty_line":first_nonempty(data.decode("utf-8","replace")) if data else "",
                    "curl_error":cerr[:500],"elapsed_s":round(elapsed,3)
                })
                if rec is not None:
                    out = outdir / f"{sid}.gbff"
                    write_genbank(rec, out)
                    ob = out.read_bytes()
                    selected = (route, len(ob), sha256_bytes(ob))
                    break
                last_err = f"{route}:{verr}:{cerr[:160]}"
            time.sleep(max(0.0,a.sleep))
        if selected:
            route, n, h = selected
            mrows.append({"sequence_id":sid,"status":"VALID","bytes":n,"sha256":h,"selected_route":route,
                          "provider":"EMBL-EBI_ENA","source_format":"EMBL","output_format":"GenBank_converted_by_Biopython"})
            print(f"[{i}/{len(ids)}] {sid} FETCHED route={route}", flush=True)
        else:
            failed.append((sid,last_err))
            mrows.append({"sequence_id":sid,"status":"FAILED","bytes":0,"sha256":"","selected_route":"NONE",
                          "provider":"EMBL-EBI_ENA","source_format":"EMBL","output_format":"GenBank_converted_by_Biopython"})
            print(f"[{i}/{len(ids)}] {sid} FAILED route=NONE", flush=True)

    with diagnostics.open("w", newline="") as h:
        w=csv.DictWriter(h, fieldnames=dfields, delimiter="\t"); w.writeheader(); w.writerows(drows)
    with manifest.open("w", newline="") as h:
        w=csv.DictWriter(h, fieldnames=mfields, delimiter="\t"); w.writeheader(); w.writerows(mrows)
    summary={
        "requested":len(ids),"valid":len(ids)-len(failed),"failed":len(failed),"failed_ids":[x[0] for x in failed],
        "transport_hotfix":"v1.0.2: switch blocked NCBI retrieval to EMBL-EBI ENA Dbfetch raw EMBL; validate and convert to GenBank before unchanged Stage05E3A analysis",
        "provider":"EMBL-EBI ENA / INSDC mirrored sequence records",
        "scientific_settings_changed":False
    }
    print(json.dumps(summary,indent=2))
    if failed:
        raise SystemExit(f"ENA retrieval incomplete after versioned/versionless routes: {failed[:5]}")

if __name__ == "__main__":
    main()
