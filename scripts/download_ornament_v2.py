"""Download official oRNAment V2 exports, retaining provenance and checksums.

Re-run to resume partial files. Coordinates are preserved in the provider assembly;
this download does not perform liftover or import sites into model training.
"""
from __future__ import annotations
import concurrent.futures as cf
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import threading
import time
from urllib.parse import urlparse
import requests

ROOT = Path(__file__).resolve().parents[1] / "data/external/ornament_v2_20261115"
LOCK = threading.RLock()
STATE = {}
RECORDS = {}
START = time.time()
CHUNK = 32 * 1024 * 1024

def stamp():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

def local_path(row):
    return ROOT / "raw" / row["url"].split("/static/download/", 1)[1]

def save(status="running"):
    with LOCK:
        result = {"database": "oRNAment V2", "updated_at": stamp(), "status": status,
                  "scope": "all 15 species V2 exports and shared metadata; V1 legacy excluded",
                  "elapsed_seconds": round(time.time()-START),
                  "expected_files": len(STATE), "files": list(STATE.values())}
        for key, value in [("download_status.json", result), ("download_manifest.json", list(RECORDS.values()))]:
            p = ROOT/key
            temp = p.with_suffix(p.suffix+".tmp")
            temp.write_text(json.dumps(value, indent=2))
            temp.replace(p)

def update(row, **extra):
    with LOCK:
        STATE[row["url"]].update(extra)

def validate(path, row, expected_md5):
    official_name = local_path(row).name
    md5, sha = hashlib.md5(), hashlib.sha256()
    with path.open("rb") as f:
        prefix = f.read(128)
        if prefix.lstrip().lower().startswith((b"<!doctype", b"<html")):
            raise ValueError("HTML returned instead of data")
        if official_name.endswith(".xz") and not prefix.startswith(b"\xfd7zXZ\x00"):
            raise ValueError("Not an XZ stream")
        f.seek(0)
        for block in iter(lambda:f.read(4*1024*1024), b""):
            md5.update(block); sha.update(block)
    if row.get("bytes") is not None and path.stat().st_size != row["bytes"]:
        raise ValueError("File length differs from official HTTP Content-Length")
    actual = md5.hexdigest()
    expected = expected_md5.get(official_name)
    if expected and actual != expected:
        raise ValueError(f"Official MD5 mismatch: {actual} != {expected}")
    return {"md5": actual, "official_md5": expected,
            "official_md5_match": actual == expected if expected else None,
            "sha256": sha.hexdigest(), "bytes": path.stat().st_size}

def download(row, expected_md5):
    url = row["url"]
    path = local_path(row)
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name+".part")
    session = requests.Session()
    session.headers.update({"User-Agent":"SpliceFate-resource-download/1.0", "Accept-Encoding":"identity"})
    try:
        with LOCK:
            previous = RECORDS.get(url, {})
        if path.exists() and previous.get("sha256") and path.stat().st_size == previous.get("bytes"):
            checked = validate(path,row,expected_md5)
        else:
            if path.exists():
                checked = validate(path,row,expected_md5)
            else:
                update(row,status="downloading")
                failures = 0
                total = row.get("bytes")
                while total is None or (part.stat().st_size if part.exists() else 0) < total:
                    offset = part.stat().st_size if part.exists() else 0
                    if total is not None and offset > total:
                        raise ValueError("Partial file is larger than provider file")
                    end = min(offset+CHUNK-1,total-1) if total is not None else offset+CHUNK-1
                    headers = {"Range":f"bytes={offset}-{end}"}
                    if row.get("etag"): headers["If-Range"] = row["etag"]
                    try:
                        with session.get(url, headers=headers, stream=True, timeout=(20,45)) as r:
                            if r.status_code == 403:
                                # The public HTTP URL is also linked by the provider.
                                with session.get(url.replace("https://", "http://", 1), headers=headers,
                                                 stream=True, timeout=(20,45)) as alt:
                                    alt.raise_for_status()
                                raise RuntimeError("HTTPS endpoint returned 403; HTTP did not provide a downloadable response")
                            r.raise_for_status()
                            if r.status_code == 206:
                                m = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)",r.headers.get("Content-Range",""))
                                if not m or int(m[1]) != offset:
                                    raise ValueError("Unexpected HTTP Content-Range")
                                if total is not None and int(m[3]) != total:
                                    raise ValueError("Provider file size changed")
                                total = int(m[3])
                                response_length = int(m[2])-offset+1
                            elif r.status_code == 200 and offset == 0:
                                response_length = int(r.headers["Content-Length"]) if r.headers.get("Content-Length") else None
                                if total is not None and response_length is not None and total != response_length:
                                    raise ValueError("Provider file size changed")
                                if total is None: total=response_length
                            else:
                                raise ValueError("Provider did not honor resume request")
                            n=0
                            with part.open("ab") as f:
                                for b in r.iter_content(256*1024):
                                    if b:
                                        f.write(b); n+=len(b)
                                        update(row,downloaded_bytes=offset+n,total_bytes=total)
                            if response_length is not None and n != response_length:
                                raise ValueError("Truncated response")
                            if total is None: total=part.stat().st_size
                        failures=0
                    except Exception as exc:
                        failures+=1
                        update(row,last_error=str(exc),retry=failures)
                        if failures >= 8 or (isinstance(exc, requests.HTTPError) and exc.response.status_code in (403,404)):
                            raise
                        time.sleep(min(2**failures,20))
                update(row,status="checking",downloaded_bytes=part.stat().st_size)
                checked = validate(part,row,expected_md5)
                part.replace(path)
        record = {"url":url,"path":str(path.relative_to(ROOT)),"downloaded_at":stamp(),
                  "provider_etag":row.get("etag"), "provider_last_modified":row.get("modified"), **checked}
        with LOCK:
            RECORDS[url]=record
        update(row,status="verified",downloaded_bytes=checked["bytes"],official_md5_match=checked["official_md5_match"])
        print(json.dumps({"verified":str(path.relative_to(ROOT)),"bytes":checked["bytes"],"official_md5_match":checked["official_md5_match"]}),flush=True)
    except Exception as exc:
        update(row,status="failed",error=f"{type(exc).__name__}: {exc}")
        print(json.dumps({"failed":url,"error":str(exc)}),flush=True)
    finally:
        session.close(); save()

def main():
    rows=json.loads((ROOT/"discovery/official_files.json").read_text())
    manifest=ROOT/"download_manifest.json"
    if manifest.exists():
        for rec in json.loads(manifest.read_text()): RECORDS[rec["url"]]=rec
    expected_md5={}
    for file in (ROOT/"discovery").glob("*.md5"):
        for line in file.read_text().splitlines():
            if line.strip():
                digest,name=line.split(maxsplit=1)
                expected_md5[Path(name.strip().lstrip("*")).name]=digest
    for row in rows:
        STATE[row["url"]]={"url":row["url"],"path":str(local_path(row).relative_to(ROOT)),
                           "status":"pending","total_bytes":row.get("bytes"),"downloaded_bytes":0}
    required=sum(r.get("bytes") or 0 for r in rows)
    if shutil.disk_usage(ROOT).free < required + 10*1024**3:
        raise RuntimeError("Insufficient free space for planned exports plus 10 GiB margin")
    # Shared metadata first, then mouse; all other current-version species follow.
    def order(row):
        u=row["url"]
        return (0 if "/data/" not in u or u.endswith(".md5") else
                1 if "mus_musculus" in u else 2 if "homo_sapiens" in u else 3,
                row.get("bytes") or 0)
    rows.sort(key=order)
    save()
    print(json.dumps({"started":stamp(),"files":len(rows),"known_bytes":required,"workers":6}),flush=True)
    with cf.ThreadPoolExecutor(max_workers=6) as pool:
        pending=[pool.submit(download,row,expected_md5) for row in rows]
        while pending:
            done,pending=cf.wait(pending,timeout=20,return_when=cf.FIRST_COMPLETED)
            for future in done: future.result()
            save()
            with LOCK:
                progress={"elapsed_seconds":round(time.time()-START),
                          "verified":sum(r["status"]=="verified" for r in STATE.values()),
                          "failed":sum(r["status"]=="failed" for r in STATE.values()),
                          "downloaded_bytes":sum(r["downloaded_bytes"] for r in STATE.values())}
            print(json.dumps(progress),flush=True)
    failed=sum(r["status"]=="failed" for r in STATE.values())
    save("complete" if not failed else "complete_with_provider_limitations")
    print(json.dumps({"finished":stamp(),"failed":failed}),flush=True)

if __name__ == "__main__":
    main()
