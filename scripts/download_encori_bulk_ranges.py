"""Resume ENCORI public bulk archives with verified byte ranges."""
from __future__ import annotations
import concurrent.futures as cf
import hashlib
import json
import re
import shutil
import time
from pathlib import Path
import requests
from collect_rbp_resources import ROOT, record


def chunk(task):
    url, start, end, path, total = task
    expected = end - start + 1
    if path.exists() and path.stat().st_size == expected:
        return path
    for attempt in range(4):
        try:
            current = path.stat().st_size if path.exists() else 0
            if current == expected:
                return path
            if current > expected:
                path.unlink(); current = 0
            left = start + current
            started = time.monotonic()
            with requests.get(url, headers={"Range": f"bytes={left}-{end}"}, stream=True, timeout=(15, 30)) as r:
                r.raise_for_status()
                if r.status_code != 206 or r.headers.get("Content-Range") != f"bytes {left}-{end}/{total}":
                    raise ValueError("Server did not honor the requested byte range")
                with path.open("ab") as f:
                    for block in r.iter_content(64 * 1024):
                        f.write(block)
                        if time.monotonic() - started > 120:
                            raise TimeoutError("Resume slow connection with a fresh range request")
            if path.stat().st_size != expected:
                raise ValueError("Byte range length mismatch")
            return path
        except Exception:
            if attempt == 3:
                raise
            time.sleep(attempt + 1)


def main():
    page = (ROOT / "discovery/encori_api.html").read_text()
    urls = sorted(set(re.findall(r"https://rnasysu\.com/encori/download/[^\s<>]+\.RBP_[^\s<>]+\.tar\.gz", page)))
    archives, jobs = [], []
    for url in urls:
        path = ROOT / "raw/encori" / url.rsplit("/", 1)[-1]
        if path.exists():
            continue
        h = requests.head(url, timeout=20); h.raise_for_status()
        total = int(h.headers["Content-Length"])
        prefix = path.with_name(path.name + ".part")
        start = prefix.stat().st_size if prefix.exists() else 0
        assert start <= total
        pieces = []
        for left in range(start, total, 8 * 1024 * 1024):
            right = min(total - 1, left + 8 * 1024 * 1024 - 1)
            piece = path.with_name(path.name + f".range_{left}_{right}")
            task = (url, left, right, piece, total)
            jobs.append(task); pieces.append(piece)
        archives.append((url, path, prefix, pieces, total))
    print(json.dumps({"archives": len(archives), "range_jobs": len(jobs), "remaining_bytes": sum(e-s+1 for _,s,e,_,_ in jobs)}), flush=True)
    with cf.ThreadPoolExecutor(max_workers=16) as pool:
        for n, _ in enumerate(pool.map(chunk, jobs), 1):
            print(json.dumps({"range_completed": n, "total": len(jobs)}), flush=True)
    for url, path, prefix, pieces, total in archives:
        assembled = path.with_name(path.name + ".assembling")
        with assembled.open("wb") as dest:
            for part in ([prefix] if prefix.exists() else []) + pieces:
                with part.open("rb") as source:
                    shutil.copyfileobj(source, dest, 1024 * 1024)
        assert assembled.stat().st_size == total
        assert assembled.open("rb").read(2) == b"\x1f\x8b"
        assembled.replace(path)
        sha = hashlib.file_digest(path.open("rb"), "sha256").hexdigest()
        record({"source": "ENCORI", "url": url, "path": str(path.relative_to(ROOT)), "status": "downloaded",
                "bytes": total, "sha256": sha, "method": "verified_HTTP_ranges_with_existing_prefix"})
        for part in ([prefix] if prefix.exists() else []) + pieces:
            part.unlink()
        print(json.dumps({"archive_ready": path.name, "bytes": total}), flush=True)
    (ROOT / "encori_bulk_download_status.json").write_text(json.dumps({"expected_archives": len(urls), "complete_archives": sum((ROOT/'raw/encori'/u.rsplit('/',1)[-1]).exists() for u in urls)}, indent=2))


if __name__ == "__main__":
    main()
