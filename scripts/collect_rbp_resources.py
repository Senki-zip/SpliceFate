"""Download public RBP resource exports; preserve provenance and resumable files."""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import gzip
import hashlib
import json
import re
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin

import requests

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "data/external/rbp_atlas_20261115"
LOCK = threading.Lock()


def record(item):
    item["recorded_at"] = datetime.now(timezone.utc).isoformat()
    with LOCK:
        with (ROOT / "downloads.jsonl").open("a") as f:
            f.write(json.dumps(item, ensure_ascii=False) + "\n")


def download(source, url, relative, compress=False, attempts=3, form=None):
    path = ROOT / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size:
        return {"source": source, "path": relative, "status": "cached"}
    part = path.with_name(path.name + ".part")
    error = ""
    for attempt in range(attempts):
        try:
            request = requests.get if form is None else requests.post
            kwargs = {} if form is None else {"data": form}
            with request(url, stream=True, timeout=(20, 90), headers={"User-Agent": "Mozilla/5.0"}, **kwargs) as response:
                response.raise_for_status()
                size = 0
                opener = gzip.open if compress else open
                with opener(part, "wb") as f:
                    for chunk in response.iter_content(1024 * 1024):
                        if chunk:
                            f.write(chunk)
                            size += len(chunk)
                if not size:
                    raise ValueError("Empty response")
                with (gzip.open(part, "rb") if compress else part.open("rb")) as f:
                    head = f.read(512).lstrip().lower()
                if head.startswith((b"<!doctype", b"<html")):
                    raise ValueError("HTML returned instead of data")
                part.rename(path)
                sha = hashlib.file_digest(path.open("rb"), "sha256").hexdigest()
                item = {"source": source, "url": url, "final_url": response.url, "path": relative,
                        "status": "downloaded", "bytes": path.stat().st_size,
                        "response_bytes": size, "sha256": sha,
                        "content_type": response.headers.get("Content-Type", ""),
                        "content_disposition": response.headers.get("Content-Disposition", ""),
                        "local_gzip": compress}
                if form is not None:
                    item["form"] = form
                record(item)
                return item
        except Exception as exc:
            error = f"{type(exc).__name__}: {exc}"
            if part.exists():
                part.unlink()
            if attempt + 1 < attempts:
                time.sleep(2 * (attempt + 1))
    item = {"source": source, "url": url, "path": relative, "status": "failed", "error": error}
    record(item)
    return item


def basic_exports():
    jobs = [("ATtRACT", "https://attract.cnic.es/attract/static/ATtRACT.zip", "raw/attract/ATtRACT.zip", False),
            ("ENCORI", "https://rnasysu.com/encori/api/ref/ENCORI_referenceData.zip", "raw/encori/ENCORI_referenceData.zip", False),
            ("CISBP-RNA", "https://cisbp-rna.ccbr.utoronto.ca/data/2_00/DataFiles/SQLDumps/SQLArchive_cisbp-rna_2_00.zip", "raw/cisbp/SQLArchive_cisbp-rna_2_00.zip", False)]
    page = (ROOT / "discovery/encori_api.html").read_text()
    for url in sorted(set(re.findall(r"https://rnasysu\.com/encori/download/[^\s<>]+\.RBP_[^\s<>]+\.tar\.gz", page))):
        jobs.append(("ENCORI", url, "raw/encori/" + url.rsplit("/", 1)[-1], False))
    with cf.ThreadPoolExecutor(max_workers=4) as pool:
        for result in pool.map(lambda args: download(*args), jobs):
            print(json.dumps(result, ensure_ascii=False), flush=True)


def cisbp_export():
    base = "https://cisbp-rna.ccbr.utoronto.ca/"
    destination = ROOT / "raw/cisbp/All_species_RBP_information_PWMs.zip"
    if destination.exists():
        print("CISBP all-species archive already downloaded", flush=True)
        return
    form = [("selSpec", "0"), ("selFamily", "0"),
            ("EntireDataset[]", "RBP_Information"), ("EntireDataset[]", "PWMs"),
            ("submit", "Download Entire Archive!")]
    response = requests.post(base + "bulk_archive.php", data=form, timeout=(20, 300))
    response.raise_for_status()
    (ROOT / "discovery/cisbp_all_export.html").write_text(response.text)
    match = re.search(r'href="(tmp/[^\"]+\.zip)"', response.text)
    if not match:
        record({"source": "CISBP-RNA", "status": "failed", "error": "No archive link in all-species response"})
        return
    print(json.dumps(download("CISBP-RNA", urljoin(base, match[1]), str(destination.relative_to(ROOT)))), flush=True)


def encori_peaks():
    jobs = []
    for path in sorted((ROOT / "raw/encori").glob("*_experiments.json")):
        assembly = path.name.split("_experiments")[0]
        entries = json.loads(path.read_text())
        if "data" not in entries:
            continue
        assert len(entries["data"]) == int(entries["recordsTotal"])
        for entry in entries["data"]:
            accession = entry["dataSetId"]
            url = "https://rnasysu.com/encori/moduleDownload.php?source=rbpbed&type=bed&value=" + assembly + ";" + accession
            jobs.append(("ENCORI", url, f"raw/encori/peaks/{assembly}/{accession}.bed.gz", True))
    failed = []
    with cf.ThreadPoolExecutor(max_workers=5) as pool:
        for n, item in enumerate(pool.map(lambda args: download(*args), jobs), 1):
            if item["status"] == "failed":
                failed.append(item)
                print(json.dumps(item), flush=True)
            if n % 50 == 0 or n == len(jobs):
                print(json.dumps({"task": "encori_peaks", "completed": n, "planned": len(jobs), "failed": len(failed)}), flush=True)
    (ROOT / "encori_peak_download_status.json").write_text(json.dumps({"planned": len(jobs), "failed": failed}, indent=2))


def postar_export():
    import ast
    jobs = []
    catalog = []
    for species in ["mouse", "human", "fly", "worm", "arabidopsis", "yeast", "zebrafish"]:
        relative = f"raw/postar3/{species}_RBP_genelist.txt.js"
        url = f"http://111.198.139.65/data/postar3/{species}_RBP_genelist.txt.js"
        item = download("POSTAR3", url, relative)
        if item["status"] == "failed":
            continue
        text = (ROOT / relative).read_text()
        names = ast.literal_eval(text.split("=", 1)[1].strip().rstrip(";"))
        for name in sorted(set(names)):
            catalog.append({"species": species, "RBP": name})
            for operation in ["RBP_info", "RBP_GO", "search_RBP"]:
                relative = f"raw/postar3/queries/{species}/{name}_{operation}.json.gz"
                form = {"species": species, "RBP_name": name, "searchID": operation}
                jobs.append((relative, form))
    (ROOT / "raw/postar3/factor_catalog.json").write_text(json.dumps(catalog, indent=2))
    failures = []
    def run(args):
        relative, form = args
        return download("POSTAR3", "http://111.198.139.65/script/ajax/RBP.php", relative, True, form=form)
    with cf.ThreadPoolExecutor(max_workers=3) as pool:
        for n, item in enumerate(pool.map(run, jobs), 1):
            if item["status"] == "failed":
                failures.append(item)
            if n % 50 == 0 or n == len(jobs):
                print(json.dumps({"task": "postar_public_queries", "completed": n, "planned": len(jobs), "failed": len(failures)}), flush=True)
    (ROOT / "postar_download_status.json").write_text(json.dumps({"public_queries": len(jobs), "failures": failures,
        "whole_database_binding_coordinates": "requires provider request; not downloaded",
        "provider_request_page": "https://wj.qq.com/s2/10477617/4f84/v"}, indent=2))


def rbp2go_export():
    state = json.loads((ROOT / "discovery/rbp2go_browser_state.json").read_text())
    identifiers = {"all-download_species_rbp", "all-download_species_non_listed", "all-download_species_homologs",
        "download_all_studies", "download_interpro_data", "Down.Table_RBDs_TXT", "Down.Table_fam_IDs_TXT",
        "Down.High_conf_RBPs_TXT", "download_peptide_data", "download_corum_complexes",
        "download_app_scripts", "download_analysis_script"}
    jobs = [("RBP2GO", r["href"], "raw/rbp2go/" + r["id"] + ".export", False)
            for r in state["links"] if r["id"] in identifiers]
    assert len(jobs) == len(identifiers), "Missing public download links"
    with cf.ThreadPoolExecutor(max_workers=3) as pool:
        for item in pool.map(lambda args: download(*args), jobs):
            print(json.dumps(item), flush=True)


def encori_motifs():
    """Public HOMER combined de novo motifs; absent runs are explicitly logged."""
    jobs = []
    for path in sorted((ROOT / "raw/encori").glob("*_experiments.json")):
        assembly = path.name.split("_experiments")[0]
        for entry in json.loads(path.read_text()).get("data", []):
            for length in ["short", "long"]:
                dataset = entry["dataSetId"]
                url = f"https://rnasysu.com/encori/motif_{length}/{assembly}/{dataset}/homerMotifs.all.motifs"
                relative = f"raw/encori/motifs/{assembly}/{dataset}_{length}.motifs"
                jobs.append((url, relative))
    results = []
    def run(job):
        url, relative = job
        # A 404 means no public motif run, rather than an incomplete download.
        result = download("ENCORI", url, relative, attempts=1)
        if result.get("status") == "failed" and "404" in result.get("error", ""):
            result["status"] = "not_published"
        return result
    with cf.ThreadPoolExecutor(max_workers=6) as pool:
        for n, item in enumerate(pool.map(run, jobs), 1):
            results.append(item)
            if n % 250 == 0 or n == len(jobs):
                from collections import Counter
                print(json.dumps({"task": "encori_motifs", "completed": n, "planned": len(jobs),
                    "statuses": dict(Counter(r["status"] for r in results))}), flush=True)
    (ROOT / "encori_motif_download_status.json").write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    ROOT.mkdir(parents=True, exist_ok=True)
    parser = argparse.ArgumentParser()
    parser.add_argument("task", choices=["basic", "cisbp", "encori-peaks", "postar", "rbp2go", "encori-motifs"])
    args = parser.parse_args()
    {"basic": basic_exports, "cisbp": cisbp_export, "encori-peaks": encori_peaks,
     "postar": postar_export, "rbp2go": rbp2go_export, "encori-motifs": encori_motifs}[args.task]()
