"""Recover one failed ENCORI download through its documented bindingSite API."""
import argparse
import csv
import gzip
import hashlib
import json
import sqlite3
from collect_rbp_resources import ROOT, download


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--dataset", required=True)
    args = parser.parse_args()
    asm, ds = args.assembly, args.dataset
    relative = f"raw/encori/peaks/{asm}/{ds}.api.raw.bed.gz"
    url = f"https://rnasysu.com/encori/api/bindingSite/?assembly={asm}&datasetID={ds}"
    result = download("ENCORI", url, relative, compress=True)
    if result["status"] == "failed":
        raise RuntimeError(result)
    path = ROOT / f"integrated/clean_peaks/{asm}/{ds}.bed.gz"
    path.parent.mkdir(parents=True, exist_ok=True)
    output = path.with_name(path.name + ".part")
    discarded, count, sha = [], 0, hashlib.sha256()
    with gzip.open(ROOT / relative, "rt") as source, gzip.open(output, "wt") as dest:
        for number, line in enumerate(source, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            if line.startswith("SELECT DISTINCT `DataSetId` FROM") and "#please cite:" in line:
                discarded.append((number, "upstream_SQL_diagnostic", line.rstrip())); continue
            fields = line.rstrip().split("\t")
            assert len(fields) >= 6 and int(fields[1]) >= 0 and int(fields[2]) > int(fields[1]) and fields[5] in {"+", "-", "."}, repr(line)
            assert fields[3].startswith(ds + "-"), "API returned a different experiment"
            dest.write(line); sha.update(line.encode()); count += 1
    assert count > 0
    output.replace(path)
    with (ROOT / f"integrated/{ds}_discarded_api_preface.tsv").open("w") as f:
        w=csv.writer(f, delimiter="\t");w.writerow(["line","reason","original_text"]);w.writerows(discarded)
    metadata = {"experiment": ds,"assembly": asm,"original_payload": f"raw/encori/peaks/{asm}/{ds}.bed.gz",
        "original_problem":"Provider moduleDownload.php returned PHP memory exhaustion instead of BED data",
        "api_url":url,"api_raw_file":relative,"usable_peak_file":str(path.relative_to(ROOT)),
        "valid_rows":count,"record_sha256":sha.hexdigest(),"removed_preface_lines":len(discarded)}
    overridepath = ROOT / "peak_overrides.json"
    alloverrides = json.loads(overridepath.read_text()) if overridepath.exists() else {}
    alloverrides[ds] = metadata; overridepath.write_text(json.dumps(alloverrides,indent=2))
    building = ROOT / "integrated/rbp_atlas.sqlite.building"
    if building.exists():
        db = sqlite3.connect(building)
        db.execute("UPDATE experiments SET peak_file=?,peak_sha256=?,peak_rows=?,invalid_rows=0 WHERE experiment_id=?",
            (metadata["usable_peak_file"],metadata["record_sha256"],count,ds))
        db.commit()
        statspath = ROOT / "integrated/core_stats.json"
        stats = json.loads(statspath.read_text())
        # Repairs may precede the first experiment import. In that case the
        # override is consumed by encori_experiments(), which builds its stats.
        if stats.get("ENCORI"):
            stats["ENCORI"]["peak_rows"],stats["ENCORI"]["invalid_peak_rows"] = db.execute("SELECT SUM(peak_rows),SUM(invalid_rows) FROM experiments").fetchone()
            stats["ENCORI"]["provider_error_payloads_recovered_via_API"] = len(alloverrides)
        statspath.write_text(json.dumps(stats,indent=2));db.close()
    print(json.dumps(metadata,indent=2),flush=True)


if __name__ == "__main__":
    main()
