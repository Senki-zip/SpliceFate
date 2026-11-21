"""Finish integration only after all expected public downloads are present."""
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from collect_rbp_resources import ROOT, PROJECT


def status(stage, **extra):
    (ROOT / "pipeline_status.json").write_text(json.dumps({"stage": stage,
        "updated_at": datetime.now(timezone.utc).isoformat(), **extra}, indent=2))


def main():
    while True:
        peak = ROOT / "encori_peak_download_status.json"
        motif = ROOT / "encori_motif_download_status.json"
        archives = len(list((ROOT / "raw/encori").glob("*.RBP_*.tar.gz")))
        if peak.exists() and json.loads(peak.read_text())["failed"]:
            raise RuntimeError("Peak download errors need repair")
        if peak.exists() and motif.exists() and archives == 12 and not (ROOT / "metadata_integration_inprogress").exists():
            if any(r["status"] == "failed" for r in json.loads(motif.read_text())):
                raise RuntimeError("Motif download errors need repair")
            break
        status("downloading", peak_files=len(list((ROOT / "raw/encori/peaks").rglob("*.gz"))),
               bulk_archives_complete=archives, bulk_archives_planned=12)
        time.sleep(15)
    status("integrating")
    with (ROOT / "integration_finish.log").open("w") as log:
        subprocess.run([sys.executable, "-u", "scripts/integrate_rbp_resources.py", "--phase", "finish"],
                       cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT, check=True)
    status("validating")
    with (ROOT / "validation.log").open("w") as log:
        subprocess.run([sys.executable, "-u", "scripts/validate_rbp_collection.py"],
                       cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT, check=True)
    status("complete_with_provider_limitations", validation="passed",
           provider_limitations=["POSTAR3 whole-database coordinate export requires provider request; public queries were downloaded"])


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        status("failed", error=str(exc))
        raise
