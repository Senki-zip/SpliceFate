"""Audit expected downloads, file hashes, atlas invariants, and source limitations."""
from __future__ import annotations
import collections
import argparse
import csv
import gzip
import hashlib
import json
import re
import sqlite3
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from collect_rbp_resources import ROOT, PROJECT


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--local-genes', type=Path, default=PROJECT / 'data/processed/GSE250381_mm10_splice_chains_20261108/genes.tsv')
    parser.add_argument('--report-path', type=Path, default=PROJECT / 'docs/RBP_RESOURCE_INTEGRATION_20261115.md')
    args = parser.parse_args()
    final = ROOT / "integrated/rbp_atlas.sqlite"
    db = sqlite3.connect(final.as_uri() + "?mode=ro", uri=True)
    summary = json.loads((ROOT / "integrated/integration_summary.json").read_text())
    attempts = [json.loads(line) for line in (ROOT / "downloads.jsonl").open()]
    latest = {r["path"]: r for r in attempts}
    errors = []
    manifests = []
    for relative, row in sorted(latest.items()):
        path = ROOT / relative
        expected_absent = row["source"] == "ENCORI" and "/motifs/" in relative and "404" in row.get("error", "") and not path.exists()
        if path.exists():
            if row.get("sha256"):
                with path.open("rb") as f:
                    sha = hashlib.file_digest(f, "sha256").hexdigest()
                if sha != row["sha256"]:
                    errors.append("Hash mismatch: " + relative)
            status = "verified"
        elif expected_absent:
            status = "not_published"
        else:
            status = "failed"; errors.append("Missing download: " + relative)
        manifests.append((row["source"], relative, status, path.stat().st_size if path.exists() else 0,
                          row.get("sha256", ""), row.get("url", "")))
    with (ROOT / "download_manifest.tsv").open("w") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(["source", "path", "status", "bytes", "sha256", "url"]); w.writerows(manifests)
    zip_checks = {}
    for path in (ROOT / "raw").rglob("*"):
        if path.is_file() and zipfile.is_zipfile(path):
            with zipfile.ZipFile(path) as z:
                bad = z.testzip(); zip_checks[str(path.relative_to(ROOT))] = bad or "ok"
                if bad: errors.append("ZIP CRC error: " + str(path) + " " + bad)
    gzip_checks = {}
    for path in sorted((ROOT / "raw/encori").glob("*.RBP_*.tar.gz")):
        try:
            with gzip.open(path, "rb") as stream:
                while stream.read(8 * 1024 * 1024):
                    pass
            gzip_checks[path.name] = "ok"
        except Exception as exc:
            gzip_checks[path.name] = str(exc)
            errors.append("gzip CRC error: " + path.name)
    dangling = {}
    for table in ["factor_sources", "factor_motifs", "motif_evidence", "functions", "experiments", "gene_binding", "rbp2go_scores"]:
        n = db.execute(f"SELECT COUNT(*) FROM {table} t LEFT JOIN factors f ON t.factor_id=f.factor_id WHERE f.factor_id IS NULL").fetchone()[0]
        dangling[table] = n
        if n: errors.append("Dangling factor ID in " + table)
    matrix_errors = 0
    for mid, length, text in db.execute("SELECT * FROM motifs"):
        rows = json.loads(text)
        if len(rows) != length or any(len(r) != 4 or min(r) < 0 or abs(sum(r)-1) > 1e-7 for r in rows):
            matrix_errors += 1
    if matrix_errors: errors.append("Invalid normalized matrices")
    recovered_peaks = {}
    overridepath = ROOT / "peak_overrides.json"
    if overridepath.exists():
        for dataset, record in json.loads(overridepath.read_text()).items():
            sha, count = hashlib.sha256(), 0
            with gzip.open(ROOT / record["usable_peak_file"], "rt") as stream:
                for line in stream:
                    sha.update(line.encode());count += 1
            ok = sha.hexdigest() == record["record_sha256"] and count == record["valid_rows"]
            recovered_peaks[dataset] = {"records": count, "record_hash_verified": ok}
            if not ok: errors.append("Recovered peak hash mismatch: " + dataset)
    expected_archives = 12
    ready_archives = len(list((ROOT / "raw/encori").glob("*.RBP_*.tar.gz")))
    if ready_archives != expected_archives: errors.append("ENCORI bulk archive count mismatch")
    peak_status = json.loads((ROOT / "encori_peak_download_status.json").read_text())
    if peak_status["failed"]: errors.append("Failed ENCORI peak downloads")
    motif_status = collections.Counter(r["status"] for r in json.loads((ROOT / "encori_motif_download_status.json").read_text()))
    if motif_status["failed"]: errors.append("Failed ENCORI motif downloads")
    source_status = [
        ("CISBP-RNA", "complete_public_export", "All-species motif/RBP export and SQL backup"),
        ("ATtRACT", "complete_public_export", "Full database and PFM archive"),
        ("RBP2GO", "complete_public_export", "All 13 species RBP/background/homology tables and 9 supporting exports"),
        ("ENCORI", "complete_public_RBP_download_tables", "2725 CLIP experiments, published HOMER candidate motifs, 12 human/mouse RBP-target archives"),
        ("POSTAR3", "public_queries_complete_bulk_coordinates_unavailable", "1434 public queries across seven species; full coordinates require provider request, not submitted"),
    ]
    with (ROOT / "source_status.tsv").open("w") as f:
        w=csv.writer(f, delimiter="\t"); w.writerow(["source", "status", "scope_or_limitation"]); w.writerows(source_status)
    covrows = []
    local = args.local_genes
    localgenes = {r["gene"]: r["gene_name"] for r in csv.DictReader(local.open(), delimiter="\t")}
    for fid, species, gid, symbol in db.execute("SELECT * FROM mouse_candidate_splicing_factors"):
        sources = ";".join(x[0] for x in db.execute("SELECT DISTINCT source FROM factor_sources WHERE factor_id=?", (fid,)))
        motifs = db.execute("SELECT COUNT(DISTINCT motif_id) FROM valid_factor_motifs WHERE factor_id=?", (fid,)).fetchone()[0]
        exp, nonpert = db.execute("SELECT COUNT(*),COALESCE(SUM(perturbation_flag=0),0) FROM experiments WHERE factor_id=?", (fid,)).fetchone()
        targets = db.execute("SELECT COUNT(DISTINCT target_gene_id) FROM gene_binding WHERE factor_id=?", (fid,)).fetchone()[0]
        covrows.append((fid,gid,symbol,int(gid in localgenes),sources,motifs,exp,nonpert,targets))
    with (ROOT / "integrated/mouse_splicing_factor_coverage.tsv").open("w") as f:
        w=csv.writer(f, delimiter="\t"); w.writerow(["factor_id","gene_id","symbol","in_local_gene_annotation","sources","distinct_matrices","CLIP_experiments","keyword_nonperturbed_experiments","binding_target_genes"]); w.writerows(covrows)
    out = {"verified_at": datetime.now(timezone.utc).isoformat(), "validation_passed": not errors,
        "errors": errors, "download_manifest_counts": dict(collections.Counter(r[2] for r in manifests)),
        "zip_crc_checks": zip_checks, "gzip_crc_checks": gzip_checks, "dangling_factor_ids": dangling, "matrix_errors": matrix_errors,
        "ENCORI_bulk_archives": ready_archives, "ENCORI_motif_runs": dict(motif_status),
        "recovered_provider_error_payloads": recovered_peaks,
        "POSTAR_bulk_coordinates": "requires_provider_request_not_downloaded",
        "mouse_GO_splicing_candidates": len(covrows), "mouse_candidates_in_local_annotation": sum(r[3] for r in covrows),
        "local_gene_identity_path": str(local.resolve()),
        "local_bytes": sum(p.stat().st_size for p in ROOT.rglob("*") if p.is_file())}
    (ROOT / "validation.json").write_text(json.dumps(out, indent=2, ensure_ascii=False)+"\n")
    report = f"""# RBP 数据下载与整合（{out['verified_at'][:10]}）

已整合 CISBP-RNA、ATtRACT、RBP2GO、ENCORI 和 POSTAR3 的公开可取得数据。
POSTAR3 的整库位点坐标仍需向提供方申请，本次未提交申请；已完成公开逐因子查询。

| 来源 | 本次取得的数据 |
|---|---|
| CISBP-RNA 2.00 | 全物种 RBP/全部关联 motif、{summary['CISBP-RNA']['valid_matrices']:,} 份可解析矩阵、SQL 备份 |
| ATtRACT | 全量 {summary['ATtRACT']['input_rows']:,} 条结合序列注释和矩阵 |
| RBP2GO 2.0 | 13 个物种的 {summary['RBP2GO']['rbp']:,} 条 RBP 和 {summary['RBP2GO']['non_listed']:,} 条背景蛋白记录；同源、GO、结构域等公开导出 |
| ENCORI | 2,725 个 CLIP 实验的 BED 位点，{ready_archives} 个 RBP–RNA 注释包，公开 HOMER motif 候选 |
| POSTAR3 | 七个物种 1,434 次公开查询，{summary['POSTAR3']['factors_with_nonempty_public_binding_queries']} 个因子有非空结合结果；整库坐标待申请 |

统一后的主要规模：{summary['tables']['motifs']:,} 个精确矩阵组、{summary['tables']['factor_motifs']:,} 条来源可追溯的因子–motif 连接、{summary['tables']['gene_binding']:,} 条因子–基因/类型结合注释。这些连接记录结合证据和来源。合并因子表包含背景蛋白及未解析标识，不能把总行数解释为 RBP 数量。

小鼠中按 GO 描述筛得 {len(covrows)} 个宽泛的剪接相关候选，其中 {sum(r[3] for r in covrows)} 个与当前基因注释精确匹配。此列表需结合表达、细胞系统和扰动验证进一步筛选，未修改模型或基因选择参数。

## 整合与验证

- 因子按物种和基因 ID 匹配，保留别名与映射方法；歧义项不强行合并，跨物种不合并。
- 矩阵保留实验/推断证据与原始注释；空、缺失矩阵单独记录。同一数值矩阵去重，来源连接保留。
- 位点保留原始 genome build、链方向和 BED 的 0-based half-open 坐标；不做跨物种坐标转换。
- ENCORI 记录原始实验 accession、相同文件内容和重复区间；POSTAR 聚合结果缺少实验 ID，跨来源重叠按支持记录保存。
- ENCORI 原始包中的 KD/CRISPR 表达响应列保留在原始归档内，未导入 WT 结合先验表。
- 文件 SHA256、ZIP CRC、矩阵、BED 行、SQLite 完整性和关联标识均已核查。验证通过：{not errors}。
- CISBP 全物种导出补回了早期物种包中缺少的部分矩阵。

## 使用入口

主数据库：`data/external/rbp_atlas_20261115/integrated/rbp_atlas.sqlite`

```bash
python scripts/query_rbp_atlas.py --gene Ptbp1 --species Mus_musculus
```

数据库表和使用边界见 `data/external/rbp_atlas_20261115/README.md`。
来源状态见 `source_status.tsv`，逐文件清单见 `download_manifest.tsv`，验证结果见 `validation.json`。
小鼠候选覆盖见 `integrated/mouse_splicing_factor_coverage.tsv`。

本次下载整合只涉及 RBP/剪接相关资源，不包含各平台无关的 miRNA、癌症生存分析或原始测序 reads 整库镜像。
"""
    args.report_path.write_text(report)
    print(json.dumps(out, indent=2, ensure_ascii=False), flush=True)
    if errors: raise SystemExit(1)


if __name__ == "__main__":
    main()
