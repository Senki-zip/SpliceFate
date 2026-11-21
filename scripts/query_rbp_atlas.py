"""Small read-only interface to the integrated RBP evidence atlas."""
import argparse
import json
import sqlite3
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene", required=True, help="Gene symbol, gene ID or source protein ID")
    parser.add_argument("--species", default="Mus_musculus")
    args = parser.parse_args()
    path = Path(__file__).resolve().parents[1] / "data/external/rbp_atlas_20261115/integrated/rbp_atlas.sqlite"
    db = sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)
    db.row_factory = sqlite3.Row
    rows = db.execute("""SELECT DISTINCT f.* FROM factors f LEFT JOIN factor_sources fs ON fs.factor_id=f.factor_id
        WHERE f.species=? AND (f.gene_symbol=? COLLATE NOCASE OR f.gene_id=? OR fs.gene_symbol=? COLLATE NOCASE OR fs.source_id=?)""",
        (args.species.replace(" ", "_"), args.gene, args.gene, args.gene, args.gene)).fetchall()
    result = []
    for row in rows:
        fid = row["factor_id"]
        item = dict(row)
        item["sources"] = [dict(r) for r in db.execute("SELECT * FROM factor_sources WHERE factor_id=?", (fid,))]
        item["motif_evidence_summary"] = [dict(r) for r in db.execute("SELECT source,evidence,COUNT(*) AS links,COUNT(DISTINCT motif_id) AS distinct_matrices FROM valid_factor_motifs WHERE factor_id=? GROUP BY source,evidence", (fid,))]
        item["CLIP_experiments"] = [dict(r) for r in db.execute("SELECT experiment_id,assembly,assay,cell_tissue,treatment,source_accession,peak_rows,perturbation_flag,content_duplicate_of FROM experiments WHERE factor_id=?", (fid,))]
        item["gene_binding_summary"] = [dict(r) for r in db.execute("SELECT source,COUNT(DISTINCT target_gene_id) AS target_genes FROM gene_binding WHERE factor_id=? GROUP BY source", (fid,))]
        item["splicing_function_terms"] = [dict(r) for r in db.execute("SELECT * FROM functions WHERE factor_id=? AND lower(label) LIKE '%splic%'", (fid,))]
        result.append(item)
    print(json.dumps({"matches": result, "interpretation": "Binding and motif evidence; no signed KO effects inferred. Source overlap is not independent replication."}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
