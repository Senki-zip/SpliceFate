"""Build a source-traceable RBP atlas without mixing perturbation readouts into priors.

Run after collect_rbp_resources.py finishes. Source evidence types remain attached to each record.
"""
from __future__ import annotations

import collections
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import re
import sqlite3
import tarfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path

csv.field_size_limit(100_000_000)

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "data/external/rbp_atlas_20261115"
OUT = ROOT / "integrated"
SPECIES = {"hg38": "Homo_sapiens", "mm10": "Mus_musculus", "dm6": "Drosophila_melanogaster",
           "ce10": "Caenorhabditis_elegans", "sacCer3": "Saccharomyces_cerevisiae",
           "human": "Homo_sapiens", "mouse": "Mus_musculus", "fly": "Drosophila_melanogaster",
           "worm": "Caenorhabditis_elegans", "yeast": "Saccharomyces_cerevisiae",
           "arabidopsis": "Arabidopsis_thaliana", "zebrafish": "Danio_rerio"}
PERTURBATION = re.compile(r"knock.?down|knock.?out|\bKD\b|\bKO\b|\bsiRNA|\bshRNA|CRISPR|mutant|\bmut\b|over.?express|deplet|treat", re.I)


def species_name(value):
    return SPECIES.get(value, value.strip().replace(" ", "_"))


def clean(value):
    return "" if value is None or str(value).strip() in {".", "NA", "NULL", "None"} else str(value).strip()


def number(value, integer=False):
    value = clean(value)
    if not value:
        return None
    return int(float(value)) if integer else float(value)


def unique_headers(header):
    seen = collections.Counter()
    result = []
    for field in header:
        seen[field] += 1
        result.append(field if seen[field] == 1 else field + "_" + str(seen[field]))
    return result


def tab_rows(stream):
    reader = csv.reader(stream, delimiter="\t")
    header = unique_headers(next(reader))
    for row in reader:
        if row:
            if len(row) != len(header):
                raise ValueError(f"Invalid TSV width: {len(row)} != {len(header)}")
            yield dict(zip(header, row))


def rbp2go_rows(path):
    with path.open() as f:
        # Official exports begin with four attribution/date comment lines.
        first = next(line for line in f if line.lstrip().startswith('"Entry_Name"'))
        header = next(csv.reader([first]))
        for row in csv.reader(f):
            if row:
                if len(row) != len(header):
                    raise ValueError("Invalid RBP2GO CSV width")
                yield dict(zip(header, row))


def normalize_matrix(rows):
    result = []
    for row in rows:
        values = [float(v) for v in row]
        if len(values) != 4 or not all(math.isfinite(v) and v >= 0 for v in values) or sum(values) <= 0:
            raise ValueError("Invalid motif matrix")
        total = sum(values)
        result.append([round(v / total, 9) for v in values])
    if not result:
        raise ValueError("Empty motif matrix")
    return result


def fasta_matrix_blocks(text):
    header = None
    rows = []
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, rows
            header, rows = line[1:].split("\t"), []
        else:
            rows.append(line.split())
    if header is not None:
        yield header, rows


def fasta_matrices(text):
    for header, rows in fasta_matrix_blocks(text):
        yield header, normalize_matrix(rows)


def write_tsv(path, fields, rows):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for row in rows:
            writer.writerow(row)


class Atlas:
    def __init__(self, resume=False, local_genes=None):
        local = Path(local_genes) if local_genes else PROJECT / "data/processed/GSE250381_mm10_splice_chains_20261108/genes.tsv"
        # Validate identity inputs before touching a resumable database.
        with local.open() as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if not {"gene", "gene_name"}.issubset(reader.fieldnames or []):
                raise ValueError("Local gene table must contain gene and gene_name: " + str(local))
            local_rows = list(reader)
        OUT.mkdir(parents=True, exist_ok=True)
        self.path = OUT / "rbp_atlas.sqlite"
        self.temp = OUT / "rbp_atlas.sqlite.building"
        if self.temp.exists() and not resume:
            self.temp.unlink()
        self.db = sqlite3.connect(self.temp)
        if not resume:
            self.db.executescript("""
        PRAGMA journal_mode=OFF;
        PRAGMA synchronous=OFF;
        CREATE TABLE factors(factor_id TEXT PRIMARY KEY,species TEXT,gene_id TEXT,gene_symbol TEXT);
        CREATE TABLE factor_sources(factor_id TEXT,source TEXT,source_id TEXT,gene_symbol TEXT,evidence TEXT,mapping_method TEXT,
          PRIMARY KEY(factor_id,source,source_id));
        CREATE TABLE motifs(motif_id TEXT PRIMARY KEY,length INTEGER,matrix_json TEXT);
        CREATE TABLE motif_sources(source TEXT,source_motif_id TEXT,motif_id TEXT,status TEXT,detail TEXT,
          PRIMARY KEY(source,source_motif_id));
        CREATE TABLE factor_motifs(factor_id TEXT,source TEXT,source_motif_id TEXT,evidence TEXT,experiment_type TEXT,pmid TEXT,
          PRIMARY KEY(factor_id,source,source_motif_id));
        CREATE TABLE functions(factor_id TEXT,source TEXT,term_id TEXT,label TEXT,
          PRIMARY KEY(factor_id,source,term_id,label));
        CREATE TABLE experiments(experiment_id TEXT PRIMARY KEY,assembly TEXT,species TEXT,rbp_label TEXT,factor_id TEXT,
          assay TEXT,cell_tissue TEXT,treatment TEXT,source_accession TEXT,original_source TEXT,pmid TEXT,
          perturbation_flag INTEGER,independent_evidence_key TEXT,peak_file TEXT,peak_sha256 TEXT,
          peak_rows INTEGER,invalid_rows INTEGER,duplicate_rows INTEGER,content_duplicate_of TEXT);
        CREATE TABLE gene_binding(factor_id TEXT,target_species TEXT,target_gene_id TEXT,target_symbol TEXT,target_type TEXT,
          source TEXT,protocol TEXT,support_count INTEGER,reference TEXT,independent_experiment_count INTEGER,
          PRIMARY KEY(factor_id,target_species,target_gene_id,source,protocol));
        CREATE TABLE rbp2go_scores(factor_id TEXT,uniprot_id TEXT PRIMARY KEY,rbp_status TEXT,score REAL,composite_score REAL,
          listed_in_datasets INTEGER,has_rbd TEXT,source_role TEXT);
        CREATE TABLE mapping_issues(source TEXT,species TEXT,label TEXT,detail TEXT);
        """)
        self.symbols = collections.defaultdict(set)
        self.direct_aliases = collections.defaultdict(set)
        self.factor_records = {}
        self.stats = {}
        self.files_qc = []
        for r in local_rows:
            self.direct_aliases[("Mus_musculus", r["gene_name"].upper())].add(r["gene"])
        if resume:
            self.stats = json.loads((OUT / "core_stats.json").read_text())
            if self.stats.get("local_gene_identity_path", str(local.resolve())) != str(local.resolve()):
                raise ValueError("Finish phase must use the same local gene table as core")
            if self.stats.get("local_gene_identity_sha256", hashlib.sha256(local.read_bytes()).hexdigest()) != hashlib.sha256(local.read_bytes()).hexdigest():
                raise ValueError("Local gene table changed after core integration")
            for fid, species, gid, label in self.db.execute("SELECT * FROM factors"):
                self.factor_records[fid] = (species, gid, label)
                self.symbols[(species, label.upper())].add(fid)
                if gid.startswith("ENS"):
                    self.direct_aliases[(species, label.upper())].add(gid)
        self.stats["local_gene_identity_path"] = str(local.resolve())
        self.stats["local_gene_identity_sha256"] = hashlib.sha256(local.read_bytes()).hexdigest()

    def factor(self, species, label, gene_id="", source="", source_id="", evidence=""):
        species, label, gene_id = species_name(species), clean(label), clean(gene_id)
        method = "provided_gene_id"
        if not gene_id:
            candidates = self.direct_aliases.get((species, label.upper()), set())
            if len(candidates) == 1:
                gene_id = next(iter(candidates)); method = "unique_species_symbol_to_gene_id"
            else:
                known = self.symbols.get((species, label.upper()), set())
                if len(known) == 1:
                    fid = next(iter(known)); method = "unique_species_symbol"
                    gene_id = self.factor_records[fid][1]
                else:
                    gene_id = "symbol:" + label
                    method = "unresolved_symbol" if not known else "ambiguous_symbol_kept_separate"
                    if len(known) > 1:
                        self.db.execute("INSERT INTO mapping_issues VALUES(?,?,?,?)", (source, species, label, json.dumps(sorted(known))))
        fid = species + ":" + gene_id
        if fid not in self.factor_records:
            self.factor_records[fid] = (species, gene_id, label)
            self.db.execute("INSERT INTO factors VALUES(?,?,?,?)", (fid, species, gene_id, label))
        if label:
            self.symbols[(species, label.upper())].add(fid)
        if gene_id.startswith("ENS") and label:
            self.direct_aliases[(species, label.upper())].add(gene_id)
        if source:
            self.db.execute("INSERT OR IGNORE INTO factor_sources VALUES(?,?,?,?,?,?)",
                (fid, source, source_id or gene_id, label, evidence, method))
        return fid

    def motif(self, source, original_id, matrix=None, status="valid", detail=""):
        mid = ""
        if matrix is not None:
            text = json.dumps(matrix, separators=(",", ":"))
            mid = "PWM_" + hashlib.sha256(text.encode()).hexdigest()[:24]
            self.db.execute("INSERT OR IGNORE INTO motifs VALUES(?,?,?)", (mid, len(matrix), text))
        self.db.execute("INSERT OR IGNORE INTO motif_sources VALUES(?,?,?,?,?)", (source, original_id, mid, status, detail))
        return mid

    def cisbp(self):
        path = ROOT / "raw/cisbp/All_species_RBP_information_PWMs.zip"
        with zipfile.ZipFile(path) as z:
            qc = collections.Counter()
            for name in z.namelist():
                if name.endswith(".txt") and name.startswith("pwms/"):
                    original = Path(name).stem
                    lines = z.read(name).decode().splitlines()
                    try:
                        if not lines:
                            raise ValueError("Empty file")
                        assert lines[0].split() == ["Pos", "A", "C", "G", "U"]
                        matrix = normalize_matrix([r.split()[1:] for r in lines[1:] if r.strip()])
                        self.motif("CISBP-RNA", original, matrix); qc["valid_matrices"] += 1
                    except Exception as exc:
                        self.motif("CISBP-RNA", original, status="unusable", detail=str(exc)); qc["unusable_matrices"] += 1
            total = 0
            with z.open("RBP_Information_all_motifs.txt") as f:
                for r in tab_rows(io.TextIOWrapper(f)):
                    fid = self.factor(r["RBP_Species"], r["RBP_Name"], r["DBID"], "CISBP-RNA", r["RBP_ID"], r["RBP_Status"])
                    if r["Motif_ID"] != ".":
                        self.db.execute("INSERT OR IGNORE INTO factor_motifs VALUES(?,?,?,?,?,?)",
                            (fid, "CISBP-RNA", r["Motif_ID"], {"D": "direct", "I": "inferred"}.get(r["RBP_Status"], "unknown"),
                             clean(r["MSource_Type"]), clean(r["PMID"])))
                    total += 1
            self.stats["CISBP-RNA"] = {**qc, "factor_motif_input_rows": total}
        self.db.commit()

    def rbp2go(self):
        counts = collections.Counter()
        for category in ["rbp", "non_listed"]:
            path = ROOT / f"raw/rbp2go/all-download_species_{category}.export"
            for r in rbp2go_rows(path):
                species = species_name(r["Species"])
                primary = clean(r["Gene_Name"]) or r["Uniprot_ID"]
                fid = self.factor(species, primary, source="RBP2GO", source_id=r["Uniprot_ID"], evidence=r["RBP_status"])
                # Resolve an alias only if there is no competing source gene mapping.
                for alias in r["Alias_Names"].split():
                    key = (species, alias.upper())
                    candidates = self.direct_aliases.get(key, set())
                    gid = self.factor_records[fid][1]
                    if not candidates or candidates == {gid}:
                        self.symbols[key].add(fid)
                self.db.execute("INSERT OR REPLACE INTO rbp2go_scores VALUES(?,?,?,?,?,?,?,?)",
                    (fid, r["Uniprot_ID"], r["RBP_status"], number(r["RBP2GO_Score"]), number(r["Composite_Score"]),
                     number(r["Listing_Count"], integer=True), r["has_RBD"], "RBP_catalog" if category == "rbp" else "background_non_listed"))
                for aspect in ["BP", "MF", "CC"]:
                    ids, labels = r.get("GO_" + aspect + "_IDs", "").split(";"), r.get("GO_" + aspect + "_Names", "").split(";")
                    if len(ids) != len(labels):
                        self.db.execute("INSERT INTO mapping_issues VALUES(?,?,?,?)", ("RBP2GO", species, primary, "GO IDs/labels count mismatch: " + aspect))
                    for n, term in enumerate(ids):
                        if term.startswith("GO:"):
                            self.db.execute("INSERT OR IGNORE INTO functions VALUES(?,?,?,?)",
                                (fid, "RBP2GO", term, labels[n] if n < len(labels) else ""))
                counts[category] += 1
            self.db.commit()
        self.stats["RBP2GO"] = dict(counts)

    def attract(self):
        with zipfile.ZipFile(ROOT / "raw/attract/ATtRACT.zip") as z:
            for header, matrix in fasta_matrices(z.read("pwm.txt").decode()):
                self.motif("ATtRACT", header[0], matrix, detail="Matrix columns A,C,G,U (T equivalent); retained as source PFM")
            total = 0
            for r in tab_rows(io.StringIO(z.read("ATtRACT_db.txt").decode())):
                gid = r["Gene_id"] if r["Gene_id"].startswith("ENS") else ""
                fid = self.factor(r["Organism"], r["Gene_name"], gid, "ATtRACT", r["Gene_id"], "mutated=" + r["Mutated"])
                evidence = "reported_experimental_binding" if r["Mutated"].lower() == "no" else "mutant_binding_exclude_from_WT_default"
                self.db.execute("INSERT OR IGNORE INTO factor_motifs VALUES(?,?,?,?,?,?)",
                    (fid, "ATtRACT", r["Matrix_id"], evidence, r["Experiment_description"], r["Pubmed"]))
                total += 1
            self.stats["ATtRACT"] = {"input_rows": total}
        self.db.commit()

    def postar(self):
        rows = counts = errors = 0
        for entry in json.loads((ROOT / "raw/postar3/factor_catalog.json").read_text()):
            species, label = entry["species"], entry["RBP"]
            paths = {op: ROOT / f"raw/postar3/queries/{species}/{label}_{op}.json.gz" for op in ["RBP_info", "RBP_GO", "search_RBP"]}
            data = {}
            for op, path in paths.items():
                try:
                    obj = json.load(gzip.open(path, "rt")); assert obj.get("code") == 1
                    data[op] = obj.get("data", [])
                except Exception as exc:
                    data[op] = []; errors += 1
                    self.db.execute("INSERT INTO mapping_issues VALUES(?,?,?,?)", ("POSTAR3", species, label, op + ": " + str(exc)))
            info = data["RBP_info"]
            gid = info[0].get("geneID", "") if info else ""
            fid = self.factor(species, label, gid, "POSTAR3", label,
                              "public_binding_rows" if data["search_RBP"] else "catalog_only_no_public_binding_rows")
            for r in data["RBP_GO"]:
                self.db.execute("INSERT OR IGNORE INTO functions VALUES(?,?,?,?)", (fid, "POSTAR3", r["GOid"], r["GOterm"]))
            for r in data["search_RBP"]:
                self.db.execute("INSERT OR IGNORE INTO gene_binding VALUES(?,?,?,?,?,?,?,?,?,?)",
                    (fid, species_name(species), r["geneID"], r["geneName"], r["geneType"], "POSTAR3", r["protocol"],
                     int(r["num"]), str(paths["search_RBP"].relative_to(ROOT)), None))
                rows += 1
            counts += bool(data["search_RBP"])
        self.stats["POSTAR3"] = {"factor_catalog_entries": len(json.loads((ROOT / "raw/postar3/factor_catalog.json").read_text())),
            "factors_with_nonempty_public_binding_queries": counts, "gene_protocol_binding_rows_input": rows,
            "invalid_queries": errors, "whole_database_coordinate_export": "not_available_without_provider_request",
            "independence": "Aggregates have no experiment accessions; may overlap ENCORI and must not count as independent replication"}
        self.db.commit()

    def encori_experiments(self):
        existing = self.db.execute("SELECT COUNT(*) FROM experiments").fetchone()[0]
        if existing == 2725 and self.stats.get("ENCORI") and self.stats["ENCORI"].get("missing_peak_files", 0) == 0:
            print("Reusing completed ENCORI peak import", flush=True)
            return
        if existing:
            self.db.execute("DELETE FROM experiments")
        stats = collections.Counter()
        same_content = {}
        overridepath = ROOT / "peak_overrides.json"
        overrides = json.loads(overridepath.read_text()) if overridepath.exists() else {}
        for metadata in sorted((ROOT / "raw/encori").glob("*_experiments.json")):
            asm = metadata.name.split("_experiments")[0]
            obj = json.loads(metadata.read_text())
            assert len(obj["data"]) == int(obj["recordsTotal"])
            for r in obj["data"]:
                ds, rbp = r["dataSetId"], r["RBP"]
                fid = self.factor(SPECIES[asm], rbp, source="ENCORI", source_id=ds, evidence="CLIP_experiment")
                path = ROOT / f"raw/encori/peaks/{asm}/{ds}.bed.gz"
                if ds in overrides:
                    path = ROOT / overrides[ds]["usable_peak_file"]
                sha = hashlib.sha256(); n = bad = duplicates = 0
                seen = set()
                if path.exists():
                    with gzip.open(path, "rt") as f:
                        for line in f:
                            if not line.strip() or line.lstrip().startswith("#"):
                                continue
                            sha.update(line.encode())
                            fields = line.rstrip().split("\t")
                            try:
                                if len(fields) < 6 or not fields[0] or int(fields[1]) < 0 or int(fields[2]) <= int(fields[1]) or fields[5] not in {"+", "-", "."}:
                                    raise ValueError("Invalid BED6")
                                key = (fields[0], fields[1], fields[2], fields[5])
                                duplicates += key in seen; seen.add(key); n += 1
                            except Exception:
                                bad += 1
                else:
                    stats["missing_peak_files"] += 1
                sha = sha.hexdigest() if n else ""
                content_key = (asm, sha)
                duplicate_of = same_content.get(content_key, "") if sha else ""
                if sha and not duplicate_of:
                    same_content[content_key] = ds
                accessions = sorted(set(re.findall(r"(?:ENCSR\w+|GSM\d+|SRR\d+|ERR\d+|DRR\d+)", r["SubAccession"])))
                independent = "|".join(accessions) or "ENCORI:" + ds
                self.db.execute("INSERT INTO experiments VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    (ds, asm, SPECIES[asm], rbp, fid, r["SeqType"], r["CellTissue"], r["Treatment"], r["SubAccession"],
                     r["Source"], r["PubMedID"], int(bool(PERTURBATION.search(r["Treatment"]))), independent,
                     str(path.relative_to(ROOT)) if path.exists() else "", sha, n, bad, duplicates, duplicate_of))
                stats[asm + "_experiments"] += 1; stats["peak_rows"] += n; stats["invalid_peak_rows"] += bad
                stats["duplicate_interval_rows_within_experiment"] += duplicates
                stats["content_duplicate_experiments"] += bool(duplicate_of)
                if stats["hg38_experiments"] % 100 == 0:
                    self.db.commit()
        self.stats["ENCORI"] = dict(stats)
        self.stats["ENCORI"]["RBP_labels_by_assembly"] = dict(self.db.execute("SELECT assembly,COUNT(DISTINCT rbp_label) FROM experiments GROUP BY assembly"))
        self.db.commit()

    def evidence_details(self):
        """Keep every distinct measurement annotation, not just one motif link."""
        self.db.execute("CREATE TABLE IF NOT EXISTS motif_evidence(factor_id TEXT,source TEXT,source_motif_id TEXT,record_json TEXT,PRIMARY KEY(factor_id,source,source_motif_id,record_json))")
        with zipfile.ZipFile(ROOT / "raw/cisbp/All_species_RBP_information_PWMs.zip") as z:
            with z.open("RBP_Information_all_motifs.txt") as f:
                for r in tab_rows(io.TextIOWrapper(f)):
                    if r["Motif_ID"] == ".":
                        continue
                    fid = species_name(r["RBP_Species"]) + ":" + r["DBID"]
                    self.db.execute("INSERT OR IGNORE INTO motif_evidence VALUES(?,?,?,?)",
                        (fid, "CISBP-RNA", r["Motif_ID"], json.dumps(r, sort_keys=True)))
        with zipfile.ZipFile(ROOT / "raw/attract/ATtRACT.zip") as z:
            for r in tab_rows(io.StringIO(z.read("ATtRACT_db.txt").decode())):
                gid = r["Gene_id"] if r["Gene_id"].startswith("ENS") else ""
                fid = self.factor(r["Organism"], r["Gene_name"], gid)
                self.db.execute("INSERT OR IGNORE INTO motif_evidence VALUES(?,?,?,?)",
                    (fid, "ATtRACT", r["Matrix_id"], json.dumps(r, sort_keys=True)))
                if r["Mutated"].lower() == "no":
                    self.db.execute("UPDATE factor_motifs SET evidence='reported_experimental_binding' WHERE factor_id=? AND source='ATtRACT' AND source_motif_id=?",
                        (fid, r["Matrix_id"]))
        self.db.commit()

    def encori_motifs(self):
        paths = sorted((ROOT / "raw/encori/motifs").rglob("*.motifs"))
        if self.stats.get("ENCORI_motifs", {}).get("public_motif_runs") == len(paths):
            print("Reusing completed ENCORI motif import", flush=True)
            return
        experiments = {r[0]: r[1] for r in self.db.execute("SELECT experiment_id,factor_id FROM experiments")}
        counts = collections.Counter()
        for path in paths:
            dataset, mode = path.stem.rsplit("_", 1)
            for header, raw_rows in fasta_matrix_blocks(path.read_text()):
                sourceid = path.parent.name + ":" + dataset + ":" + mode + ":" + header[1]
                detail = {"file": str(path.relative_to(ROOT)), "header": header,
                    "stage": "HOMER_combined_de_novo_candidates_before_final_optimization"}
                try:
                    matrix = normalize_matrix(raw_rows)
                    self.motif("ENCORI", sourceid, matrix, detail=json.dumps(detail))
                except ValueError as exc:
                    detail["matrix_error"] = str(exc)
                    self.motif("ENCORI", sourceid, status="unusable", detail=json.dumps(detail))
                    counts["unusable_candidate_matrices"] += 1
                self.db.execute("INSERT OR IGNORE INTO factor_motifs VALUES(?,?,?,?,?,?)",
                    (experiments[dataset], "ENCORI", sourceid, "CLIP_enrichment_candidate_not_unique_factor_specificity", "HOMER_de_novo", ""))
                counts["candidate_matrices"] += 1
            counts["public_motif_runs"] += 1
        self.stats["ENCORI_motifs"] = dict(counts)
        self.db.commit()

    def encori_gene_binding(self):
        """Import binding annotations only. KD and cancer columns stay in raw archives."""
        rows = 0
        archive_stats = {}
        for path in sorted((ROOT / "raw/encori").glob("*.RBP_*.tar.gz")):
            asm = path.name.split(".")[0]
            n = 0
            with tarfile.open(path, "r:gz") as archive:
                for member in archive:
                    if not member.isfile() or not member.name.endswith(".txt") or "_cluster" in member.name:
                        continue
                    with archive.extractfile(member) as stream:
                        for r in tab_rows(io.TextIOWrapper(stream)):
                            if not r.get("RBP") or not r.get("geneID"):
                                continue
                            fid = self.factor(SPECIES[asm], r["RBP"], source="ENCORI", source_id="binding:" + r["RBP"], evidence="CLIP_gene_aggregate")
                            self.db.execute("INSERT OR IGNORE INTO gene_binding VALUES(?,?,?,?,?,?,?,?,?,?)",
                                (fid, SPECIES[asm], r["geneID"], r["geneName"], r["geneType"], "ENCORI", "CLIP_aggregate:" + path.name.split(".RBP_")[1].split(".")[0],
                                 int(r.get("clipIDnum", 0)), path.name + ":" + member.name, int(r.get("clipExpNum", 0))))
                            n += 1
            rows += n; archive_stats[path.name] = n; self.db.commit()
        self.stats["ENCORI_gene_binding"] = {"input_rows": rows, "archives": archive_stats,
            "perturbation_columns_excluded": ["HepG2_sh", "K562_sh", "HepG2_crispr", "K562_crispr", "pancancerNum"],
            "coordinate_policy": "gene aggregates are not peaks; exact per-experiment BED files indexed separately"}

    def reconcile_identities(self):
        """Resolve symbols for which a later source supplied one unique gene ID."""
        canonical = collections.defaultdict(set)
        for fid, (species, gid, label) in self.factor_records.items():
            if not gid.startswith("symbol:"):
                canonical[(species, label.upper())].add(fid)
        mappings = []
        for fid, (species, gid, label) in list(self.factor_records.items()):
            candidates = canonical.get((species, label.upper()), set())
            if gid.startswith("symbol:") and len(candidates) == 1:
                mappings.append((fid, next(iter(candidates))))
        for old, new in mappings:
            for table in ["factor_sources", "factor_motifs", "motif_evidence", "functions", "gene_binding"]:
                fields = [r[1] for r in self.db.execute(f"PRAGMA table_info({table})")]
                other = [field for field in fields if field != "factor_id"]
                names = ",".join(other)
                self.db.execute(f"INSERT OR IGNORE INTO {table}(factor_id,{names}) SELECT ?,{names} FROM {table} WHERE factor_id=?", (new, old))
                self.db.execute(f"DELETE FROM {table} WHERE factor_id=?", (old,))
            for table in ["rbp2go_scores", "experiments"]:
                self.db.execute(f"UPDATE {table} SET factor_id=? WHERE factor_id=?", (new, old))
            self.db.execute("DELETE FROM factors WHERE factor_id=?", (old,))
            self.factor_records.pop(old)
        write_tsv(OUT / "late_identity_merges.tsv", ["original_factor_id", "canonical_factor_id"], mappings)
        self.stats["identity_reconciliation"] = {"late_unique_symbol_merges": len(mappings),
            "policy": "Same species, case-insensitive exact symbol, unique explicit gene ID; no homology merge"}
        self.db.commit()

    def finalize(self):
        self.db.executescript("""
        CREATE INDEX factors_species_symbol ON factors(species,gene_symbol);
        CREATE INDEX factor_motifs_factor ON factor_motifs(factor_id);
        CREATE INDEX motif_sources_matrix ON motif_sources(motif_id);
        CREATE INDEX functions_factor ON functions(factor_id);
        CREATE INDEX gene_binding_target ON gene_binding(target_species,target_gene_id);
        CREATE INDEX gene_binding_factor ON gene_binding(factor_id);
        CREATE INDEX experiment_factor ON experiments(factor_id);
        CREATE INDEX experiment_accession ON experiments(independent_evidence_key);
        CREATE VIEW valid_factor_motifs AS
          SELECT fm.*,ms.motif_id FROM factor_motifs fm JOIN motif_sources ms
          ON fm.source=ms.source AND fm.source_motif_id=ms.source_motif_id WHERE ms.status='valid';
        CREATE VIEW mouse_candidate_splicing_factors AS
          SELECT DISTINCT f.* FROM factors f JOIN functions fn ON f.factor_id=fn.factor_id
          WHERE f.species='Mus_musculus' AND (lower(fn.label) LIKE '%splic%' OR fn.term_id='GO:0008380');
        CREATE VIEW experiment_dedup_groups AS
          SELECT assembly,independent_evidence_key,COUNT(*) AS database_entries,
          GROUP_CONCAT(experiment_id) AS encori_ids FROM experiments GROUP BY assembly,independent_evidence_key;
        CREATE VIEW canonical_binding_experiments AS
          SELECT e.* FROM experiments e WHERE e.content_duplicate_of='' AND NOT EXISTS
          (SELECT 1 FROM experiments b WHERE b.assembly=e.assembly
           AND b.independent_evidence_key=e.independent_evidence_key AND b.factor_id=e.factor_id
           AND b.assay=e.assay AND b.content_duplicate_of='' AND b.experiment_id<e.experiment_id);
        CREATE VIEW candidate_nonperturbed_binding_experiments AS
          SELECT * FROM canonical_binding_experiments WHERE perturbation_flag=0;
        """)
        # Missing motif references are explicit rows rather than silently dropped.
        missing = list(self.db.execute("SELECT DISTINCT fm.source,fm.source_motif_id FROM factor_motifs fm LEFT JOIN motif_sources ms ON fm.source=ms.source AND fm.source_motif_id=ms.source_motif_id WHERE ms.source IS NULL"))
        for source, mid in missing:
            self.motif(source, mid, status="missing_from_export")
        for table in ["factors", "factor_sources", "motifs", "motif_sources", "factor_motifs", "motif_evidence", "functions", "experiments", "gene_binding", "mapping_issues", "rbp2go_scores"]:
            self.stats.setdefault("tables", {})[table] = self.db.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        checks = {"sqlite_integrity": self.db.execute("PRAGMA integrity_check").fetchone()[0]}
        checks["dangling_factor_motif_links"] = self.db.execute("SELECT COUNT(*) FROM factor_motifs fm LEFT JOIN factors f ON fm.factor_id=f.factor_id WHERE f.factor_id IS NULL").fetchone()[0]
        checks["invalid_peak_rows"] = self.db.execute("SELECT COALESCE(SUM(invalid_rows),0) FROM experiments").fetchone()[0]
        checks["missing_peak_files"] = self.db.execute("SELECT COUNT(*) FROM experiments WHERE peak_file=''").fetchone()[0]
        checks["perturbation_fields_absent_from_gene_binding"] = not any(re.search("crispr|_sh|logfc|delta_psi", r[1], re.I) for r in self.db.execute("PRAGMA table_info(gene_binding)"))
        if checks["sqlite_integrity"] != "ok" or checks["dangling_factor_motif_links"] or checks["invalid_peak_rows"] or checks["missing_peak_files"]:
            raise ValueError("Atlas validation failed: " + str(checks))
        self.stats["checks"] = checks
        for name, query in {
            "factors.tsv.gz": "SELECT * FROM factors",
            "factor_sources.tsv.gz": "SELECT * FROM factor_sources",
            "factor_motifs.tsv.gz": "SELECT * FROM valid_factor_motifs",
            "experiments.tsv": "SELECT * FROM experiments",
            "motif_qc.tsv.gz": "SELECT * FROM motif_sources",
            "mapping_issues.tsv": "SELECT * FROM mapping_issues",
            "mouse_candidate_splicing_factors.tsv": "SELECT * FROM mouse_candidate_splicing_factors",
        }.items():
            cur = self.db.execute(query); write_tsv(OUT / name, [c[0] for c in cur.description], cur)
        self.stats["created_at"] = datetime.now(timezone.utc).isoformat()
        self.stats["scope"] = "All publicly exported species from five RBP resources; POSTAR coordinate bulk export unavailable without request"
        self.db.commit(); self.db.close(); self.temp.replace(self.path)
        (OUT / "integration_summary.json").write_text(json.dumps(self.stats, indent=2, ensure_ascii=False) + "\n")
        print(json.dumps(self.stats, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["all", "core", "finish"], default="all")
    parser.add_argument("--local-genes", type=Path, help="TSV with gene and gene_name used for mouse identity mapping")
    args = parser.parse_args()
    atlas = Atlas(resume=args.phase == "finish", local_genes=args.local_genes)
    core = ["cisbp", "rbp2go", "attract", "postar"]
    finish = ["evidence_details", "encori_experiments", "encori_motifs", "encori_gene_binding", "reconcile_identities", "finalize"]
    for operation in core if args.phase == "core" else finish if args.phase == "finish" else core + finish:
        print("Integrating " + operation, flush=True)
        getattr(atlas, operation)()
    if args.phase == "core":
        atlas.db.commit(); atlas.db.close()
        (OUT / "core_stats.json").write_text(json.dumps(atlas.stats, indent=2))
        print("Core atlas ready; peak downloads can now be finalized separately", flush=True)
