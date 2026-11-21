"""Build a provenance-preserving mouse RBP database and Parquet site store.

Run with .venv-rbp-db/bin/python. Stages are restartable; raw sources are read-only.
"""
from __future__ import annotations
import argparse
import collections
import csv
import gzip
import hashlib
import json
import lzma
import math
import os
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import time

import duckdb
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / 'data/external/mouse_rbp_unified_20261116'
ATLAS = PROJECT / 'data/external/rbp_atlas_20261115'
ORN = PROJECT / 'data/external/ornament_v2_20261115'
REF = PROJECT.parent / 'reference_documents'
REFERENCES = {
    'mm10': REF / 'refdata-cellranger-arc-mm10-2020-A-2.0.0',
    'mm39': REF / 'refdata-cellranger-arc-GRCm39-2024-A',
}
DBPATH = ROOT / 'mouse_rbp.duckdb'
csv.field_size_limit(100_000_000)


def log(stage, **data):
    rec = {'time': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'stage': stage, **data}
    print(json.dumps(rec, ensure_ascii=False), flush=True)
    with (ROOT / 'build_events.jsonl').open('a') as f:
        f.write(json.dumps(rec, ensure_ascii=False) + '\n')


def sqlstr(value):
    return "'" + str(value).replace("'", "''") + "'"


def sha256(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def create_rows(db, name, columns, rows):
    """Typed Arrow batches, with no pandas or SQLite extension dependency."""
    db.execute(f'CREATE OR REPLACE TABLE {name} ({columns})')
    names, types = zip(*[(r[0], r[1]) for r in db.execute(f'DESCRIBE {name}').fetchall()])
    numeric = {'BIGINT': pa.int64(), 'INTEGER': pa.int32(), 'DOUBLE': pa.float64(), 'BOOLEAN': pa.bool_()}
    schema = pa.schema([(n, numeric.get(t, pa.string())) for n, t in zip(names, types)])
    batch = []
    count = 0
    def insert():
        nonlocal count
        table = pa.Table.from_pylist([dict(zip(names, r)) for r in batch], schema=schema)
        db.register('_batch', table)
        db.execute(f'INSERT INTO {name} SELECT * FROM _batch')
        db.unregister('_batch')
        count += len(batch)
        batch.clear()
    for row in rows:
        batch.append(tuple(row))
        if len(batch) == 20000:
            insert()
    if batch:
        insert()
    return count


class Fasta:
    def __init__(self, assembly):
        self.path = REFERENCES[assembly] / 'fasta/genome.fa'
        self.file = self.path.open('rb')
        self.index = {r.split()[0]: tuple(map(int, r.split()[1:5]))
                      for r in Path(str(self.path)+'.fai').read_text().splitlines()}

    def fetch(self, chrom, start, end):
        if chrom not in self.index:
            if chrom=='chrMT' and 'chrM' in self.index:chrom='chrM'
            elif chrom.startswith('chr') and chrom[3:] in self.index:chrom=chrom[3:]
        n, offset, bases, width = self.index[chrom]
        if not 0 <= start < end <= n:
            raise ValueError((chrom, start, end, n))
        left = offset + start//bases*width + start%bases
        right = offset + (end-1)//bases*width + (end-1)%bases+1
        self.file.seek(left)
        return self.file.read(right-left).replace(b'\n', b'').replace(b'\r', b'').decode().upper()


def matrices():
    text = (ORN / 'raw/oRNAment_motif_ppms.meme').read_text()
    result = {}
    for mid, name, body in re.findall(r'MOTIF (\d+) (\S+)(.*?)(?=\*+\nMOTIF|\Z)', text, re.S):
        raw = [list(map(float, line.split())) for line in body.splitlines()
               if re.fullmatch(r'[0-9eE.+\- ]+', line.strip()) and len(line.split()) == 4]
        if len(raw) != 7 or any(not all(math.isfinite(v) and v >= 0 for v in r) or sum(r) <= 0 for r in raw):
            raise ValueError('Invalid oRNAment matrix ' + mid)
        normalized = [[round(v/sum(row), 9) for v in row] for row in raw]
        txt = json.dumps(normalized, separators=(',', ':'))
        result[int(mid)] = (mid, name, raw, normalized, 'PWM_' + hashlib.sha256(txt.encode()).hexdigest()[:24], txt)
    if len(result) != 492:
        raise ValueError('Unexpected oRNAment motif count')
    return result


def reference_stage(db):
    out = ROOT / 'reference'
    out.mkdir(exist_ok=True)
    fields = 'assembly VARCHAR,gene_id VARCHAR,transcript_id VARCHAR,feature_type VARCHAR,chromosome VARCHAR,start0 BIGINT,end0 BIGINT,strand VARCHAR,gene_symbol VARCHAR,gene_type VARCHAR,transcript_type VARCHAR,exon_number INTEGER,annotation_source VARCHAR'
    def rows():
        for assembly, directory in REFERENCES.items():
            path = directory / 'genes/genes.gtf.gz'
            with gzip.open(path, 'rt') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    c = line.rstrip('\n').split('\t')
                    if c[2] not in {'gene', 'transcript', 'exon'}:
                        continue
                    a = dict(re.findall(r'(\w+) "([^"]*)"', c[8]))
                    strip = lambda x: x.split('.')[0]
                    yield (assembly, strip(a['gene_id']), strip(a.get('transcript_id', '')), c[2], c[0],
                           int(c[3])-1, int(c[4]), c[6], a.get('gene_name', ''), a.get('gene_type', ''),
                           a.get('transcript_type', ''), int(a['exon_number']) if 'exon_number' in a else None,
                           'local_filtered_GENCODE:' + directory.name)
            log('reference_annotation', assembly=assembly)
    count = create_rows(db, 'reference_features', fields, rows())
    db.execute("CREATE OR REPLACE VIEW genes AS SELECT * FROM reference_features WHERE feature_type='gene'")
    db.execute("CREATE OR REPLACE VIEW transcripts AS SELECT * FROM reference_features WHERE feature_type='transcript'")
    db.execute("CREATE OR REPLACE VIEW exons AS SELECT * FROM reference_features WHERE feature_type='exon'")
    chromosomes = []
    refs = []
    for assembly, directory in REFERENCES.items():
        for line in (directory / 'fasta/genome.fa.fai').read_text().splitlines():
            c, n, *_ = line.split()
            chromosomes.append((assembly, c, int(n)))
        for rel in ['reference.json', 'genes/genes.gtf.gz', 'fasta/genome.fa', 'fasta/genome.fa.fai']:
            p = directory/rel
            refs.append((assembly, str(p), p.stat().st_size, sha256(p)))
    create_rows(db, 'chromosomes', 'assembly VARCHAR,chromosome VARCHAR,length BIGINT', chromosomes)
    create_rows(db, 'reference_files', 'assembly VARCHAR,path VARCHAR,bytes BIGINT,sha256 VARCHAR', refs)
    log('reference_complete', features=count)


def source_rows(path, delimiter=','):
    with path.open() as f:
        first = next(line for line in f if line.startswith('"'))
        # CORUM export uses backslash-escaped quotes in free-text comments.
        reader = csv.reader(f, delimiter=delimiter,
                            escapechar='\\' if path.name=='download_corum_complexes.export' else None)
        header = next(csv.reader([first], delimiter=delimiter))
        for row in reader:
            if row:
                if len(row) != len(header):
                    raise ValueError(f'Width mismatch {path}')
                yield dict(zip(header, row))


def catalog_stage(db):
    old = sqlite3.connect((ATLAS/'integrated/rbp_atlas.sqlite').as_uri()+'?mode=ro', uri=True)
    filt = "factor_id IN (SELECT factor_id FROM factors WHERE species='Mus_musculus')"
    tables = {
        'factors': "species='Mus_musculus'",
        'factor_sources': filt, 'factor_motifs': filt, 'functions': filt,
        'rbp2go_scores': filt, 'motif_evidence': filt,
        'experiments': "species='Mus_musculus'",
        'gene_binding': "target_species='Mus_musculus'",
        'mapping_issues': "species='Mus_musculus'",
        'motif_sources': "(source,source_motif_id) IN (SELECT fm.source,fm.source_motif_id FROM factor_motifs fm JOIN factors f USING(factor_id) WHERE f.species='Mus_musculus')",
        'motifs': "motif_id IN (SELECT ms.motif_id FROM motif_sources ms JOIN factor_motifs fm USING(source,source_motif_id) JOIN factors f USING(factor_id) WHERE f.species='Mus_musculus')",
    }
    for name, where in tables.items():
        schema = old.execute(f'PRAGMA table_info({name})').fetchall()
        typemap = {'TEXT': 'VARCHAR', 'INTEGER': 'BIGINT', 'REAL': 'DOUBLE'}
        columns = ','.join('"'+r[1]+'" '+typemap[r[2]] for r in schema)
        n = create_rows(db, name, columns, old.execute(f'SELECT * FROM {name} WHERE {where}'))
        log('catalog_import', table=name, rows=n)
    old.close()
    db.execute('CREATE OR REPLACE TABLE identity_mapping (source VARCHAR,label VARCHAR,factor_id VARCHAR,method VARCHAR,candidates_json VARCHAR)')
    refs = collections.defaultdict(set)
    for gid, symbol in db.execute('SELECT DISTINCT gene_id,gene_symbol FROM genes').fetchall():
        if symbol:
            refs[symbol.upper()].add(gid)
    factors = {r[0]: r for r in db.execute('SELECT * FROM factors').fetchall()}
    symbols = collections.defaultdict(set)
    for fid, species, gid, symbol in factors.values():
        if symbol:
            symbols[symbol.upper()].add(fid)
    for fid, symbol in db.execute('SELECT DISTINCT factor_id,gene_symbol FROM factor_sources').fetchall():
        if symbol:
            symbols[symbol.upper()].add(fid)
    def resolve(symbol):
        known = symbols.get(symbol.upper(), set())
        genes = refs.get(symbol.upper(), set())
        if len(known) == 1:
            return next(iter(known)), 'unique_existing_mouse_symbol', []
        if not known and len(genes) == 1:
            gid = next(iter(genes)); fid = 'Mus_musculus:' + gid
            if fid not in factors:
                db.execute('INSERT INTO factors VALUES (?,?,?,?)', [fid,'Mus_musculus',gid,symbol])
                factors[fid]=(fid,'Mus_musculus',gid,symbol)
            symbols[symbol.upper()].add(fid)
            return fid, 'unique_local_reference_symbol', []
        fid = 'Mus_musculus:symbol:' + symbol
        if fid not in factors:
            db.execute('INSERT INTO factors VALUES (?,?,?,?)',[fid,'Mus_musculus','symbol:'+symbol,symbol])
            factors[fid]=(fid,'Mus_musculus','symbol:'+symbol,symbol)
        return fid, 'unresolved_or_ambiguous_kept_separate', sorted(known or genes)
    motifdata = matrices()
    for mid, name, raw, norm, motif_id, txt in motifdata.values():
        if not db.execute('SELECT 1 FROM motifs WHERE motif_id=?',[motif_id]).fetchone():
            db.execute('INSERT INTO motifs VALUES (?,?,?)',[motif_id,len(norm),txt])
        detail=json.dumps({'name':name,'original_matrix':raw,'alphabet':'ACGU','coordinate_source':'oRNAment V2 Ensembl115'})
        db.execute('INSERT INTO motif_sources VALUES (?,?,?,?,?)',['oRNAmentV2',mid,motif_id,'valid',detail])
    rels = list(csv.DictReader((ORN/'raw/Motif.Info.tsv').open(), delimiter='\t'))
    create_rows(db, 'ornament_motif_relationships',
                'source_motif_id VARCHAR,motif_name VARCHAR,rbp_label VARCHAR,species VARCHAR,assay VARCHAR,motif_status VARCHAR,motif_cluster VARCHAR,rbp_clipdb VARCHAR',
                ((str(int(r['motif_id'])).zfill(3),r['motif_name'],r['rbp'],r['species'],r['assay'],r['Motif Status'],r['motif_cluster'],r['rbp_clipdb']) for r in rels))
    for r in rels:
        if r['species'] != 'Mus_musculus':
            continue
        fid, method, candidates = resolve(r['rbp'])
        mid = str(int(r['motif_id'])).zfill(3)
        db.execute('INSERT INTO identity_mapping VALUES (?,?,?,?,?)',['oRNAmentV2',r['rbp'],fid,method,json.dumps(candidates)])
        db.execute('INSERT INTO factor_sources VALUES (?,?,?,?,?,?)',[fid,'oRNAmentV2',mid+':'+r['rbp'],r['rbp'],r['Motif Status'],method])
        db.execute('INSERT INTO factor_motifs VALUES (?,?,?,?,?,?)',[fid,'oRNAmentV2',mid,r['Motif Status'],r['assay'],'10.1093/nar/gkag1165'])
        db.execute('INSERT INTO motif_evidence VALUES (?,?,?,?)',[fid,'oRNAmentV2',mid,json.dumps(r,ensure_ascii=False)])
    # Include complete mouse RBP2GO records, rather than losing domains, aliases and PPI metadata.
    folder = ATLAS/'raw/rbp2go'
    uid = {u:f for f,u in db.execute('SELECT factor_id,uniprot_id FROM rbp2go_scores').fetchall()}
    def profile_rows():
        for role,file in [('RBP_catalog','all-download_species_rbp.export'),('background_non_listed','all-download_species_non_listed.export')]:
            for r in source_rows(folder/file):
                if r.get('Species') == 'Mus musculus':
                    yield (uid.get(r['Uniprot_ID']),r['Uniprot_ID'],role,json.dumps(r,ensure_ascii=False))
    create_rows(db,'rbp2go_profiles','factor_id VARCHAR,uniprot_id VARCHAR,source_role VARCHAR,record_json VARCHAR',profile_rows())
    def support_rows():
        files=['all-download_species_homologs.export','download_corum_complexes.export',
               'download_interpro_data.export','download_peptide_data.export','download_all_studies.export',
               'Down.High_conf_RBPs_TXT.export','Down.Table_RBDs_TXT.export','Down.Table_fam_IDs_TXT.export']
        for file in files:
            count=0
            for r in source_rows(folder/file, '\t' if file.startswith('Down.') else ','):
                symbol = r.get('UniProt_ID',r.get('Uniprot_ID',''))
                members = re.split('[;, ]+',r.get('Cluster_Members',''))
                organism = r.get('Organism',r.get('Species',''))
                universal = file in ['download_interpro_data.export','Down.Table_RBDs_TXT.export','Down.Table_fam_IDs_TXT.export']
                if universal or symbol in uid or any(m in uid for m in members) or organism in {'Mouse','Mus musculus'}:
                    yield (file,uid.get(symbol),symbol,'shared_metadata' if universal else 'mouse_or_mouse_member',json.dumps(r,ensure_ascii=False))
                    count+=1
            log('rbp2go_support',file=file,rows=count)
    create_rows(db,'rbp2go_support_records','source_file VARCHAR,factor_id VARCHAR,source_protein_id VARCHAR,scope VARCHAR,record_json VARCHAR',support_rows())
    for name in ['gene_biotype','transcript_biotype']:
        path=ORN/'raw'/f'{name}.csv'
        db.execute(f'CREATE OR REPLACE TABLE ornament_{name} AS SELECT * FROM read_csv({sqlstr(path)},header=true,all_varchar=true)')
    db.execute("""CREATE OR REPLACE TABLE gene_identities AS
        SELECT DISTINCT gene_id,gene_symbol,assembly,annotation_source AS source FROM genes""")
    # Stable Ensembl IDs first; exact symbol is used only if it has one reference ID.
    db.execute("""CREATE OR REPLACE TABLE target_identity_mapping AS
        WITH by_symbol AS (SELECT upper(gene_symbol) symbol,min(gene_id) gene_id,count(DISTINCT gene_id) n
                           FROM gene_identities GROUP BY upper(gene_symbol)),
        targets AS (SELECT DISTINCT target_gene_id,target_symbol,source FROM gene_binding)
        SELECT t.*,CASE WHEN starts_with(target_gene_id,'ENSMUSG') THEN split_part(target_gene_id,'.',1)
                         WHEN s.n=1 THEN s.gene_id ELSE NULL END AS canonical_gene_id,
            CASE WHEN starts_with(target_gene_id,'ENSMUSG') THEN 'source_stable_ensembl_id'
                 WHEN s.n=1 THEN 'unique_reference_exact_symbol' ELSE 'unresolved_or_ambiguous' END mapping_method
        FROM targets t LEFT JOIN by_symbol s ON upper(t.target_symbol)=s.symbol""")
    log('catalog_complete')


def stream_csv(db, path, expression, columns=None, **copy_options):
    """xz -> FIFO -> DuckDB, without a >100 GB uncompressed CSV on disk."""
    fifo=ROOT/'_stream.csv'
    if fifo.exists():fifo.unlink()
    os.mkfifo(fifo)
    # Popen's shell opens the FIFO asynchronously; no shell interpolation of filenames.
    process=subprocess.Popen(['bash','-c','exec xz -dc -- "$1" > "$2"','bash',str(path),str(fifo)],stderr=subprocess.PIPE)
    opts="header=true,parallel=false,strict_mode=true,delim=',',quote='\"',escape='\"',nullstr=['','NA']"
    if columns:
        opts+=',columns={'+','.join(sqlstr(k)+':'+sqlstr(v) for k,v in columns.items())+'}'
    else:
        opts+=',all_varchar=true'
    try:
        db.execute(expression.format(reader=f'read_csv({sqlstr(fifo)},{opts})'))
        _,err=process.communicate(timeout=120)
        if process.returncode:
            raise RuntimeError(err.decode())
    except BaseException:
        process.terminate();process.wait();raise
    finally:
        fifo.unlink(missing_ok=True)


def coordinate_audit(db):
    matrix=matrices();fa=Fasta('mm39');counts=collections.Counter();examples=[]
    path=ORN/'raw/data/tables/mus_musculus.csv.xz'
    with lzma.open(path,'rt') as f:
        for i,r in enumerate(csv.DictReader(f)):
            if i>=3000:break
            blocks=sorted(zip(map(int,re.split('[,;]',r['genomic_start'])),map(int,re.split('[,;]',r['genomic_end']))))
            m=matrix[int(r['motif'])][2]
            lo,hi=math.prod(min(x) for x in m),math.prod(max(x) for x in m)
            for shift in [0,-1]:
                seq=''.join(fa.fetch(r['chromosome'],s+shift,e+shift+1) for s,e in blocks)
                if r['strand']=='-1':seq=seq.translate(str.maketrans('ACGT','TGCA'))[::-1]
                if len(seq)!=len(m):raise ValueError('Unexpected motif length')
                score=(math.prod(row['ACGT'.index(b)] if b in 'ACGT' else 0 for row,b in zip(m,seq))-lo)/(hi-lo)
                counts[shift]+=int(abs(score-float(r['score']))<1e-7)
            if i<3:examples.append(r)
    report={'sample_size':i,'sampling':'first 3000 rows of motif 001; wider motif and split-site validation runs after import',
            'score_reconstruction_matches_by_start_shift':dict(counts),
            'conclusion':'provider genomic coordinates are 1-based inclusive; normalize each block [start-1,end)',
            'split_coordinate_delimiters':'comma and semicolon accepted; provider sample uses comma',
            'transcriptional_order':'retain original block order; interval overlaps use each block separately',
            'examples':examples}
    (ROOT/'coordinate_audit.json').write_text(json.dumps(report,indent=2))
    if counts[-1]!=3000 or counts[0]>100:
        raise ValueError('Coordinate convention could not be verified')
    fa.file.close()
    log('coordinate_audit',matches=dict(counts))


def ornament_stage(db):
    coordinate_audit(db)
    folder=ORN/'raw/data/tables'
    stream_csv(db,folder/'mus_musculus.ids.csv.xz','CREATE OR REPLACE TABLE ornament_features AS SELECT * FROM {reader}')
    stream_csv(db,folder/'mus_musculus.structure.csv.xz','CREATE OR REPLACE TABLE rna_structure AS SELECT * FROM {reader}')
    db.execute("""INSERT INTO gene_identities SELECT DISTINCT ensembl_gene_id,gene_name,'mm39','oRNAmentV2_Ensembl115'
                  FROM ornament_features WHERE ensembl_gene_id IS NOT NULL""")
    columns={k:pa.string() for k in ['gene_id','transcript_id','feature_id','motif','chromosome','region',
                                    'genomic_start','genomic_end','age_median','age_interval','ancestor']}
    columns.update({k:pa.uint8() for k in ['gene_biotype','transcript_biotype','percentile','CLIP_support']})
    columns.update(transcript_position=pa.uint32(),score=pa.float64(),unpaired_probability=pa.float64(),
                   strand=pa.int8(),species_specific=pa.bool_())
    staging=ROOT/'sites.building'
    if staging.exists():shutil.rmtree(staging)
    query="""SELECT (row_number() OVER ()+{offset})::UBIGINT site_id,*,
       list_transform(regexp_split_to_array(genomic_start,'[,;]'),lambda x: CAST(x AS BIGINT)-1) block_starts0,
       list_transform(regexp_split_to_array(genomic_end,'[,;]'),lambda x: CAST(x AS BIGINT)) block_ends0,
       list_min(list_transform(regexp_split_to_array(genomic_start,'[,;]'),lambda x: CAST(x AS BIGINT)-1)) start0,
       list_max(list_transform(regexp_split_to_array(genomic_end,'[,;]'),lambda x: CAST(x AS BIGINT))) end0
       FROM _site_batch"""
    log('ornament_sites_started',compressed_bytes=(folder/'mus_musculus.csv.xz').stat().st_size)
    staging.mkdir()
    writers={};buffers=collections.defaultdict(list);buffer_rows=collections.Counter();offset=0
    def flush(chrom):
        if not buffers[chrom]:return
        table=pa.concat_tables(buffers[chrom])
        if chrom not in writers:
            dest=staging/('chromosome='+chrom);dest.mkdir()
            writers[chrom]=pq.ParquetWriter(dest/'sites.parquet',table.schema,compression='zstd',compression_level=3,
                                           use_dictionary=True,write_statistics=True)
        writers[chrom].write_table(table,row_group_size=131072)
        buffers[chrom].clear();buffer_rows[chrom]=0
    try:
        with (folder/'mus_musculus.csv.xz').open('rb') as raw:
            with lzma.LZMAFile(raw,'rb') as stream:
                reader=pacsv.open_csv(stream,read_options=pacsv.ReadOptions(block_size=16*1024**2,use_threads=False),
                    convert_options=pacsv.ConvertOptions(column_types=columns,null_values=['','NA'],strings_can_be_null=True,
                                                         true_values=['1'],false_values=['0']))
                for batch_id,batch in enumerate(reader,1):
                    db.register('_site_batch',pa.Table.from_batches([batch]))
                    table=db.execute(query.format(offset=offset)).to_arrow_table()
                    db.unregister('_site_batch')
                    for chrom in pc.unique(table['chromosome']).to_pylist():
                        subset=table.filter(pc.equal(table['chromosome'],chrom))
                        buffers[chrom].append(subset);buffer_rows[chrom]+=len(subset)
                        if buffer_rows[chrom]>=131072:flush(chrom)
                    offset+=len(table)
                    if batch_id%50==0:log('ornament_batch',batches=batch_id,rows=offset,compressed_bytes_read=raw.tell())
                for chrom in list(buffers):flush(chrom)
    finally:
        for writer in writers.values():writer.close()
    final=ROOT/'sites'
    if final.exists():shutil.rmtree(final)
    staging.rename(final)
    # Paths are relative to this database root; query helper sets file_search_path.
    db.execute("""CREATE OR REPLACE VIEW motif_sites AS SELECT 'oRNAmentV2'::VARCHAR source,'mm39'::VARCHAR assembly,
                  motif AS source_motif_id,* EXCLUDE(motif) FROM read_parquet('sites/**/*.parquet',hive_partitioning=true)""")
    db.execute("""CREATE OR REPLACE VIEW mouse_factor_motif_sites AS
        SELECT s.*,fm.factor_id FROM motif_sites s JOIN factor_motifs fm
          ON fm.source=s.source AND fm.source_motif_id=s.source_motif_id
        WHERE s.species_specific""")
    db.execute("""CREATE OR REPLACE VIEW motif_site_blocks AS
        SELECT source,assembly,site_id,chromosome,strand,source_motif_id,gene_id,transcript_id,feature_id,
               unnest(block_starts0) start0,unnest(block_ends0) end0 FROM motif_sites""")
    log('ornament_sites_complete',rows=db.execute('SELECT count(*) FROM motif_sites').fetchone()[0])


def peaks_stage(db):
    files=ROOT/'peaks';files.mkdir(exist_ok=True)
    raw=files/'mm10.original.bed'
    experiments=db.execute('SELECT experiment_id,peak_file FROM experiments ORDER BY experiment_id').fetchall()
    counts=collections.Counter()
    with raw.open('w') as out:
        for eid,file in experiments:
            with gzip.open(ATLAS/file,'rt') as f:
                for line in f:
                    if not line.strip() or line.lstrip().startswith(('#','track','browser')):continue
                    c=line.rstrip('\n').split('\t')
                    if len(c)!=6:raise ValueError(f'Unexpected peak columns {eid}: {len(c)}')
                    chrom,start,end,name,score,strand=c
                    n=counts[eid];counts[eid]+=1
                    # Unique original row identifier; source names and score stored separately.
                    out.write('\t'.join((chrom,start,end,eid+':'+str(n),'0',strand,name,score))+'\n')
            if len(counts)%100==0:log('peak_read',experiments=len(counts),peaks=sum(counts.values()))
    db.execute(f"""CREATE OR REPLACE TABLE clip_peaks_original AS SELECT 'mm10'::VARCHAR assembly,
        split_part(peak_id,':',1) experiment_id,* FROM read_csv({sqlstr(raw)},header=false,delim='\t',quote='',escape='',
          columns={{'chromosome':'VARCHAR','start0':'BIGINT','end0':'BIGINT','peak_id':'VARCHAR',
          'placeholder':'VARCHAR','strand':'VARCHAR','source_name':'VARCHAR','source_score':'VARCHAR'}})""")
    bed6=files/'mm10.input.bed'
    db.execute(f"COPY (SELECT chromosome,start0,end0,peak_id,0,strand FROM clip_peaks_original) TO {sqlstr(bed6)} (FORMAT CSV,HEADER false,DELIMITER '\t',QUOTE '')")
    lifted=files/'mm39.lifted.bed';unmapped=files/'mm39.unmapped.bed'
    tool=ROOT/'reference/liftOver';chain=ROOT/'reference/mm10ToMm39.over.chain.gz'
    with (files/'liftover.log').open('w') as err:
        subprocess.run([str(tool),'-minMatch=1','-minBlocks=1',str(bed6),str(chain),str(lifted),str(unmapped)],check=True,stderr=err)
    db.execute(f"""CREATE OR REPLACE TABLE clip_peaks_lifted_raw AS SELECT * FROM read_csv({sqlstr(lifted)},header=false,delim='\t',quote='',escape='',
      columns={{'chromosome':'VARCHAR','start0':'BIGINT','end0':'BIGINT','peak_id':'VARCHAR','placeholder':'VARCHAR','strand':'VARCHAR'}})""")
    reasons={};reason=''
    for line in unmapped.open():
        if line.startswith('#'):reason=line.strip('# \n')
        elif line.strip():reasons[line.split('\t')[3]]=reason
    create_rows(db,'clip_unmapped_reasons','peak_id VARCHAR,reason VARCHAR',reasons.items())
    db.execute("""CREATE OR REPLACE TABLE clip_peak_mapping AS
        SELECT o.peak_id,CASE WHEN r.peak_id IS NULL THEN 'unmapped'
            WHEN r.end0-r.start0 != o.end0-o.start0 THEN 'mapped_length_changed_excluded'
            ELSE 'mapped_full_length' END status,
            u.reason, r.chromosome,r.start0,r.end0,r.strand
        FROM clip_peaks_original o LEFT JOIN clip_peaks_lifted_raw r USING(peak_id)
        LEFT JOIN clip_unmapped_reasons u USING(peak_id)""")
    db.execute("""CREATE OR REPLACE TABLE clip_peaks_mm39 AS
        SELECT 'mm39'::VARCHAR assembly,o.experiment_id,m.chromosome,m.start0,m.end0,o.peak_id,
               m.strand,o.source_name,o.source_score FROM clip_peak_mapping m JOIN clip_peaks_original o USING(peak_id)
        WHERE m.status='mapped_full_length'""")
    db.execute("""CREATE OR REPLACE VIEW clip_peaks AS SELECT assembly,experiment_id,chromosome,start0,end0,
        peak_id,strand,source_name,source_score FROM clip_peaks_original UNION ALL SELECT * FROM clip_peaks_mm39""")
    for file in [raw,bed6]:file.unlink()
    log('peaks_complete',mapping=db.execute('SELECT status,count(*) FROM clip_peak_mapping GROUP BY status').fetchall())


def final_stage(db):
    # Ensembl, UCSC and Cell Ranger use different names for the same unplaced contigs.
    # Canonical query names match the supplied FASTA; retain the provider's original name.
    local_chromosomes=[]
    for assembly,directory in REFERENCES.items():
        for line in (directory/'fasta/genome.fa.fai').read_text().splitlines():
            c,n,*_=line.split();local_chromosomes.append((assembly,c,int(n)))
    create_rows(db,'local_reference_chromosomes','assembly VARCHAR,chromosome VARCHAR,length BIGINT',local_chromosomes)
    known={(a,c):n for a,c,n in local_chromosomes};aliases=[];official=[]
    for assembly in REFERENCES:
        path=ROOT/'reference'/f'{assembly}.chrom.sizes'
        for line in path.read_text().splitlines():
            c,n=line.split();n=int(n)
            if (assembly,c) in known and known[assembly,c]!=n:raise ValueError('Reference chromosome length conflict')
            official.append((assembly,c,n))
            known.setdefault((assembly,c),n)
    for (assembly,c),n in list(known.items()):
        aliases.append((assembly,c,c,'local_or_UCSC'))
        if c.startswith(('GL','JH','MU')):aliases.append((assembly,'chr'+c,c,'Ensembl_chr_prefix'))
        if c=='chrM':aliases.append((assembly,'chrMT',c,'mitochondrial_alias'))
    create_rows(db,'chromosome_aliases','assembly VARCHAR,source_chromosome VARCHAR,chromosome VARCHAR,mapping_method VARCHAR',aliases)
    create_rows(db,'chromosomes','assembly VARCHAR,chromosome VARCHAR,length BIGINT',((a,c,n) for (a,c),n in known.items()))
    db.execute("""CREATE OR REPLACE VIEW chromosome_reference_scope AS SELECT c.*,
        l.chromosome IS NOT NULL AS sequence_available_in_local_fasta FROM chromosomes c
        LEFT JOIN local_reference_chromosomes l USING(assembly,chromosome)""")
    db.execute("""CREATE OR REPLACE VIEW motif_sites AS
        SELECT 'oRNAmentV2'::VARCHAR source,'mm39'::VARCHAR assembly,motif AS source_motif_id,
          chromosome AS source_chromosome,
          CASE WHEN chromosome='chrMT' THEN 'chrM'
               WHEN starts_with(chromosome,'chrGL') OR starts_with(chromosome,'chrJH') THEN substr(chromosome,4)
               ELSE chromosome END AS chromosome,* EXCLUDE(motif,chromosome)
        FROM read_parquet('sites/**/*.parquet',hive_partitioning=true)""")
    # Provider tables repeat some mouse RBPs under both upper-case and mouse-case labels.
    # Preserve raw labels/measurement payloads, but never multiply factor-site rows by aliases.
    db.execute('CREATE OR REPLACE TABLE factor_motifs AS SELECT DISTINCT * FROM factor_motifs')
    db.execute('CREATE OR REPLACE TABLE identity_mapping AS SELECT DISTINCT * FROM identity_mapping')
    db.execute("""CREATE OR REPLACE VIEW mouse_factor_motif_sites AS SELECT s.*,fm.factor_id
        FROM motif_sites s JOIN (SELECT DISTINCT factor_id,source,source_motif_id FROM factor_motifs) fm
        ON fm.source=s.source AND fm.source_motif_id=s.source_motif_id WHERE s.species_specific""")
    db.execute('CREATE OR REPLACE TABLE gene_identities AS SELECT DISTINCT * FROM gene_identities')
    db.execute("""CREATE OR REPLACE TABLE target_identity_mapping AS
        WITH by_symbol AS (SELECT upper(gene_symbol) symbol,min(gene_id) gene_id,count(DISTINCT gene_id) n
                           FROM gene_identities GROUP BY upper(gene_symbol)),
        targets AS (SELECT DISTINCT target_gene_id,target_symbol,source FROM gene_binding)
        SELECT t.*,CASE WHEN starts_with(target_gene_id,'ENSMUSG') THEN split_part(target_gene_id,'.',1)
                         WHEN s.n=1 THEN s.gene_id ELSE NULL END AS canonical_gene_id,
            CASE WHEN starts_with(target_gene_id,'ENSMUSG') THEN 'source_stable_ensembl_id'
                 WHEN s.n=1 THEN 'unique_reference_exact_symbol' ELSE 'unresolved_or_ambiguous' END mapping_method
        FROM targets t LEFT JOIN by_symbol s ON upper(t.target_symbol)=s.symbol""")
    db.execute("""CREATE OR REPLACE TABLE motif_gene_loci AS SELECT gene_id,chromosome,min(start0) start0,max(end0) end0
        FROM motif_sites GROUP BY gene_id,chromosome""")
    def alias_rows():
        seen=set()
        def emit(fid,alias,source,kind):
            record=(fid,alias,source,kind)
            if alias and alias not in {'NA','None','NULL'} and record not in seen:
                seen.add(record)
                return record
        for fid,gid,symbol in db.execute('SELECT factor_id,gene_id,gene_symbol FROM factors').fetchall():
            for alias,kind in [(gid,'stable_or_unresolved_id'),(symbol,'catalog_symbol')]:
                row=emit(fid,alias,'integrated_catalog',kind)
                if row:yield row
        for fid,source,sourceid,symbol in db.execute('SELECT factor_id,source,source_id,gene_symbol FROM factor_sources').fetchall():
            for alias,kind in [(sourceid,'source_identifier'),(symbol,'source_symbol')]:
                row=emit(fid,alias,source,kind)
                if row:yield row
        for fid,uniprot,payload in db.execute('SELECT factor_id,uniprot_id,record_json FROM rbp2go_profiles').fetchall():
            if not fid:continue
            r=json.loads(payload)
            for alias in [uniprot,r.get('Entry_Name','')]+re.split(r'[;,\s]+',r.get('Alias_Names','')):
                row=emit(fid,alias,'RBP2GO','protein_alias')
                if row:yield row
    create_rows(db,'factor_aliases','factor_id VARCHAR,alias VARCHAR,source VARCHAR,alias_type VARCHAR',alias_rows())
    db.execute("""CREATE OR REPLACE VIEW valid_factor_motifs AS
      SELECT fm.*,ms.motif_id FROM factor_motifs fm JOIN motif_sources ms USING(source,source_motif_id) WHERE ms.status='valid'""")
    db.execute("""CREATE OR REPLACE VIEW canonical_binding_experiments AS
      SELECT * FROM experiments WHERE coalesce(content_duplicate_of,'')=''
      QUALIFY row_number() OVER(PARTITION BY assembly,independent_evidence_key,factor_id,assay ORDER BY experiment_id)=1""")
    db.execute("""CREATE OR REPLACE VIEW candidate_nonperturbed_binding_experiments AS
      SELECT * FROM canonical_binding_experiments WHERE perturbation_flag=0""")
    db.execute("""CREATE OR REPLACE VIEW candidate_nonperturbed_clip_peaks AS
      SELECT p.*,e.factor_id,e.cell_tissue,e.assay,e.source_accession,e.independent_evidence_key
      FROM clip_peaks p JOIN candidate_nonperturbed_binding_experiments e USING(experiment_id)""")
    db.execute("""CREATE OR REPLACE VIEW candidate_splicing_factors AS SELECT DISTINCT f.* FROM factors f
      JOIN functions fn USING(factor_id) WHERE lower(fn.label) LIKE '%splic%' OR fn.term_id='GO:0008380'""")
    db.execute("""CREATE OR REPLACE VIEW factor_catalog AS SELECT f.*,
      EXISTS(SELECT 1 FROM rbp2go_scores r WHERE r.factor_id=f.factor_id AND r.source_role='RBP_catalog') rbp2go_listed,
      EXISTS(SELECT 1 FROM factor_sources s WHERE s.factor_id=f.factor_id AND s.source!='RBP2GO') other_rbp_source,
      EXISTS(SELECT 1 FROM candidate_splicing_factors s WHERE s.factor_id=f.factor_id) broad_splicing_candidate,
      EXISTS(SELECT 1 FROM valid_factor_motifs v WHERE v.factor_id=f.factor_id AND v.source!='ENCORI'
          AND v.evidence!='mutant_binding_exclude_from_WT_default') has_curated_motif,
      EXISTS(SELECT 1 FROM experiments e WHERE e.factor_id=f.factor_id) has_mouse_clip
      FROM factors f""")
    db.execute("""CREATE OR REPLACE VIEW resolved_gene_binding AS SELECT b.*,m.canonical_gene_id,m.mapping_method
      FROM gene_binding b JOIN target_identity_mapping m ON b.source=m.source AND b.target_gene_id=m.target_gene_id AND b.target_symbol=m.target_symbol""")
    snapshots={}
    for line in (ATLAS/'downloads.jsonl').read_text().splitlines():
        item=json.loads(line)
        if item.get('status')=='downloaded' and item.get('recorded_at'):
            snapshots[item['source']]=item['recorded_at'][:10].replace('-','')
    sources=[
      ('CISBP-RNA','2.00','motif/protein identity','direct and inferred motifs; no factor-specific occupancy'),
      ('ATtRACT','download_'+snapshots.get('ATtRACT','unknown'),'motif/measurement','partly reuses CISBP/RBPDB; not independent validation'),
      ('RBP2GO','2.0','protein/function/domain/complex/homology','background proteins are not confirmed RBPs; peptides are protein coordinates'),
      ('ENCORI','download_'+snapshots.get('ENCORI','unknown'),'CLIP peaks/gene aggregates/HOMER candidates','candidate motifs are before final optimization; signed KD columns excluded'),
      ('POSTAR3','download_'+snapshots.get('POSTAR3','unknown'),'factor catalog/gene binding annotations','full coordinate export unavailable; no fabricated site coordinates'),
      ('oRNAmentV2','Ensembl115_2026-11-11','predicted sites/motif relationships/RNA structure','1-based source positions normalized; mouse mm39; CLIP support not signed regulation'),
    ]
    create_rows(db,'sources','source VARCHAR,version VARCHAR,evidence_types VARCHAR,limitations VARCHAR',sources)
    inputs=[]
    ornament_manifest = ORN/'mouse_build_input_manifest.json'
    if not ornament_manifest.exists():
        ornament_manifest = ORN/'download_manifest.json'
    for p in [ATLAS/'integrated/rbp_atlas.sqlite',ornament_manifest,ATLAS/'download_manifest.tsv',
              ROOT/'reference/liftOver',ROOT/'reference/mm10ToMm39.over.chain.gz',
              ROOT/'reference/mm10.chrom.sizes',ROOT/'reference/mm39.chrom.sizes']:
        inputs.append((str(p),p.stat().st_size,sha256(p)))
    for r in json.loads(ornament_manifest.read_text()):
        if '/mus_musculus.' in r['path'] or '/data/' not in r['path']:
            inputs.append((str(ORN/r['path']),r['bytes'],r['sha256']))
    create_rows(db,'input_assets','path VARCHAR,bytes BIGINT,sha256 VARCHAR',inputs)
    create_rows(db,'build_metadata','key VARCHAR,value VARCHAR',[
        ('schema_version','1.0'),('species','Mus_musculus'),('primary_assembly','mm39'),
        ('coordinate_convention','0-based half-open, per block; original oRNAment coordinates retained'),
        ('parquet_root','sites'),('duckdb_version',duckdb.__version__),('pyarrow_version',pa.__version__),
        ('identity_policy','mouse only; stable Ensembl IDs or unique exact species symbol; ambiguous kept separate'),
        ('liftover_policy','UCSC mm10ToMm39; minMatch=1 minBlocks=1; single result and equal interval length required'),
        ('causal_policy','no inferred regulatory sign or KO effect; motif and CLIP are separate evidence'),
        ('training_policy','no perturbation effect columns imported; nonperturbed CLIP view is only a metadata keyword screen'),
        ('cross_source_policy','shared accessions / matrices are not independent replicates'),
    ])
    db.execute('CHECKPOINT')
    log('final_complete')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=['all','references','catalog','ornament','peaks','final'],default='all')
    args=parser.parse_args()
    ROOT.mkdir(parents=True,exist_ok=True)
    db=duckdb.connect(str(DBPATH))
    db.execute("SET memory_limit='16GB'")
    db.execute('SET threads=4')
    db.execute('SET preserve_insertion_order=true')
    db.execute('SET partitioned_write_max_open_files=64')
    db.execute('SET file_search_path='+sqlstr(ROOT))
    db.execute('SET temp_directory='+sqlstr(ROOT/'_tmp'))
    stages={'references':reference_stage,'catalog':catalog_stage,'ornament':ornament_stage,'peaks':peaks_stage,'final':final_stage}
    statefile=ROOT/'build_status.json'
    state=json.loads(statefile.read_text()) if statefile.exists() else {'completed_stages':[]}
    try:
        for name,fun in stages.items():
            if args.stage!='all' and name!=args.stage:continue
            if args.stage=='all' and name in state['completed_stages']:continue
            state.pop('error',None)
            state.pop('validation_passed',None)
            state.update(status='running',current_stage=name,updated_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
            statefile.write_text(json.dumps(state,indent=2))
            log('stage_started',name=name)
            fun(db)
            if name not in state['completed_stages']:state['completed_stages'].append(name)
            db.execute('CHECKPOINT')
            state.update(current_stage=None,updated_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
            statefile.write_text(json.dumps(state,indent=2))
        state['status']='built_pending_validation' if len(state['completed_stages'])==len(stages) else 'partial'
        statefile.write_text(json.dumps(state,indent=2))
    except BaseException as exc:
        state.update(status='failed',error=str(exc))
        statefile.write_text(json.dumps(state,indent=2));raise
    finally:
        db.close()


if __name__=='__main__':main()
