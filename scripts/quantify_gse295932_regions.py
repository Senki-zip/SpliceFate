"""Quantify E/I region incidence directly from author integer-cluster structures."""
from __future__ import annotations
import argparse
import collections
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
import duckdb
import numpy as np
import pandas as pd
from scipy import sparse

from splicefate.io import digest, dump
from splicefate.regions import reference_regions, classify_blocks
from convert_gse295932 import parse_blocks, import_source, condition_samples, condition_output

PROJECT = Path(__file__).resolve().parents[1]


def quantify(source, converted, gtf, output, condition='both'):
    output.mkdir(parents=True, exist_ok=True)
    if (output/'quantification.json').exists():
        raise FileExistsError('Region quantification already exists: '+str(output))
    started = time.time()
    cells = pd.read_csv(converted/'genes/cells.tsv', sep='\t', dtype={'cell':str})
    genes = pd.read_csv(converted/'genes/genes.tsv', sep='\t', dtype={'gene':str})
    cell_index = {cell:i for i, cell in enumerate(cells.cell)}
    G = sparse.load_npz(converted/'genes/gene_counts.npz').tocsr()
    library = np.asarray(G.sum(1)).ravel().astype(np.int64)
    db = duckdb.connect()
    db.execute('SET threads=4')
    db.execute('CREATE TABLE observations(sample VARCHAR,cell VARCHAR,gene VARCHAR,isoform VARCHAR,count BIGINT)')
    sources = []
    for gsm, name, _ in condition_samples(condition):
        path = source/'WT/Nontarget_controls'/gsm/name
        print('Reading', gsm, flush=True)
        import_source(db, path, gsm)
        sources.append(path)
    print('Reading reference regions', flush=True)
    reference = reference_regions(gtf)
    # Match the integer conversion on its full source axis before classifying.
    gene_index = {g:j for j,g in enumerate(genes.gene)}
    values = db.execute('SELECT cell,gene,sum(count) FROM observations GROUP BY cell,gene').fetchall()
    imported = sparse.coo_matrix(([int(v) for _,_,v in values],
        ([cell_index[cell] for cell,_,_ in values], [gene_index[gene] for _,gene,_ in values])),
        shape=G.shape, dtype=np.int64).tocsr()
    if (imported != G).nnz:
        raise ValueError('Author tables differ from converted integer G input')
    db.execute('CREATE INDEX gene_index ON observations(gene)')
    triples = {key:[[],[],[]] for key in ['exon','intron','boundary','retained','intron_without_full_retention']}
    metrics = collections.Counter()
    audit = []
    for j, gene in enumerate(genes.gene):
        if gene not in reference:
            continue
        ref = reference[gene]
        cache = {}
        rows = db.execute('SELECT cell,isoform,sum(count) FROM observations WHERE gene=? GROUP BY cell,isoform', [gene]).fetchall()
        summary = collections.Counter()
        for cell, text, count in rows:
            count = int(count)
            if not count:
                continue
            if text not in cache:
                try:
                    cache[text] = classify_blocks(parse_blocks(text), ref)
                except ValueError:
                    cache[text] = None
            result = cache[text]
            if result is None:
                summary['invalid_structure_clusters'] += count
                continue
            E,I,boundary,retained,_,_ = result
            summary['classified_clusters'] += count
            summary['mixed_clusters' if E and I else 'exon_only_clusters' if E else 'intron_only_clusters' if I else 'neither_clusters'] += count
            for key, include in [('exon',E),('intron',I),('boundary',boundary),
                                 ('retained',retained),('intron_without_full_retention',I and not retained)]:
                if include:
                    row, col, val = triples[key]
                    row.append(cell_index[cell]); col.append(j); val.append(count)
        metrics.update(summary)
        audit.append({'gene':gene, 'gene_name':genes.gene_name.iloc[j],
                      'strict_intron_segments':len(ref['introns'].intervals),
                      'annotated_introns':len(ref['annotated_introns']), **summary})
        if (j+1)%2000 == 0:
            print('Classified genes', j+1, flush=True)
    db.close()
    matrices = {}
    for key,(row,col,value) in triples.items():
        matrix = sparse.coo_matrix((value,(row,col)), shape=G.shape,dtype=np.int64).tocsr()
        matrices[key] = matrix
        sparse.save_npz(output/(key+'_counts.npz'), matrix)
    genes.to_csv(output/'genes.tsv',sep='\t',index=False)
    cells['normalization_denominator'] = library
    cells.to_csv(output/'cells.tsv',sep='\t',index=False)
    pd.DataFrame(audit).fillna(0).to_csv(output/'gene_region_audit.tsv',sep='\t',index=False)
    checks = {'author_integer_G_reproduced':True, 'all_cells_preserved':len(cells)==G.shape[0],
              'nonnegative_integer_counts':all(m.dtype==np.int64 and np.all(m.data>=0) for m in matrices.values()),
              'region_counts_not_above_assigned_gene_clusters':all((G-m).data.min(initial=0)>=0 for m in matrices.values()),
              'intron_partition_conserved':(matrices['retained']+matrices['intron_without_full_retention']!=matrices['intron']).nnz==0}
    per_cell = cells[['cell','stage','batch']].copy()
    for key,m in matrices.items():
        per_cell[key+'_clusters'] = np.asarray(m.sum(1)).ravel()
        per_cell[key+'_genes'] = np.diff(m.indptr)
    joint = matrices['exon'].multiply(matrices['intron']).astype(bool)
    per_cell['joint_EI_genes'] = np.asarray(joint.sum(1)).ravel()
    per_cell.to_csv(output/'cell_region_support.tsv',sep='\t',index=False)
    region_summary = {key:{'cluster_incidence':int(m.sum()),'genes':int(np.count_nonzero(np.asarray(m.sum(0)))),
                           'cells':int(np.count_nonzero(np.asarray(m.sum(1))))} for key,m in matrices.items()}
    metrics['reference_unresolved_clusters'] = int(G.sum())-metrics['classified_clusters']-metrics['invalid_structure_clusters']
    report = {'passed':all(checks.values()), 'checks':checks, 'cells':len(cells),'genes':len(genes),
              'full_G_clusters':int(G.sum()),'classified_source_categories':dict(metrics),'regions':region_summary,
              'joint_EI_genes_per_cell_quantiles':np.quantile(per_cell.joint_EI_genes,[0,.25,.5,.75,1]).tolist(),
              'seconds':time.time()-started,'scope':'author_corrected_structure_region_incidence',
              'config_changed':False,'model_fit_started':False}
    metadata = {'assembly':'hg38','count_unit':'umi_region_incidence',
                'source_count_unit':'author_integer_UMI_clusters_pre_scatter_filter','source_count_field':'cluster',
                'exon_definition':'assigned gene annotated exon union; >=10 aligned block bp per corrected structure',
                'intron_definition':'union of within-transcript introns minus all reference exons on either strand; >=20 aligned block bp',
                'ambiguity_policy':'any annotated exon masks intronic sequence; invalid structures excluded; N gaps supply no coverage',
                'normalization_denominator_definition':'complete_longread_assigned_gene_integer_UMI_cluster_library',
                'annotation_source':{'reference':'GENCODE_v44_CellRanger_ARC_GRCh38_2024_A','sha256':digest(gtf)},
                'classification_validated':all(checks.values()),'classification_scope':'corrected_author_block_geometry',
                'molecule_states':'E/I overlapping region incidence; precursor/mature identity unresolved',
                'upstream_modifications':'LongcellPre splice/end correction and UMI structure correction precede released table',
                'sources':{p.relative_to(PROJECT).as_posix() if p.is_relative_to(PROJECT) else p.name:digest(p) for p in sources}}
    dump(output/'quantification.json',metadata)
    dump(output/'region_audit.json',report)
    print(report,flush=True)
    if not report['passed']:
        raise ValueError('Region count validation failed')


if __name__=='__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--source',type=Path,default=PROJECT/'test/GSE295932')
    p.add_argument('--processed-dir','--converted',dest='converted',type=Path)
    p.add_argument('--condition',choices=['stimulated','unstimulated','both'],default='both')
    p.add_argument('--gtf',type=Path,default=PROJECT.parent/'reference_documents/refdata-cellranger-arc-GRCh38-2024-A/genes/genes.gtf.gz')
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    a.converted=a.converted or condition_output(a.condition)
    a.output=a.output or a.converted/'regions'
    quantify(a.source.resolve(),a.converted.resolve(),a.gtf.resolve(),a.output.resolve(),a.condition)
