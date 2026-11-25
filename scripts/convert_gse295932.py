"""Convert Nontarget observations using integer author UMI-cluster counts.

Reuse SpliceFate's event discovery, common-anchor counting and chain identity.
Neither KO cells nor targeted WT/validation libraries enter these WT inputs.
"""
from pathlib import Path
import argparse
import collections
import gzip
import json
import re
import time

import duckdb
import numpy as np
import pandas as pd
from scipy import sparse

from splicefate.convert import (reference_genes,chain_record,discover_events,
    observation_path_groups,count_event_observations)
from splicefate.io import dump,digest

PROJECT=Path(__file__).resolve().parents[1]
DEFAULT_GTF=PROJECT.parent/'reference_documents/refdata-cellranger-arc-GRCh38-2024-A/genes/genes.gtf.gz'
DEFAULT_OUTPUT=PROJECT/'data/processed/GSE295932_Nontarget_hg38_cluster_20261121'


def parse_blocks(text):
    # LongcellPre uses Rsamtools pos (1-based) and inclusive end-1.
    blocks=sorted((int(a)-1,int(b)) for a,b in (x.split(',') for x in text.split('|')))
    if not blocks or any(a<0 or a>=b for a,b in blocks) or any(a[1]>=b[0] for a,b in zip(blocks,blocks[1:])):
        raise ValueError('Invalid author exon blocks')
    return blocks


def import_source(db,path,gsm):
    """Select pre-scatter-filter UMI clusters, never reads or weighted count."""
    db.execute("CREATE OR REPLACE TEMP TABLE raw AS SELECT * FROM read_csv(?,header=true,delim=' ')",[str(path)])
    columns={r[0] for r in db.execute('DESCRIBE raw').fetchall()}
    if not {'cell','gene','isoform','cluster','count','size'}<=columns:
        raise ValueError('Longcell input must provide cell/gene/isoform/cluster/count/size')
    invalid=db.execute('''SELECT count(*) FROM (SELECT try_cast(cluster AS DOUBLE) value FROM raw)
        WHERE value IS NULL OR NOT isfinite(value) OR value<0 OR value!=floor(value)''').fetchone()[0]
    if invalid:raise ValueError('cluster must contain nonnegative integer UMI-cluster counts')
    db.execute('INSERT INTO observations SELECT ?,?||\':\'||cell,split_part(gene,\'.\',1),isoform,CAST(cluster AS BIGINT) FROM raw',[gsm,gsm])
    return {'sample':gsm,'source_rows':db.execute('SELECT count(*) FROM raw').fetchone()[0],
        'selected_count_field':'cluster','selected_UMI_clusters':int(db.execute('SELECT sum(cluster) FROM raw').fetchone()[0]),
        'weighted_count_sum_for_comparison':float(db.execute('SELECT sum(count) FROM raw').fetchone()[0]),
        'fractional_weighted_count_rows_not_used':db.execute('SELECT count(*) FROM raw WHERE abs(count-round(count))>1e-9').fetchone()[0],
        'zero_weighted_count_rows_with_positive_cluster':db.execute('SELECT count(*) FROM raw WHERE count=0 AND cluster>0').fetchone()[0]}


SAMPLES = [('GSM8961552','GSM8961552_FT_sc_iso_count.Nontarget.txt.gz','unstimulated'),
           ('GSM8961554','GSM8961554_FT_stim_sc_iso_count.Nontarget.txt.gz','stimulated')]


def condition_samples(condition):
    if condition not in ('stimulated','unstimulated','both'):
        raise ValueError('Unknown condition: '+condition)
    return [row for row in SAMPLES if condition=='both' or row[2]==condition]


def condition_output(condition):
    return DEFAULT_OUTPUT if condition=='both' else PROJECT/f'data/processed/GSE295932_Nontarget_{condition}_hg38_cluster'


def convert(source,output,gtf,condition='both'):
    output.mkdir(parents=True,exist_ok=True)
    if (output/'conversion_status.json').exists() and json.loads((output/'conversion_status.json').read_text()).get('status')=='complete':
        raise FileExistsError('Complete conversion already exists; preserving it')
    state={'status':'running','started_at':time.strftime('%Y-%m-%dT%H:%M:%S%z'),'assembly':'hg38'}
    dump(output/'conversion_status.json',state)
    ref=reference_genes(gtf)
    db=duckdb.connect(str(output/'conversion_staging.duckdb'));db.execute('SET threads=4')
    db.execute("SET memory_limit='8GB'")
    audits=[];sources=[]
    db.execute('CREATE OR REPLACE TABLE observations(sample VARCHAR,cell VARCHAR,gene VARCHAR,isoform VARCHAR,count BIGINT)')
    frames=[]
    for gsm,name,sample_condition in condition_samples(condition):
        folder=source/'WT/Nontarget_controls'/gsm;path=folder/name
        axis=pd.read_csv(folder/'cells.tsv',sep='\t',dtype=str)
        frames.append(pd.DataFrame({'cell':gsm+':'+axis.barcode,'barcode':axis.barcode,
            'sample_id':gsm,'batch':gsm,'stage':sample_condition,'stage_order':0,
            'trajectory_direction_available':'false','perturbation':'Nontarget','species':'Homo_sapiens'}))
        audit=import_source(db,path,gsm)
        expected=set(axis.barcode)
        actual={x[0] for x in db.execute('SELECT DISTINCT cell FROM raw').fetchall()}
        if actual!=expected:raise ValueError('Control file cell axis differs from separated whitelist')
        audits.append({**audit,'condition':sample_condition,'cells':len(actual)})
        sources.extend([path,folder/'cells.tsv'])
    cells=pd.concat(frames,ignore_index=True);cell_index={s:i for i,s in enumerate(cells.cell)}
    found={x[0] for x in db.execute('SELECT DISTINCT gene FROM observations').fetchall()}
    unknown=sorted(found-set(ref));dump(output/'genes_not_in_reference.json',unknown)
    # Keep all full-transcriptome assigned genes on G; only reference-resolved
    # genes can contribute genomic splice events.
    gids=sorted(found);gene_index={s:i for i,s in enumerate(gids)}
    genes=pd.DataFrame([{'gene':g,**ref.get(g,{'gene_name':g,'chromosome':'','strand':'','gene_type':'unresolved_in_local_reference'})} for g in gids])
    grow=[];gcol=[];gvals=[]
    for cell,gene,count in db.execute('SELECT cell,gene,sum(count) FROM observations GROUP BY cell,gene').fetchall():
        grow.append(cell_index[cell]);gcol.append(gene_index[gene]);gvals.append(int(count))
    G=sparse.coo_matrix((gvals,(grow,gcol)),shape=(len(cells),len(gids)),dtype=np.int64).tocsr()
    library=np.asarray(G.sum(1)).ravel()
    if np.any(library<=0):raise ValueError('Empty control libraries')
    cells['complete_gene_UMI_library']=library
    db.execute('CREATE INDEX observations_gene ON observations(gene)')
    chains=[];events=[];ir=[];ic=[];iv=[];er=[];ec=[];ev=[];offsets=[0];qc=collections.Counter()
    for n,gene in enumerate(gids):
        if gene not in ref:continue
        reference=ref[gene];chrom=reference['chromosome'];strand=reference['strand']
        raw=db.execute('SELECT cell,isoform,sum(count) FROM observations WHERE gene=? GROUP BY cell,isoform',[gene]).fetchall()
        observations=[];chain_groups=collections.defaultdict(list);parsed={}
        for cell,text,count in raw:
            if count<=0:continue
            if text not in parsed:
                try:parsed[text]=parse_blocks(text)
                except ValueError:parsed[text]=None
            blocks=parsed[text]
            if blocks is None:qc['invalid_structure_weight']+=count;continue
            junctions=[(a[1],b[0]) for a,b in zip(blocks,blocks[1:])]
            if not junctions:qc['single_block_weight']+=count;continue
            ci=cell_index[cell];info=(chrom,strand,junctions,blocks)
            observations.append((ci,count,info));chain_groups[tuple(junctions)].append((ci,count,blocks));qc['valid_spliced_weight']+=count
        for path,members in sorted(chain_groups.items()):
            blocks=[list(x) for x in members[0][2]]
            for _,_,variant in members[1:]:
                blocks[0][0]=max(blocks[0][0],variant[0][0]);blocks[-1][1]=min(blocks[-1][1],variant[-1][1])
            if any(a>=b for a,b in blocks):qc['chain_common_coverage_unresolved']+=1;continue
            record=chain_record(gene,chrom,strand,list(path),blocks);record['gene_name']=reference['gene_name']
            col=len(chains);chains.append(record);sums=collections.Counter()
            for ci,count,_ in members:sums[ci]+=count
            for ci,count in sums.items():ir.append(ci);ic.append(col);iv.append(count)
        if observations:
            groups=observation_path_groups(observations)
            for event in discover_events(gene,chrom,strand,set(chain_groups),10):
                counts=count_event_observations(groups,event);totals=collections.Counter();observed=set()
                for (ci,branch),count in counts.items():
                    if count>0:totals[branch]+=count;observed.add(ci)
                if sum(totals.values())<10 or len(observed)<3 or len(totals)<2:continue
                event.update(gene_name=reference['gene_name'],branch_count=len(event['branches']),
                             total_UMIs=int(sum(totals.values())),informative_cells=len(observed),
                             count_unit='integer_umi',source_count_field='cluster')
                events.append(event);base=offsets[-1]
                for (ci,branch),count in counts.items():er.append(ci);ec.append(base+branch);ev.append(count)
                offsets.append(base+len(event['branches']))
        if (n+1)%1000==0:print(json.dumps({'genes':n+1,'total':len(gids),'chains':len(chains),'events':len(events)}),flush=True)
    I=sparse.coo_matrix((iv,(ir,ic)),shape=(len(cells),len(chains)),dtype=np.int64).tocsr()
    C=sparse.coo_matrix((ev,(er,ec)),shape=(len(cells),offsets[-1]),dtype=np.int64).tocsr()
    gp=output/'genes';ep=output/'events';gp.mkdir(exist_ok=True);ep.mkdir(exist_ok=True)
    for folder in [gp,ep]:cells.to_csv(folder/'cells.tsv',sep='\t',index=False)
    genes.to_csv(gp/'genes.tsv',sep='\t',index=False);sparse.save_npz(gp/'gene_counts.npz',G)
    sparse.save_npz(gp/'isoform_counts.npz',I);sparse.save_npz(ep/'branch_counts.npz',C)
    np.save(ep/'branch_offsets.npy',np.asarray(offsets,dtype=np.int64));dump(ep/'groups.json',events)
    pd.DataFrame([{**r,'structure':json.dumps(r['structure']),'exons':json.dumps(r['exons'])} for r in chains]).to_csv(gp/'isoforms.tsv',sep='\t',index=False)
    import anndata as ad
    import scanpy as sc
    a=ad.AnnData(G.astype(float),var=genes.set_index('gene'))
    sc.pp.highly_variable_genes(a,flavor='seurat_v3',n_top_genes=3000,subset=False)
    hv=a.var.highly_variable.to_numpy();genes.loc[hv].to_csv(ep/'genes.tsv',sep='\t',index=False)
    pd.DataFrame({'cell':cells.cell,'complete_gene_UMI_library':library,'scale_to_100000':1e5/library}).to_csv(ep/'normalization_factors.tsv',sep='\t',index=False)
    dump(gp/'quantification.json',{'assembly':'hg38','validated':True,
        'G_count_unit':'integer_gene_UMI_clusters','S_count_unit':'integer_umi',
        'source_count_field':'cluster','counting_stage':'author_corrected_UMI_clusters_before_scatter_small_cluster_filter',
        'source_weighted_count_field_used':False,'source_read_size_field_used':False,'S_rounding':False,
        'source':'Longcell author full-transcriptome Nontarget tables',
        'S_likelihood':'Dirichlet_multinomial_composite_counts_fixed_shared_concentration',
        'normalization_denominator':'complete_longread_assigned_gene_UMI_cluster_sum; no paired shortread G matrix provided in GEO attachments',
        'limitations':'Author-clustered/corrected structures, not a raw per-UMI whitelist; pre-scatter-filter counts may include erroneous small clusters'})
    sums=np.asarray(C.sum(0)).ravel();checks={
        'all_condition_control_cells_preserved':len(cells)==sum(a['cells'] for a in audits),'unique_cells':cells.cell.is_unique,
        'no_KO_or_targeted_WT_in_input':all(s['sample'] in ['GSM8961552','GSM8961554'] for s in audits),
        'G_integer':bool(np.issubdtype(G.dtype,np.integer)),
        'S_chain_integer':bool(np.issubdtype(I.dtype,np.integer)),
        'S_event_integer':bool(np.issubdtype(C.dtype,np.integer)),
        'G_count_conservation':int(G.sum())==sum(s['selected_UMI_clusters'] for s in audits),
        'chain_count_conservation':int(I.sum())==qc['valid_spliced_weight'],
        'finite_nonnegative_S':bool(np.isfinite(C.data).all() and np.all(C.data>=0)),
        'positive_G_libraries':bool(np.all(library>0)),'HVG3000':int(hv.sum())==3000,
        'events_nonempty':len(events)>0,'chain_ids_unique':len({r['isoform'] for r in chains})==len(chains),
        'event_ids_unique':len({e['event'] for e in events})==len(events),
        'event_weight_conservation':all(np.isclose(sums[offsets[i]:offsets[i+1]].sum(),e['total_UMIs']) for i,e in enumerate(events))}
    report={'status':'complete','assembly':'hg38','condition':condition,'cells':len(cells),'genes':len(genes),'chains':len(chains),'events':len(events),
        'branches':offsets[-1],'HVGs':int(hv.sum()),'samples':audits,'G_source':'author_longread_gene_UMI_totals',
        'S_count_unit':'integer_umi','source_count_field':'cluster',
        'counting_stage':'before_author_scatter_small_cluster_filter','QC':dict(qc),'reference_unresolved_genes':len(unknown),
        'coordinates':'source_1based_inclusive_to_0based_half_open; independently checked against GENCODE44 splice junctions',
        'full_length_claim':False,'all_WT_cells':True,'gene_holdout':False,'S_completion':False,
        'sources':{str(p):digest(p) for p in sources+[gtf]},'no_training_started':True}
    dump(output/'conversion_report.json',report);dump(output/'validation.json',{'passed':all(checks.values()),'checks':{k:bool(v) for k,v in checks.items()}})
    db.close();state.update(status='complete' if all(checks.values()) else 'validation_failed',finished_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
    dump(output/'conversion_status.json',state)
    with (output/'SHA256SUMS').open('w') as stream:
        for p in sorted(list(gp.iterdir())+list(ep.iterdir())):
            if p.is_file():stream.write(digest(p)+'  '+str(p.relative_to(output))+'\n')
    print(json.dumps({**report,'sources':'saved_in_report','checks':checks},indent=2),flush=True)
    if not all(checks.values()):raise ValueError('Conversion validation failed')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,default=PROJECT/'test/GSE295932')
    p.add_argument('--output',type=Path)
    p.add_argument('--condition',choices=['stimulated','unstimulated','both'],default='both')
    p.add_argument('--gtf',type=Path,default=DEFAULT_GTF);args=p.parse_args()
    args.output=args.output or condition_output(args.condition)
    try:convert(args.source,args.output,args.gtf,args.condition)
    except BaseException as exc:
        # Preserve a previously completed output when refusing to overwrite it.
        status=args.output/'conversion_status.json'
        if not status.exists() or json.loads(status.read_text()).get('status')!='complete':
            dump(status,{'status':'failed','error':str(exc)})
        raise
