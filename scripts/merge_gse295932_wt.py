"""Merge matched WT cell identities; retain platform-specific count channels.

The paper states the 10x cDNA library was split for ONT/PacBio sequencing.
Counts are not pooled across platforms without a cross-platform molecule map.
Non-targeting controls from separate experiments retain cohort-prefixed IDs.
"""
from pathlib import Path
import collections
import csv
import gzip
import hashlib
import json
import time
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

ROOT=Path(__file__).resolve().parents[1]/'test/GSE295932'


def digest(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    source=ROOT/'files';out=ROOT/'WT/merged';out.mkdir(parents=True,exist_ok=True)
    evidence=ROOT/'metadata/paper_PMC12276307.xml'
    paragraphs=[''.join(p.itertext()) for p in ET.parse(evidence).getroot().iter('p')]
    quote=next(t for t in paragraphs if '10× full-length cDNA library was randomly split into two parts' in t)
    (out/'library_pairing_evidence.json').write_text(json.dumps({
        'paper':'https://doi.org/10.1038/s41467-025-60902-2','figure':'2B',
        'pmc':'https://pmc.ncbi.nlm.nih.gov/articles/PMC12276307/',
        'source_xml_sha256':digest(evidence),'source_excerpt':quote,
        'conclusion':'Same 10x cDNA library split between ONT and PacBio; shared barcodes represent matched cells'},indent=2))
    ontfile=source/'GSM8961556_Jurkat_UnconjGSP1_sc_iso_count.txt.gz'
    ont=pd.read_csv(ontfile,sep='\t',dtype={'cell':str,'gene':str,'isoform':str})
    with gzip.open(source/'GSM8961557_Jurkat_UnconjGSP1.hifi.isoseq_barcodes.tsv.gz','rt') as f:
        pb_original=[line.strip() for line in f if line.strip()]
    pb=[s.removesuffix('-1') for s in pb_original]
    if len(set(pb))!=len(pb):raise ValueError('PacBio barcode normalization collision')
    on=set(ont.cell);pn=set(pb);barcodes=sorted(on|pn);indices={b:i for i,b in enumerate(barcodes)}
    cells=pd.DataFrame({'cell':['GSE295932_WT_matched:'+b for b in barcodes],
        'barcode':barcodes,'cohort':'WT_matched_cDNA','condition':'WT',
        'ONT_observed':[b in on for b in barcodes],'PacBio_observed':[b in pn for b in barcodes],
        'ONT_source_sample':['GSM8961556' if b in on else '' for b in barcodes],
        'PacBio_source_sample':['GSM8961557' if b in pn else '' for b in barcodes]})
    cells.to_csv(out/'cells.tsv',sep='\t',index=False)
    links=[]
    for sample,platform,original,normalized in [
        ('GSM8961556','ONT',sorted(on),sorted(on)),
        ('GSM8961557','PacBio',pb_original,pb)]:
        for old,new in zip(original,normalized):
            links.append({'cell':cells.cell.iloc[indices[new]],'barcode':new,
                          'source_sample':sample,'platform':platform,'source_barcode':old})
    pd.DataFrame(links).to_csv(out/'source_cell_mapping.tsv',sep='\t',index=False)
    keys=list(zip(ont.gene,ont.isoform));feature_list=sorted(set(keys));findex={k:i for i,k in enumerate(feature_list)}
    counts=ont['count'].to_numpy(float)
    if np.any(~np.isfinite(counts)) or np.any(counts<0):raise ValueError('Invalid ONT source counts')
    matrix=sparse.coo_matrix((counts,([indices[b] for b in ont.cell],[findex[k] for k in keys])),
                             shape=(len(cells),len(feature_list))).tocsr()
    folder=out/'ONT';folder.mkdir(exist_ok=True)
    sparse.save_npz(folder/'counts.npz',matrix)
    pd.DataFrame({'feature':['ONT_'+hashlib.sha256((g+'|'+s).encode()).hexdigest()[:24] for g,s in feature_list],
                  'source_gene': [g for g,s in feature_list],'source_isoform':[s for g,s in feature_list],
                  'platform':'ONT','source_sample':'GSM8961556'}).to_csv(folder/'features.tsv',sep='\t',index=False)
    # The author ONT count field may contain fractional estimates; preserve it.
    ont_matrix_sum=float(matrix.sum());ont_source_sum=float(counts.sum())
    if not np.isclose(ont_matrix_sum,ont_source_sum,rtol=1e-12):raise ValueError('ONT count conservation failure')
    del ont,keys,counts
    with gzip.open(source/'GSM8961557_Jurkat_UnconjGSP1.hifi.isoseq_matrix.mtx.gz','rb') as f:
        pbmat=mmread(f).tocoo()
    features=pd.read_csv(source/'GSM8961557_Jurkat_UnconjGSP1.hifi.isoseq_genes.tsv.gz',sep='\t',header=None,names=['source_isoform','source_gene'],dtype=str)
    if pbmat.shape!=(len(features),len(pb)):raise ValueError('PacBio source dimensions inconsistent')
    values=pbmat.data
    if np.any(~np.isfinite(values)) or np.any(values<0):raise ValueError('Invalid PacBio counts')
    pbmatrix=sparse.coo_matrix((values,([indices[pb[j]] for j in pbmat.col],pbmat.row)),
                              shape=(len(cells),len(features))).tocsr()
    folder=out/'PacBio';folder.mkdir(exist_ok=True)
    sparse.save_npz(folder/'counts.npz',pbmatrix)
    features.insert(0,'feature',['PacBio_'+v for v in features.source_isoform]);features['platform']='PacBio';features['source_sample']='GSM8961557'
    features.to_csv(folder/'features.tsv',sep='\t',index=False)
    if not np.isclose(pbmatrix.sum(),pbmat.sum()):raise ValueError('PacBio count conservation failure')
    checks={'unique_merged_cell_IDs':cells.cell.is_unique,'union_cells':len(cells)==len(on|pn),
            'ONT_cells_preserved':set(cells.barcode[cells.ONT_observed])==on,
            'PacBio_cells_preserved':set(cells.barcode[cells.PacBio_observed])==pn,
            'ONT_unobserved_rows_zero':matrix[[indices[b] for b in pn-on]].nnz==0,
            'PacBio_unobserved_rows_zero':pbmatrix[[indices[b] for b in on-pn]].nnz==0,
            'source_mapping_unique':not pd.DataFrame(links).duplicated(['source_sample','source_barcode']).any(),
            'all_outputs_reload':sparse.load_npz(out/'ONT/counts.npz').shape==matrix.shape and sparse.load_npz(out/'PacBio/counts.npz').shape==pbmatrix.shape}
    if not all(checks.values()):raise ValueError(checks)
    # A common WT/control inventory also includes the two separate CRISPR cohorts.
    # Namespace them so incidental cross-library barcode matches are not merged.
    master=cells[['cell','barcode','cohort','condition']].copy();master['stimulated']='not_applicable'
    control_counts={}
    for gsm,cohort,stim in [('GSM8961552','CRISPR_unstimulated',False),('GSM8961554','CRISPR_stimulated',True)]:
        x=pd.read_csv(ROOT/f'WT/Nontarget_controls/{gsm}/cells.tsv',sep='\t',dtype=str)
        part=pd.DataFrame({'cell':cohort+':'+x.barcode,'barcode':x.barcode,'cohort':cohort,'condition':'Nontarget','stimulated':str(stim).lower()})
        master=pd.concat([master,part],ignore_index=True);control_counts[gsm]=len(part)
    if not master.cell.is_unique:raise ValueError('Cohort namespace collision')
    master.to_csv(ROOT/'WT/all_WT_control_cells.tsv',sep='\t',index=False)
    report={'status':'complete','finished_at':time.strftime('%Y-%m-%dT%H:%M:%S%z'),
        'passed':all(checks.values()),'checks':{k:bool(v) for k,v in checks.items()},
        'explicit_WT_unique_cells':len(cells),'ONT_cells':len(on),'PacBio_cells':len(pn),
        'matched_cells':len(on&pn),'ONT_only_cells':len(on-pn),'PacBio_only_cells':len(pn-on),
        'ONT_matrix_shape':list(matrix.shape),'PacBio_matrix_shape':list(pbmatrix.shape),
        'ONT_count_sum':ont_matrix_sum,'PacBio_count_sum':float(pbmatrix.sum()),
        'matrix_orientation':'cells x platform-specific isoform features; both use merged/cells.tsv',
        'platform_unobserved_rows':'zero storage with explicit availability flags; not measured biological zeros',
        'cross_platform_counts_pooled':False,'cross_platform_UMI_deduplication_possible_with_current_ONT_table':False,
        'Nontarget_controls_kept_as_separate_cohorts':control_counts,
        'all_WT_and_Nontarget_cell_inventory':len(master),'original_files_modified':False,
        'not_model_input_conversion':'Counts and feature IDs preserved; human/hg38 adaptation and count-unit review remain separate'}
    (out/'summary.json').write_text(json.dumps(report,indent=2))
    with (out/'SHA256SUMS').open('w') as f:
        for p in sorted(out.rglob('*')):
            if p.is_file() and p.name!='SHA256SUMS':f.write(digest(p)+'  '+str(p.relative_to(out))+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
