"""Read-only evidence audit for isoform-conditioned cascade inputs.

These thresholds describe information support, not selected hyperparameters or
validated regulatory edges. No counts, model parameters or database rows change.
"""
import collections
import json
from pathlib import Path
import time

import duckdb
import numpy as np
import pandas as pd
from scipy import sparse

PROJECT=Path(__file__).resolve().parents[1]
DATA=PROJECT/'data/processed/NEMO_P21_author_20261117'
OUT=PROJECT/'runs/rbp_isoform_support_audit_20261120'


def log(**kw):
    print(json.dumps(kw),flush=True)


def classify_states(M,cells,context,threshold=3):
    """Independent composition dimensions among adequately observed cells."""
    N=np.asarray(M.sum(1)).ravel();valid=N>=threshold
    if not valid.any():return {'cells_ge_'+str(threshold):0,'variable_contexts':0,'variable_context_cells':0,'mixed_cells':0,'donors_supported':0}
    A=M[valid].toarray() if sparse.issparse(M) else M[valid]
    totals=A.sum(0);order=np.argsort(-totals);minor=order[1] if len(order)>1 else 0
    frac=A[:,minor]/N[valid];indices=np.flatnonzero(valid)
    mixed=int(((A>0).sum(1)>=2).sum());contexts=0;context_cells=0
    # At least two observed state groups in the same fixed Leiden x sample context.
    for _,ix in pd.Series(context[indices]).groupby(context[indices]).groups.items():
        ix=np.asarray(list(ix),int);p=frac[ix]
        if len(ix)>=30 and (p>=.2).sum()>=10 and (p<=.8).sum()>=10 and np.std(p)>=.1:
            contexts+=1;context_cells+=len(ix)
    donors=0
    for donor in sorted(cells.donor.unique()):
        ii=np.flatnonzero(cells.donor.to_numpy()[indices]==donor)
        if len(ii)>=20 and totals.size>1:
            sums=A[ii].sum(0)
            if (sums>=10).sum()>=2:donors+=1
    return {'cells_ge_'+str(threshold):int(valid.sum()),'variable_contexts':contexts,
            'variable_context_cells':context_cells,'mixed_cells':mixed,'donors_supported':donors}


def subsequence(short,long):
    if not short:return bool(long)
    return len(short)<len(long) and any(long[i:i+len(short)]==short for i in range(len(long)-len(short)+1))


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    d=duckdb.connect(str(PROJECT/'data/external/mouse_rbp_unified_20261116/mouse_rbp.duckdb'),read_only=True)
    mapping=dict(d.execute('SELECT original_factor_id,canonical_factor_id FROM family_sharing_factor_identity').fetchall())
    factors={f:(g,s) for f,g,s in d.execute('SELECT factor_id,gene_id,gene_symbol FROM factors').fetchall()}
    splice={mapping.get(f,f) for (f,) in d.execute('SELECT factor_id FROM candidate_splicing_factors').fetchall()}
    rbps={mapping.get(f,f) for (f,) in d.execute('SELECT factor_id FROM factor_catalog WHERE rbp2go_listed OR other_rbp_source').fetchall()};d.close()
    genes={factors[f][0]:factors[f][1] for f in rbps|splice if factors[f][0].startswith('ENSMUSG')}
    sg={factors[f][0] for f in splice};rg={factors[f][0] for f in rbps}
    prior=json.loads((PROJECT/'runs/NEMO_P21_restored_database_smoke_20261119/model/priors.json').read_text())
    focal=set(prior['gene_ids']);genes.update(dict(zip(prior['gene_ids'],prior['symbols'])))
    cells=pd.read_csv(DATA/'genes/cells.tsv',sep='\t')
    ref=pd.read_csv(PROJECT/'runs/NEMO_P21_CytoTRACE_DPT_20261120/cells.tsv',sep='\t',usecols=['cell','leiden','pseudotime_G'])
    cells=cells.merge(ref,on='cell',validate='one_to_one',sort=False)
    context=(cells.leiden.astype(str)+'|'+cells.sample_id).to_numpy()
    cells['audit_context']=context
    C=sparse.load_npz(DATA/'genes/isoform_counts.npz').tocsc();totals=np.asarray(C.sum(0)).ravel();nn=np.diff(C.indptr)
    global_info={'cells':C.shape[0],'observed_chains':C.shape[1],'molecules':int(totals.sum()),
        'single_molecule_chains':int((totals==1).sum()),'chains_ge20UMI_ge10cells':int(((totals>=20)&(nn>=10)).sum())}
    indices=collections.defaultdict(list);records={};all_gene_counts=collections.Counter();nj=collections.Counter()
    for chunk in pd.read_csv(DATA/'genes/isoforms.tsv',sep='\t',chunksize=50000):
        all_gene_counts.update(chunk.gene)
        chosen=chunk[chunk.gene.isin(genes)]
        for row in chosen.itertuples():
            j=row.Index;indices[row.gene].append(j)
            if totals[j]>=20 and nn[j]>=10:
                records[j]={'chain':row.isoform,'gene':row.gene,'structure':tuple(json.loads(row.structure)),
                            'exons':json.loads(row.exons),'full_length_validated':bool(row.full_length_validated)}
    log(phase='chain_table_loaded',genes=len(indices),supported_chains=len(records))
    rows=[];toprows=[];nested=[]
    for gid,columns in sorted(indices.items()):
        idx=np.asarray(columns);M=C[:,idx];N=np.asarray(M.sum(1)).ravel()
        supported=idx[(totals[idx]>=20)&(nn[idx]>=10)];strong=idx[(totals[idx]>=50)&(nn[idx]>=20)]
        ranked=idx[np.argsort(-totals[idx])];top2=ranked[:2]
        frac=totals[ranked]/max(1,totals[idx].sum())
        row={'gene':gid,'symbol':genes[gid],'splicing_candidate':gid in sg,'source_supported_RBP':gid in rg,'current_factor_universe':gid in focal,
             'chains':len(idx),'UMIs':int(N.sum()),'observed_cells':int((N>0).sum()),'cells_ge3_gene_UMI':int((N>=3).sum()),
             'cells_ge5_gene_UMI':int((N>=5).sum()),'supported_chains_20UMI_10cells':len(supported),
             'strong_chains_50UMI_20cells':len(strong),'top1_fraction':float(frac[0]),'top2_fraction':float(frac[:2].sum()),
             'supported_molecule_fraction':float(totals[supported].sum()/max(1,N.sum()))}
        if len(supported)>=2:
            X=C[:,supported];a=classify_states(X,cells,context,3);row.update({'chain_'+k:v for k,v in a.items()})
        else:row.update({'chain_cells_ge_3':0,'chain_variable_contexts':0,'chain_variable_context_cells':0,'chain_mixed_cells':0,'chain_donors_supported':0})
        # Audit nesting in top supported chains; suffix/prefix-only differences may be truncation.
        tops=supported[np.argsort(-totals[supported])][:32]
        for short in tops:
            shortpath=records[short]['structure'];parents=[int(long) for long in tops if short!=long and subsequence(shortpath,records[long]['structure'])]
            if parents:
                nested.append({'gene':gid,'symbol':genes[gid],'chain':records[short]['chain'],'UMIs':int(totals[short]),
                               'junctions':len(shortpath),'superchain_count':len(parents),
                               'possible_terminal_truncation':True})
        row['nested_top_supported_chains']=sum(r['gene']==gid for r in nested)
        rows.append(row)
        if gid in focal or gid in sg:
            for j in ranked[:10]:
                rec=records.get(int(j),{})
                toprows.append({'gene':gid,'symbol':genes[gid],'chain':rec.get('chain',str(j)),
                                'column_index':int(j),'UMIs':int(totals[j]),'cells':int(nn[j]),
                                'gene_UMI_fraction':float(totals[j]/max(1,N.sum())),
                                'junctions':len(rec['structure']) if rec else None})
    gene=pd.DataFrame(rows);gene.to_csv(OUT/'gene_chain_support.tsv',sep='\t',index=False)
    pd.DataFrame(toprows).to_csv(OUT/'top_chains.tsv',sep='\t',index=False)
    pd.DataFrame(nested).to_csv(OUT/'nested_chain_candidates.tsv',sep='\t',index=False)
    log(phase='gene_support_complete',genes=len(gene))
    # Anchor-eligible events provide independent evidence of internal composition differences.
    events=json.loads((DATA/'events/groups.json').read_text());offsets=np.load(DATA/'events/branch_offsets.npy')
    B=sparse.load_npz(DATA/'events/branch_counts.npz').tocsc();event_rows=[]
    for e,event in enumerate(events):
        gid=event['gene']
        if gid not in genes:continue
        M=B[:,offsets[e]:offsets[e+1]];bs=np.asarray(M.sum(0)).ravel();bn=np.diff(M.indptr)
        good=(bs>=20)&(bn>=10)
        row={'gene':gid,'symbol':genes[gid],'event':event['event'],'kind':event['kind'],
             'branches':len(bs),'supported_branches':int(good.sum()),'UMIs':int(bs.sum()),
             'minor_branch_UMIs':int(np.sort(bs)[-2]),'minor_branch_cells':int(bn[np.argsort(bs)[-2]])}
        row.update(classify_states(M,cells,context,3));event_rows.append(row)
    ef=pd.DataFrame(event_rows);ef.to_csv(OUT/'RBP_local_event_support.tsv',sep='\t',index=False)
    agg=ef.groupby('gene').agg(local_events=('event','size'),events_two_supported_branches=('supported_branches',lambda s:int((s>=2).sum())),
                               events_ge30cells_3UMI=('cells_ge_3',lambda s:int((s>=30).sum())),
                               events_variable_within_context=('variable_contexts',lambda s:int((s>0).sum())),
                               events_both_donors=('donors_supported',lambda s:int((s>=2).sum())))
    gene=gene.merge(agg,on='gene',how='left').fillna(0);gene.to_csv(OUT/'gene_chain_support.tsv',sep='\t',index=False)
    gene[gene.current_factor_universe].to_csv(OUT/'current_30_factor_support.tsv',sep='\t',index=False)
    def stats(f):
        ids=set(f.gene);ev=ef[ef.gene.isin(ids)]
        return {'genes_observed':len(f),'genes_two_supported_chains':int((f.supported_chains_20UMI_10cells>=2).sum()),
                'genes_two_strong_chains':int((f.strong_chains_50UMI_20cells>=2).sum()),
                'genes_chain_variable_within_context':int((f.chain_variable_contexts>0).sum()),
                'genes_with_local_event':int((f.local_events>0).sum()),
                'genes_with_two_branch_supported_event':int((f.events_two_supported_branches>0).sum()),
                'genes_with_event_30cells_3UMI':int((f.events_ge30cells_3UMI>0).sum()),
                'genes_event_variable_within_context':int((f.events_variable_within_context>0).sum()),
                'genes_event_both_donors':int((f.events_both_donors>0).sum()),
                'local_events':len(ev),'local_events_two_supported_branches':int((ev.supported_branches>=2).sum()),
                'median_gene_UMIs':float(f.UMIs.median()),'median_cells_ge3_UMI':float(f.cells_ge3_gene_UMI.median())}
    summary={'status':'complete','checked_at':time.strftime('%Y-%m-%dT%H:%M:%S%z'),'global':global_info,
             'thresholds':{'supported_chain_or_event_branch':'total molecules>=20 and cells>=10',
                           'strong_chain':'total molecules>=50 and cells>=20','composition_cell':'gene or event denominator>=3 molecules',
                           'variable_context':'fixed Leiden x sample; >=30 informative cells; >=10 cells with minor fraction>=0.2; >=10 with fraction<=0.8; fraction SD>=0.1',
                           'note':'Audit sensitivity thresholds for observed expression and composition support'},
             'splicing_factors':stats(gene[gene.splicing_candidate]),'all_RBPs':stats(gene[gene.source_supported_RBP]),
             'current_factor_universe':stats(gene[gene.current_factor_universe]),
             'current_seed_examples':gene[gene.symbol.isin(['Ptbp1','Rbfox2','Ddx5','Hnrnpm'])].to_dict('records'),
             'no_data_modified':True,'full_length_validated':False}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2));log(**summary)


if __name__=='__main__':main()
