"""Optional scVelo raw-layer adapter; never substitute mature AS/IR for precursors.

Calls upstream scVelo API, not copied kinetic code. Run in an environment with
scVelo installed and PYTHONPATH=src. No neighbor moment imputation is performed.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse

from .io import dump,digest


def validate_velocity_input(a):
    if not a.obs_names.is_unique or not a.var_names.is_unique:raise ValueError('Unique cell/gene axes required')
    qc=a.uns.get('velocity_input_qc',{})
    validated=qc.get('classification_validated')
    if not isinstance(validated,(bool,np.bool_)) or not bool(validated) or qc.get('unit')!='umi':
        raise ValueError('velocity_input_qc requires classification_validated=true and unit=umi; IR/branch counts are not precursor evidence')
    for key in ['assembly','annotation_source','quantification_method','precursor_definition','mature_definition']:
        if not qc.get(key):raise ValueError('Missing velocity quantification provenance: '+key)
    for layer in ['spliced','unspliced']:
        if layer not in a.layers:raise ValueError('Both genuine spliced/unspliced gene layers required')
        x=a.layers[layer];values=x.data if sparse.issparse(x) else np.asarray(x)
        if np.any(~np.isfinite(values)) or np.any(values<0) or np.any(values!=np.floor(values)):
            raise ValueError('Raw nonnegative integer UMI layers required')
    if not all(str(g).startswith('ENSMUSG') for g in a.var_names):raise ValueError('Mouse stable gene IDs required')


def main():
    p=argparse.ArgumentParser(description='Optional genuine precursor/mature RNA scVelo adapter')
    p.add_argument('--h5ad',required=True);p.add_argument('--prepared',required=True);p.add_argument('--output',required=True)
    p.add_argument('--min-shared-counts',type=int,default=20)
    p.add_argument('--min-genes',type=int,default=5);args=p.parse_args()
    import anndata as ad
    import scanpy as sc
    import scvelo as scv
    a=ad.read_h5ad(args.h5ad);validate_velocity_input(a)
    root=Path(args.prepared);frame=pd.read_csv(root/'cells.tsv',sep='\t');pcs=np.load(root/'prepared_arrays.npz')['pcs']
    if set(a.obs_names)!=set(frame.cell):raise ValueError('Velocity and prepared WT cell sets differ')
    a=a[frame.cell].copy();result=pd.DataFrame({'cell':frame.cell,'branch':frame.branch});audit=[]
    result['latent_time']=np.nan
    for branch,ix in frame.groupby('branch').indices.items():
        if len(ix)<5:audit.append({'branch':branch,'status':'too_few_cells'});continue
        sub=a[ix].copy();s=sub.layers['spliced'];u=sub.layers['unspliced']
        keep=(np.asarray(s.sum(0)).ravel()>=args.min_shared_counts)&(np.asarray(u.sum(0)).ravel()>=args.min_shared_counts)
        sub=sub[:,keep].copy()
        if sub.n_vars<args.min_genes:audit.append({'branch':branch,'status':'insufficient_precursor_gene_support'});continue
        # Dedicated kinetic normalization from complete layers, not inherited log PSI.
        denom=np.asarray((a.layers['spliced'][ix]+a.layers['unspliced'][ix]).sum(1)).ravel()
        if np.any(denom<=0):raise ValueError('Empty precursor/mature library')
        for layer in ['spliced','unspliced']:
            sub.layers[layer]=sparse.csr_matrix(sub.layers[layer]).multiply((1e4/denom)[:,None]).tocsr()
        sub.X=sub.layers['spliced'].copy();sub.obsm['X_pca']=pcs[ix]
        sc.pp.neighbors(sub,n_neighbors=min(20,len(ix)-1),use_rep='X_pca',random_state=42)
        earliest=np.argmin(pd.to_numeric(frame.stage_order.iloc[ix]).to_numpy());sub.uns['known_root']=int(earliest)
        # Explicit use_raw avoids pp.moments and neighbor-replaced molecular inputs.
        try:
            scv.tl.recover_dynamics(sub,var_names='all',use_raw=True,fit_connected_states=False,n_jobs=1)
            good=np.isfinite(sub.var['fit_likelihood'])&(sub.var['fit_likelihood']>=.1)
            if int(good.sum())<args.min_genes:raise ValueError('Insufficient kinetic model likelihood support')
            # Select by dynamical likelihood, not an additional steady-state R²
            # regression gate. The latter would reintroduce the model being replaced.
            scv.tl.velocity(sub,mode='dynamical',use_raw=True,min_r2=None,min_likelihood=.1)
            if int(sub.var['velocity_genes'].sum())<args.min_genes:raise ValueError('Insufficient selected velocity genes')
            scv.tl.velocity_graph(sub,xkey='spliced',approx=False,n_jobs=1)
            scv.tl.latent_time(sub,root_key='known_root',min_likelihood=.1)
            values=sub.obs['latent_time'].to_numpy(float)
            if not np.isfinite(values).all():raise ValueError('Nonfinite latent time')
            result.loc[ix,'latent_time']=values
            audit.append({'branch':branch,'status':'ok','kinetic_genes':int(good.sum())})
        except (ValueError,RuntimeError,KeyError,np.linalg.LinAlgError) as exc:
            audit.append({'branch':branch,'status':'kinetic_fit_unavailable','reason':str(exc)})
    result['qc_pass']=np.isfinite(result['latent_time'])
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True);result.to_csv(out,sep='\t',index=False)
    dump(str(out)+'.provenance.json',{'source':'scVelo_api','version':scv.__version__,'h5ad_sha256':digest(args.h5ad),
        'prepared_provenance_sha256':digest(root/'provenance.json'),'time_axis':'shared_WT_pseudotime_for_all_events',
        'no_neighbor_moment_imputation':True,'known_stage_root':True,'kinetic_unit_not_hours':True,'audit':audit,
        'velocity_gene_selection':'dynamical_likelihood>=0.1; no extra steady_state_R2_gate',
        'status':'usable' if bool(result.qc_pass.any()) else 'unavailable',
        'input_qc':dict(a.uns['velocity_input_qc'])})
    if not result.qc_pass.any():raise SystemExit('No cells passed kinetic calibration; inspect the provenance audit. No time values may be used.')


if __name__=='__main__':main()
