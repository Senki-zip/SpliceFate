"""Expression reference clustering and count-aware WT splicing calibration."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import pairwise_distances
import anndata as ad
import scanpy as sc

from .io import normalized,aligned_cells


def reference_clustering(data,c):
    ids=pd.Index(data['genes'].gene).get_indexer(data['hvg'])
    X=normalized(data['G'][:,ids],data['scales'])
    obs=data['cells'].set_index('cell').copy()
    a=ad.AnnData(X=X,obs=obs,var=pd.DataFrame(index=data['hvg']))
    cfg=c['clustering'];key=cfg['cluster_key']
    enabled=cfg.get('enabled',True)
    # Explicit fixed HVG and normalization contract; no target-dependent HVG retention.
    ncomp=min(max(50,cfg['n_pcs']),len(obs)-1,X.shape[1]-1)
    if ncomp<2:raise ValueError('Insufficient PCA dimensions')
    sc.pp.pca(a,n_comps=ncomp,random_state=cfg['random_state'])
    sc.pp.neighbors(a,n_neighbors=min(cfg['n_neighbors'],len(obs)-1),
        n_pcs=min(cfg['n_pcs'],ncomp),random_state=cfg['random_state'])
    reused=enabled and key in obs and not obs[key].isna().any()
    if not enabled:
        # Ignore supplied Leiden labels as well as skipping the algorithm.
        a.obs[key]=pd.Categorical(['all']*len(obs))
    elif reused:a.obs[key]=a.obs[key].astype('category')
    else:sc.tl.leiden(a,resolution=cfg['resolution'],key_added=key,random_state=cfg['random_state'],flavor='igraph',n_iterations=2,directed=False)
    sc.tl.umap(a,random_state=cfg['random_state'])
    labels=a.obs[key].astype(str).to_numpy()
    if enabled and cfg['paga']['enabled']:
        sc.tl.paga(a,groups=key)
        cats=list(a.obs[key].cat.categories)
        conn=a.uns['paga']['connectivities'].toarray()
        _,labs=connected_components(sparse.csr_matrix(conn>=cfg['paga']['connectivity_threshold']),directed=False)
        for k in sorted(set(labs)):
            own=np.flatnonzero(labs==k)
            size=np.isin(labels,np.asarray(cats)[own]).sum()
            others=np.flatnonzero(labs!=k)
            if size<cfg['paga']['min_cells_per_lineage'] and len(others):
                scores=conn[np.ix_(own,others)]
                if scores.max()>0:
                    target=others[np.unravel_index(scores.argmax(),scores.shape)[1]]
                    labs[own]=labs[target]
        labels=np.array(['L'+str(labs[cats.index(v)]) for v in labels])
    counts=pd.Series(labels).value_counts()
    if not enabled:
        cohorts=np.repeat('global',len(labels));method='dpt';fallback='clustering_disabled_single_context'
    elif (counts<cfg['min_cells_per_cluster']).all():
        cohorts=np.repeat('global',len(labels));method='dpt';fallback='all_groups_small_global_DPT'
    else:
        cohorts=np.array([v if counts[v]>=cfg['min_cells_per_cluster'] else 'shared_small' for v in labels])
        method='dpt';fallback='small_groups_shared' if 'shared_small' in cohorts else 'none'
    # A shared model retains all cells; cohort fields describe trajectory contexts, not cell deletion.
    frame=data['cells'].copy()
    frame[key]=a.obs[key].astype(str).values;frame['model_context']=cohorts
    for j in range(min(2,a.obsm['X_umap'].shape[1])):frame[f'UMAP{j+1}']=a.obsm['X_umap'][:,j]
    meta={'cluster_sizes':a.obs[key].value_counts().to_dict(),'reference_labels_reused':reused,
          'trajectory_method':method,'context_fallback':fallback,'all_cells_retained':len(frame),
          'PAGA_enabled':bool(enabled and cfg['paga']['enabled']),
          'clustering_enabled':enabled,'G_only_clustering':enabled,'S_imputed':False}
    return frame,a.obsm['X_pca'][:,:cfg['n_pcs']],a.obsp['connectivities'].tocsr(),method,meta


def splicing_distances(data,c):
    n=len(data['cells']);numer=np.zeros((n,n));mass=np.zeros((n,n));shared=np.zeros((n,n),int)
    for e,event in enumerate(data['events']):
        C=data['C'][:,data['offsets'][e]:data['offsets'][e+1]].toarray().astype(float)
        N=C.sum(1);valid=N>0
        if valid.sum()<c['priors']['min_event_cells']:continue
        prior=(C.sum(0)+.5)/(C.sum()+.5*C.shape[1]);s=c['trajectory']['prior_strength']
        P=(C+s*prior)/(N[:,None]+s)
        # Hellinger distance, only over observed pairs; gene balancing for overlapping events.
        root=np.sqrt(P);distance=np.maximum(0,1-root@root.T)
        rel=N/(N+s);pairrel=np.minimum(rel[:,None],rel[None,:])
        pairrel*=valid[:,None]&valid[None,:]
        count=sum(x['gene']==event['gene'] for x in data['events'])
        numer+=pairrel*distance/count;mass+=pairrel/count
        shared+=valid[:,None]&valid[None,:]
    d=np.divide(numer,mass,out=np.zeros_like(numer),where=mass>0)
    reliability=np.minimum(1,shared/c['trajectory']['reliability_events'])
    reliability[shared<c['trajectory']['min_shared_events']]=0
    np.fill_diagonal(d,0);np.fill_diagonal(reliability,0)
    return d,reliability


def stage_rank(values,raw):
    """Keep known stage intervals ordered; refine order within a stage only."""
    out=np.zeros_like(raw,dtype=float);unique=np.sort(np.unique(values))
    for k,v in enumerate(unique):
        ix=np.flatnonzero(values==v);x=raw[ix]
        order=pd.Series(x).rank(method='average').to_numpy()/(len(ix)+1)
        out[ix]=(k+order)/max(1,len(unique))
    return out


def trajectory_contexts(frame,graph,c):
    isolated=frame.stage.isin(c['trajectory']['isolate_stages']).to_numpy()
    scope=frame.model_context.astype(str).to_numpy(copy=True)
    if 'shared_small' in scope:scope[:]='shared_all_WT'
    scope=np.array([s+(':isolated:'+str(v) if iso else ':development') for s,v,iso in zip(scope,frame.stage,isolated)])
    for context in sorted(set(scope)):
        subset=np.flatnonzero(scope==context)
        nc,comp=connected_components(graph[subset][:,subset],directed=False)
        for component in range(nc):
            yield f'{context}:component{component}',subset[comp==component]


def select_root(frame,ix,c):
    cfg=c['trajectory'];configured=cfg.get('root_cell')
    matches=np.flatnonzero(frame.cell.iloc[ix].values==configured) if configured else np.zeros(0,int)
    if len(matches):return int(matches[0]),'configured_cell',None
    if cfg.get('root_method','stage')=='cytotrace':
        if 'cytotrace_score' not in frame or 'cytotrace_qc' not in frame:
            raise ValueError('CytoTRACE root selection requires official scores and QC')
        scores=frame.cytotrace_score.iloc[ix].to_numpy(float)
        stages=pd.to_numeric(frame.stage_order.iloc[ix]).to_numpy(float)
        valid=frame.cytotrace_qc.iloc[ix].to_numpy(bool)&np.isfinite(scores)&(stages==stages.min())
        candidates=np.flatnonzero(valid)
        if not len(candidates):return None,'cytotrace_unavailable',None
        order=np.lexsort((frame.cell.iloc[ix[candidates]].to_numpy(),-scores[candidates]))
        root=int(candidates[order[0]])
        return root,'official_cytotrace',float(scores[root])
    stage=pd.to_numeric(frame.stage_order.iloc[ix]).to_numpy(float)
    return int(np.lexsort((frame.cell.iloc[ix].values,stage))[0]),'earliest_stage',None


def initial_times(frame,pcs,graph,method,c):
    if method!='dpt':raise ValueError('Initial ordering must use DPT')
    stage=pd.to_numeric(frame.stage_order).to_numpy(float)
    n=len(frame);times=np.zeros(n);branches=np.empty(n,dtype=object);diagnostics=[]
    for branch,ix in trajectory_contexts(frame,graph,c):
        branches[ix]=branch
        root,root_method,root_score=select_root(frame,ix,c)
        status='ok';used=method
        no_direction=('trajectory_direction_available' in frame and frame.trajectory_direction_available.iloc[ix].eq('false').all())
        if root is None or (no_direction and root_method not in ('configured_cell','official_cytotrace')):
            raw=np.zeros(len(ix));status='trajectory_direction_unavailable';used='state_only'
        elif len(ix)<c['trajectory']['min_component_cells']:
            raw=np.zeros(len(ix));status='insufficient_component_cells';used='stage_only'
        else:
            a=ad.AnnData(np.zeros((len(ix),1)));a.obsm['X_pca']=pcs[ix]
            # Use the exact connected component that defined the scoring context.
            # Rebuilding KNN here can create a graph disconnected from its root.
            subgraph=sparse.csr_matrix(graph[ix][:,ix])
            a.obsp['connectivities']=subgraph
            a.obsp['distances']=subgraph.copy()
            a.uns['neighbors']={'connectivities_key':'connectivities','distances_key':'distances',
                                'params':{'method':'umap','n_neighbors':c['clustering']['n_neighbors']}}
            a.uns['iroot']=root
            try:
                sc.tl.diffmap(a,n_comps=min(10,len(ix)-1));sc.tl.dpt(a,n_dcs=min(10,len(ix)-1))
                raw=a.obs.dpt_pseudotime.to_numpy()
                if not np.isfinite(raw).all():raise ValueError('Nonfinite DPT times')
            except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
                raw=np.zeros(len(ix));status='trajectory_unavailable:'+str(exc);used='stage_only'
        times[ix]=stage_rank(stage[ix],raw) if status!='trajectory_direction_unavailable' else .5
        diagnostics.append({'context':branch,'cells':len(ix),'method':used,'status':status,
                            'root_cell':frame.cell.iloc[ix[root]] if root is not None else None,
                            'root_method':root_method,'root_cytotrace_score':root_score})
    return times,branches.astype(str),diagnostics


def calibrate(base,branches,frame,pcs,graph,dS,reliability,weight):
    if weight==0:return base.copy()
    if not np.any(reliability*(graph.toarray()>0)):return base.copy()
    dg=pairwise_distances(pcs);scale=np.median(dg[graph.toarray()>0]);dg/=max(scale,1e-12)
    # Comparable reference distance scale; unsupported S pairs contribute exactly no change.
    ss=dS[(reliability>0)&(graph.toarray()>0)]
    ds=dS/max(float(np.median(ss)) if len(ss) else 1,1e-12)
    from .intron_time import local_order_correction
    times=base.copy();stage=pd.to_numeric(frame.stage_order).to_numpy()
    for branch in sorted(set(branches)):
        ix=np.flatnonzero(branches==branch)
        if len(ix)<3:continue
        mask=graph[ix][:,ix].toarray()>0
        gd=dg[np.ix_(ix,ix)];w=reliability[np.ix_(ix,ix)]
        if not np.any(w*mask):continue
        # Composition is undirected: preserve DPT orientation while asking similar
        # observed S states to have closer local gaps, not a second global ordering.
        distance=ds[np.ix_(ix,ix)]
        magnitude=np.clip(distance/np.maximum(gd+distance,1e-12),1e-6,1.)
        score=np.sign(base[ix][None,:]-base[ix][:,None])*magnitude
        subgraph=graph[ix][:,ix]
        times[ix],_=local_order_correction(base[ix],subgraph,sparse.csr_matrix(score*mask),
            sparse.csr_matrix(w*mask),weight,stage[ix])
    return times


def blend_velocity(times,velocity,valid,stage,branches,weight):
    out=times.copy()
    if not 0<=weight<=1:raise ValueError('velocity_weight must be within [0,1]')
    for branch in sorted(set(branches)):
        full=np.flatnonzero(branches==branch);levels=np.sort(np.unique(stage[full]))
        for k,level in enumerate(levels):
            ix=np.flatnonzero((branches==branch)&(stage==level)&valid)
            if not len(ix):continue
            ranks=pd.Series(velocity[ix]).rank(method='average').to_numpy()/(len(ix)+1)
            # Same stage intervals as G+S, even if kinetic QC fails in other stages.
            vt=(k+ranks)/max(1,len(levels));out[ix]=(1-weight)*out[ix]+weight*vt
    return out


def prepare_trajectory(data,c):
    frame,pcs,graph,method,meta=reference_clustering(data,c)
    if c['trajectory'].get('root_method','stage')=='cytotrace':
        from .cytotrace import score_contexts
        meta['cytotrace']=score_contexts(data,frame,graph,c)
    base,branches,diag=initial_times(frame,pcs,graph,method,c)
    frame['branch']=branches;frame['pseudotime_G']=base
    # All events share one calibrated WT time axis and all supported S features.
    intron_counts=None;intron_audit=[]
    calibration_mode=c['trajectory'].get('S_calibration','composition')
    if c['input'].get('intron_counts') and calibration_mode in ('intron_phase','hybrid'):
        from .intron_time import load_intron_counts
        intron_counts=load_intron_counts(c['input']['intron_counts'],list(frame.cell),c['assembly'])
    velocity=None
    if c['input'].get('velocity_times'):
        velocity=aligned_cells(pd.read_csv(c['input']['velocity_times'],sep='\t'),list(frame.cell))
        if 'qc_pass' not in velocity or 'latent_time' not in velocity:
            raise ValueError('Velocity input needs shared latent_time and qc_pass columns')
    unavailable={d['context'] for d in diag if d['status']=='trajectory_direction_unavailable'}
    t=base.copy()
    if calibration_mode in ('composition','hybrid'):
        ds,rel=splicing_distances(data,c)
        for branch in unavailable:
            ix=np.flatnonzero(branches==branch);rel[ix,:]=0;rel[:,ix]=0
        t=calibrate(base,branches,frame,pcs,graph,ds,rel,c['trajectory']['splice_weight'])
    if intron_counts is not None:
        from .intron_time import calibrate_intron_time
        valid=~frame.branch.isin(unavailable).to_numpy()
        if valid.any():
            ix=np.flatnonzero(valid)
            subset={**intron_counts,'exons':intron_counts['exons'][ix],
                    'introns':intron_counts['introns'][ix],'denominator':intron_counts['denominator'][ix]}
            t[ix],audit=calibrate_intron_time(t[ix],frame.iloc[ix].reset_index(drop=True),graph[ix][:,ix],
                                            subset,c['trajectory']['intron_phase'])
            intron_audit.append(audit)
    if velocity is not None:
        valid=velocity.qc_pass.eq(True).to_numpy()&np.isfinite(velocity.latent_time.to_numpy(float))
        valid &= ~frame.branch.isin(unavailable).to_numpy()
        t=blend_velocity(t,velocity.latent_time.to_numpy(float),valid,
                         pd.to_numeric(frame.stage_order).to_numpy(),branches,
                         c['trajectory'].get('velocity_weight',0.))
    frame['pseudotime_GS']=t
    frame['time_status']=[next(d['status'] for d in diag if d['context']==b) for b in branches]
    meta['contexts']=diag
    meta['time_axis']='shared_WT_pseudotime_for_all_events'
    meta['ordering_pipeline']='DPT_initial_then_local_supported_splicing_correction'
    meta['S_changes_mean_abs']=float(np.mean(np.abs(frame.pseudotime_GS-base)))
    meta['velocity_calibration']='disabled_missing_verified_input' if velocity is None else 'verified_external_shared_latent_time'
    meta['S_calibration_mode']=calibration_mode
    meta['intron_count_calibration']='enabled_independent_phase_direction' if intron_counts is not None else 'not_run_no_input_or_other_mode'
    meta['intron_calibration_audit']=intron_audit
    meta['stage_batch_confounds']={str(stage):list(g.batch.unique()) for stage,g in frame.groupby('stage') if g.batch.nunique()==1}
    return frame,t,pcs,meta
