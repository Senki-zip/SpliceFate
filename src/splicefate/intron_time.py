"""Independent intron/exon phase-based trajectory calibration.

Borrows the residual -> directional graph principle of velocity methods, not
their implementation. No scVelo/Regionvelocity runtime or latent-time API needed.
This is an auxiliary phase direction, not fitted transcription/splicing kinetics.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import spsolve



def load_intron_counts(directory,cells,assembly):
    root=Path(directory);meta=json.loads((root/'quantification.json').read_text())
    required=['intron_definition','exon_definition','count_unit','ambiguity_policy',
              'normalization_denominator_definition','annotation_source']
    for name in required:
        if not meta.get(name):raise ValueError('Missing intron quantification provenance: '+name)
    if meta.get('assembly')!=assembly:raise ValueError('Intron-count genome assembly mismatch')
    if meta.get('classification_validated') is not True:raise ValueError('Intron/exon region quantification needs validation metadata')
    if meta['count_unit'] not in ['umi_molecules','umi_region_incidence']:
        raise ValueError('Expected UMI molecules or explicit molecule-region incidence counts')
    celltable=pd.read_csv(root/'cells.tsv',sep='\t')
    genes=pd.read_csv(root/'genes.tsv',sep='\t',dtype=str)
    if celltable.cell.duplicated().any() or genes.gene.duplicated().any():raise ValueError('Duplicate intron-count identities')
    if set(celltable.cell)!=set(cells):raise ValueError('Intron counts must cover the same WT cell set')
    if 'normalization_denominator' not in celltable:raise ValueError('Complete count denominator must be provided per cell')
    order=pd.Index(celltable.cell).get_indexer(cells)
    denom=celltable.normalization_denominator.to_numpy(float)[order]
    if np.any(~np.isfinite(denom)) or np.any(denom<=0):raise ValueError('Invalid intron/exon normalization denominator')
    matrices=[]
    for name in ['exon_counts.npz','intron_counts.npz']:
        a=sparse.load_npz(root/name).tocsr()
        if a.shape!=(len(cells),len(genes)):raise ValueError('Intron/exon count axes differ')
        if np.any(~np.isfinite(a.data)) or np.any(a.data<0) or np.any(a.data!=np.floor(a.data)):
            raise ValueError('Raw nonnegative integer region counts required')
        matrices.append(a[order])
    return {'exons':matrices[0],'introns':matrices[1],'genes':genes.gene.tolist(),
            'denominator':denom,'provenance':meta}


def phase_directions(E,I,denominator,config):
    """Per-gene balanced residual, without cell smoothing or missing-value filling.

    The I ~= slope*E reference is estimated in high-count cells. A joint shift in
    total capture is removed by a shared externally supplied denominator. Each
    cell is scored only at genes where both regions have molecule support.
    """
    if E.shape!=I.shape:raise ValueError('Phase count matrices must have the same axes')
    denominator=np.asarray(denominator,float)
    if len(denominator)!=E.shape[0] or np.any(denominator<=0):raise ValueError('Invalid phase denominator')
    E=E.toarray().astype(float) if sparse.issparse(E) else np.asarray(E,float)
    I=I.toarray().astype(float) if sparse.issparse(I) else np.asarray(I,float)
    en=E/denominator[:,None];intr=I/denominator[:,None]
    state=np.zeros_like(en);direction=np.zeros_like(en);weight=np.zeros_like(en)
    fits=[]
    for j in range(E.shape[1]):
        observed=(E[:,j]>0)&(I[:,j]>0)
        report={'gene_index':j,'joint_observed_cells':int(observed.sum())}
        if observed.sum()<config['min_gene_cells'] or E[:,j].sum()<config['min_region_counts'] or I[:,j].sum()<config['min_region_counts']:
            fits.append({**report,'status':'insufficient_region_support'});continue
        x=en[:,j];y=intr[:,j]
        magnitude=x+y;cut=np.quantile(magnitude[observed],config['reference_quantile'])
        reference=observed&(magnitude>=cut)
        if reference.sum()<3:fits.append({**report,'status':'insufficient_reference_cells'});continue
        xx=x[reference];yy=y[reference]
        if np.ptp(xx)<=1e-12 or np.dot(xx,xx)<=1e-20:
            fits.append({**report,'status':'constant_exon_reference'});continue
        slope=float(xx@yy/(xx@xx))
        if not np.isfinite(slope) or slope<=1e-12:fits.append({**report,'status':'invalid_phase_slope'});continue
        residual=y-slope*x
        scale=max(float(np.std(x[observed])),1e-12)
        rscale=float(np.std(residual[observed]));reference_error=float(np.sqrt(np.mean((yy-slope*xx)**2)))
        if rscale<=1e-12:fits.append({**report,'status':'no_phase_deviation','slope':slope});continue
        # Scaling gives every gene a comparable exon-state and phase direction.
        state[:,j]=x/scale;direction[:,j]=residual/(slope*scale)
        count=np.minimum(E[:,j],I[:,j])
        weight[:,j]=np.where(observed,count/(count+config['count_shrinkage']),0)
        fits.append({**report,'status':'ok','slope':slope,'reference_cells':int(reference.sum()),
                     'reference_error':reference_error,'residual_sd':rscale})
    return state,direction,weight,fits


def phase_edge_scores(state,direction,weight,graph,config):
    """Compare each observed intron excess with changes to candidate G neighbors."""
    graph=sparse.csr_matrix(graph);coo=graph.tocoo();values=[];reliabilities=[]
    for i,j in zip(coo.row,coo.col):
        w=np.minimum(weight[i],weight[j]);valid=w>0
        number=int(valid.sum())
        if number<config['min_shared_genes']:
            values.append(0.);reliabilities.append(0.);continue
        delta=state[j]-state[i];v=direction[i]
        den=np.sqrt(np.sum(w*delta**2)*np.sum(w*v**2))
        score=float(np.sum(w*delta*v)/den) if den>1e-15 else 0.
        rel=min(1.,float(w.sum())/config['reliability_genes']) if den>1e-15 else 0.
        values.append(float(np.clip(score,-1,1)));reliabilities.append(rel)
    shape=graph.shape
    return (sparse.csr_matrix((values,(coo.row,coo.col)),shape=shape),
            sparse.csr_matrix((reliabilities,(coo.row,coo.col)),shape=shape))


def local_order_correction(base,graph,scores,reliability,weight,stage):
    """Correct only supported local edges, keeping every other DPT value exact.

    Minimize ||delta||² + weight*sum r_ij*(delta_j-delta_i-request_ij)².
    Neighbor gaps, a frozen root and stage bounds constrain the correction scale.
    No new global path distance or second trajectory is estimated.
    """
    base=np.asarray(base,float);stage=np.asarray(stage)
    result=base.copy();n=len(base)
    if weight==0:return result,{'supported_cells':0,'changed_cells':0}
    edges=sparse.csr_matrix(graph).tocoo()
    score=scores[edges.row,edges.col].A1;rel=reliability[edges.row,edges.col].A1
    keep=(rel>0)&np.isfinite(score)&(np.abs(score)>1e-12)&(edges.row!=edges.col)&(stage[edges.row]==stage[edges.col])
    src,dst=edges.row[keep],edges.col[keep];score=score[keep];rel=rel[keep]
    if not len(src):return result,{'supported_cells':0,'changed_cells':0}
    supported=np.unique(np.r_[src,dst]);root=int(np.argmin(base))
    free=supported[supported!=root]
    if not len(free):return result,{'supported_cells':len(supported),'changed_cells':0}
    gaps=np.abs(base[dst]-base[src]);positive=gaps[gaps>1e-12]
    if not len(positive):return result,{'supported_cells':len(supported),'changed_cells':0,'status':'no_local_DPT_resolution'}
    scale=float(np.median(positive));desired=score*np.maximum(gaps,scale)
    request=desired-(base[dst]-base[src])
    row=np.repeat(np.arange(len(src)),2);col=np.column_stack([src,dst]).ravel()
    B=sparse.csr_matrix((np.tile([-1.,1.],len(src)),(row,col)),shape=(len(src),n))[:,free]
    W=sparse.diags(weight*rel)
    delta=spsolve(sparse.eye(len(free),format='csr')+B.T@W@B,B.T@(weight*rel*request))
    delta=np.clip(delta,-weight*scale,weight*scale)
    for j,d in zip(free,delta):
        lo=base[stage==stage[j]].min();hi=base[stage==stage[j]].max()
        result[j]=np.clip(base[j]+d,lo,hi)
    return result,{'supported_cells':int(len(supported)),'changed_cells':int(np.count_nonzero(result!=base)),
                   'unchanged_unsupported_cells':int(n-len(supported)),'local_DPT_gap_scale':scale,
                   'max_absolute_correction':float(np.max(np.abs(result-base))),'root_fixed':True}


def calibrate_intron_time(base,frame,graph,counts,config):
    """Use phase evidence for a small supported correction of fixed DPT ordering."""
    result=base.copy();diagnostics=[]
    for branch,ix in frame.groupby('branch',sort=True).indices.items():
        ix=np.asarray(ix);subgraph=sparse.csr_matrix(graph[ix][:,ix])
        if len(ix)<config['min_gene_cells']:
            diagnostics.append({'branch':branch,'status':'insufficient_context_cells'});continue
        state,velocity,weight,fits=phase_directions(counts['exons'][ix],counts['introns'][ix],
                                                     counts['denominator'][ix],config)
        scores,rels=phase_edge_scores(state,velocity,weight,subgraph,config)
        active=np.asarray(rels.data)>0
        diag={'branch':branch,'selected_genes':len(counts['genes']),'supported_genes':sum(f['status']=='ok' for f in fits),
              'gene_fits':[{**f,'gene_id':counts['genes'][f['gene_index']]} for f in fits],
              'supported_directed_edges':int(active.sum()),'status':'ok' if active.any() else 'no_directional_support'}
        diagnostics.append(diag)
        if not active.any() or config['weight']==0:continue
        stage=pd.to_numeric(frame.stage_order.iloc[ix]).to_numpy()
        result[ix],correction=local_order_correction(base[ix],subgraph,scores,rels,config['weight'],stage)
        diag.update(correction)
    return result,{'method':'intron_exon_phase_direction','time_axis':'shared_WT_pseudotime_for_all_events',
        'contexts':diagnostics,'count_unit':counts['provenance']['count_unit'],
        'kinetic_rate_fit':False,'global_trajectory_refit':False,'time_unit':'relative_cell_state_order',
        'interpretation':'intron/exon phase-deviation auxiliary direction; not intron_total_as_age',
        'no_neighbor_count_smoothing':True,'no_S_completion':True}
