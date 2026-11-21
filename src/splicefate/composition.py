"""Observed local isoform-state inputs and synchronous composition propagation.

Common-anchor alternatives define isoform groups without inventing full-length
transcripts. Missing observations have zero gated contribution, never imputed
fractions. Each downstream event learns its own signed composition coefficients.
"""
from __future__ import annotations
import numpy as np
from scipy.special import softmax


DEFAULTS={'enabled':True,'max_events_per_factor':1,'max_branches':4,
          'min_branch_molecules':20,'min_branch_cells':10,'min_joint_cells':20,
          'count_shrinkage':2.,'ridge':10.}


def settings(c):
    return {**DEFAULTS,**c['model'].get('isoform_state',{})}


def build_states(data,prior,X,c):
    cfg=settings(c);n=len(X);states=[];features=[];columns=[X];offset=X.shape[1]
    if not cfg['enabled']:
        prior['isoform_states']=[];prior['composition_features']=[]
        return prior,X
    by_gene={g:i for i,g in enumerate(prior['gene_ids'])};options={}
    for e,event in enumerate(data['events']):
        if event['gene'] not in by_gene or len(event['branches'])>cfg['max_branches']:continue
        C=data['C'][:,data['offsets'][e]:data['offsets'][e+1]].tocsc()
        totals=np.asarray(C.sum(0)).ravel();cells=np.diff(C.indptr)
        if ((totals>=cfg['min_branch_molecules'])&(cells>=cfg['min_branch_cells'])).sum()<2:continue
        N=np.asarray(C.sum(1)).ravel();valid=N>0
        if valid.sum()<cfg['min_joint_cells']:continue
        # A near-constant abundant branch alone does not supply composition information.
        fractions=C[valid].toarray()/N[valid,None]
        score=float(np.sum(totals)-np.max(totals))
        options.setdefault(event['gene'],[]).append((-score,event['event'],e,float(np.max(np.std(fractions,axis=0)))))
    for gene in sorted(options):
        factor=by_gene[gene]
        for _,_,e,variation in sorted(options[gene])[:cfg['max_events_per_factor']]:
            if variation<=1e-8:continue
            event=data['events'][e];C=data['C'][:,data['offsets'][e]:data['offsets'][e+1]].toarray().astype(float)
            N=C.sum(1);observed=N>0;totals=C.sum(0);pooled=(totals+.5)/(totals.sum()+.5*C.shape[1])
            P=np.divide(C,N[:,None],out=np.full_like(C,np.nan),where=N[:,None]>0)
            reliability=N/(N+cfg['count_shrinkage'])
            # K-1 contrasts avoid the sum-to-one redundancy. No per-context cell quotas.
            reference=int(np.argmax(totals));branches=[k for k in range(C.shape[1]) if k!=reference]
            mask_index=offset;columns.append(observed[:,None].astype(float));offset+=1
            sid=len(states);state={'state_index':sid,'factor_index':factor,'factor_id':prior['factor_ids'][factor],
                'gene_id':gene,'source_event_index':e,'source_event':event['event'],'reference_branch':reference,
                'branch_count':C.shape[1],'pooled_composition':pooled.tolist(),'observed':observed.tolist(),
                'counts':C.tolist(),'reliability':reliability.tolist(),'mask_index':mask_index,
                'feature_indices':[],'interpretation':'common_anchor_local_isoform_groups'}
            for k in branches:
                values=np.zeros(n);values[observed]=X[observed,factor]*(P[observed,k]-pooled[k])*reliability[observed]
                columns.append(values[:,None]);state['feature_indices'].append(offset)
                features.append({'index':offset,'factor_index':factor,'factor_id':state['factor_id'],'gene_id':gene,
                    'state_index':sid,'source_event_index':e,'source_event':event['event'],
                    'branch_index':k,'reference_branch':reference,'mask_index':mask_index,
                    'kind':'expression_times_observed_composition_contrast'})
                offset+=1
            states.append(state)
    prior['isoform_states']=states;prior['composition_features']=features
    prior['isoform_state_policy']={'unit':'local_isoform_group_fraction_times_log_normalized_gene_expression',
        'missing':'zero gated contribution; no assigned fraction','weight_learning':'signed downstream-event-specific WT regression',
        'full_length_or_protein_annotation_required':False,'settings':cfg}
    return prior,np.column_stack(columns)


def apply_intervention(X,indices,strength,prior):
    if not np.isfinite(strength) or not 0<=strength<=1:raise ValueError('Dose must be within [0,1]')
    out=X.copy()
    if strength:
        out[:,indices]=np.log1p((1-strength)*np.expm1(X[:,indices]))
        for feature in prior.get('composition_features',[]):
            j=feature['factor_index']
            if j not in indices:continue
            ratio=np.divide(out[:,j],X[:,j],out=np.zeros(len(X)),where=X[:,j]!=0)
            out[:,feature['index']]=X[:,feature['index']]*ratio
    return out


def next_inputs(X,direct,event_wt,event_new,models,prior):
    """Every round is recalculated relative to observed WT, never accumulated."""
    result=direct.copy();changes=[];by_state={}
    for feature in prior.get('composition_features',[]):by_state.setdefault(feature['state_index'],[]).append(feature)
    for state in prior.get('isoform_states',[]):
        e=str(state['source_event_index']);model=models.get(e)
        if not model or model['status']!='ok' or e not in event_new:continue
        observed=np.asarray(state['observed'],bool);C=np.asarray(state['counts'],float);N=C.sum(1)
        # Exact no-op preserves zeros and all original features byte-for-byte.
        delta=np.log(np.maximum(event_new[e],1e-300))-np.log(np.maximum(event_wt[e],1e-300))
        changed=observed&np.any(delta!=0,axis=1)
        if not changed.any():continue
        P=C[changed]/N[changed,None]
        updated=softmax(np.log(np.maximum(P,1e-12))+delta[changed],axis=1)
        rel=np.asarray(state['reliability']);factor=state['factor_index'];pooled=np.asarray(state['pooled_composition'])
        for feat in by_state.get(state['state_index'],[]):
            k=feat['branch_index']
            result[changed,feat['index']]=direct[changed,factor]*(updated[:,k]-pooled[k])*rel[changed]
        changes.append({'factor_id':state['factor_id'],'source_event':state['source_event'],
                        'observed_cells_updated':int(changed.sum()),'max_composition_delta':float(np.max(np.abs(updated-P)))})
    return result,changes


def state_diagnostics(prior,models):
    used={f['state_index'] for m in models.values() if m['status']=='ok' for f in m['design'].get('compositions',[])}
    return [{'factor_id':s['factor_id'],'gene_id':s['gene_id'],'source_event':s['source_event'],
             'observed_cells':int(np.sum(s['observed'])),'source_event_fitted':models.get(str(s['source_event_index']),{}).get('status')=='ok',
             'used_by_downstream_model':s['state_index'] in used}
            for s in prior.get('isoform_states',[])]
