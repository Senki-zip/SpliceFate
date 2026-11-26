"""Joint conditional DM count responses with frozen spline and interaction bases."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
from scipy.interpolate import BSpline
from scipy.optimize import minimize
from scipy.special import softmax,gammaln,digamma
from scipy.spatial.distance import cdist,pdist
from sklearn.neighbors import NearestNeighbors
from .composition import settings as composition_settings
from .binning import CellSupport, settings as binning_settings


HYBRID = 'single_sum_hybrid'
HYBRID_SELECTION = 'all_passing_singles_hybrid'
GROUP_DEFAULTS = dict(min_group_rows=12, min_group_cell_ess=30., min_group_effective_rows=12.,
                      max_geometry_condition=10000., min_group_oob_gain=.005,
                      min_group_positive_gain_fraction=.6)


def hybrid_selection(c):
    return c['model'].get('selection_mode', HYBRID_SELECTION) == HYBRID_SELECTION


def observed_rows(counts, frame):
    return counts.sum(1) > (1e-8 if 'bin_id' in frame else 0.)


def group_settings(c):
    return {**GROUP_DEFAULTS, **c['model'].get('group_gates', {})}



def _spline_basis(x,knots,degree=3):
    x=np.asarray(x,float);knots=np.asarray(knots,float)
    lo,hi=knots[degree],knots[-degree-1]
    clamped=np.clip(x,lo,hi)
    basis=BSpline.design_matrix(clamped,knots,degree).toarray()
    # Linear-tangent saturating tails; freeze support range and knots at WT.
    coef=np.eye(len(knots)-degree-1);spline=BSpline(knots,coef,degree)
    left=x<lo;right=x>hi
    if left.any():basis[left]=spline(lo)+np.tanh(x[left]-lo)[:,None]*spline.derivative()(lo)
    if right.any():basis[right]=spline(hi)+np.tanh(x[right]-hi)[:,None]*spline.derivative()(hi)
    return basis[:,1:]


def _factor_basis(x,spec):
    if spec['kind']=='constant':return np.zeros((len(x),1))
    if spec['kind']=='linear':return ((x-spec['mean'])/spec['scale'])[:,None]
    return _spline_basis((x-spec['mean'])/spec['scale'],spec['knots'])-np.asarray(spec['basis_mean'])


def standardized_time(frame,t,spec):
    if 'context_time' not in spec:
        return (np.asarray(t)-spec['time_mean'])/spec['time_scale']
    z=np.zeros(len(frame))
    for branch,(mean,scale) in spec['context_time'].items():
        ix=frame.branch.astype(str).to_numpy()==branch
        z[ix]=(np.asarray(t)[ix]-mean)/scale
    return z


def state_basis(frame,t,spec=None,observed=None):
    good=np.ones(len(frame),bool) if observed is None else np.asarray(observed,bool)
    fitting=spec is None
    if fitting:
        spec={'stages':sorted(frame.stage.unique()),'batches':sorted(frame.batch.unique()),
              'branches':sorted(frame.branch.unique()),'context_time':{}}
        for branch in spec['branches']:
            ix=good & (frame.branch.to_numpy()==branch)
            values=np.asarray(t)[ix]
            spec['context_time'][str(branch)]=[float(values.mean()) if len(values) else 0.,
                                              max(float(values.std()),1e-6) if len(values) else 1.]
    columns=[np.ones(len(frame))];names=['intercept']
    for key,cats in [('stage',spec['stages']),('batch',spec['batches']),('branch',spec['branches'])]:
        for v in cats[1:]:columns.append((frame[key].values==v).astype(float));names.append(key+':'+str(v))
    columns.append(standardized_time(frame,t,spec));names.append('context_time')
    matrix=np.column_stack(columns)
    if fitting:
        keep=[];rank=0
        for j in range(matrix.shape[1]):
            new_rank=np.linalg.matrix_rank(matrix[good][:,keep+[j]],tol=1e-8)
            if new_rank>rank:keep.append(j);rank=new_rank
        spec['kept_columns']=keep;spec['column_names']=[names[j] for j in keep]
    return matrix[:,spec.get('kept_columns',list(range(matrix.shape[1])))],spec


def _stable_projection(base, target, ridge=1e-2):
    """Project onto a basis with ridge stabilization outside the intercept."""
    penalty = np.full(base.shape[1], ridge, dtype=float)
    penalty[0] = 0.
    gram = base.T @ base + np.diag(penalty)
    return np.linalg.lstsq(gram, base.T @ target, rcond=1e-8)[0]


def fit_design(X,frame,t,factor_indices,pair_indices,c,composition_features=(),observed_target=None):
    good=np.ones(len(frame),bool) if observed_target is None else np.asarray(observed_target,bool)
    B,state=state_basis(frame,t,observed=good);parts=[B];penalties=[np.r_[0,np.repeat(c['model']['state_ridge'],B.shape[1]-1)]]
    spec={'state':state,'factors':[],'interactions':[],'state_width':B.shape[1],'state_dependent':c['model'].get('state_dependent',False)}
    z=standardized_time(frame,t,state)
    offset=B.shape[1]
    for j in factor_indices:
        x=X[:,j];mean=float(np.mean(x[good]));scale=max(float(np.std(x[good])),1e-6);u=(x-mean)/scale;ug=u[good]
        f={'index':int(j),'mean':mean,'scale':scale,'range':[float(x[good].min()),float(x[good].max())],'p05':float(np.quantile(x[good],.05)),
           'kind':c['model']['response']}
        if np.ptp(x[good])<1e-6:f['kind']='constant';f['weak_identifiability']=True
        if f['kind']=='spline':
            internal=np.unique(np.quantile(ug,np.linspace(0,1,c['model'].get('spline_knots',0)+2)[1:-1]))
            internal=internal[(internal>ug.min()+1e-8)&(internal<ug.max()-1e-8)]
            if np.ptp(ug)<1e-6:f['kind']='linear'
            else:
                f['knots']=np.r_[np.repeat(ug.min(),4),internal,np.repeat(ug.max(),4)].tolist()
                temp=_spline_basis(u,f['knots']);f['basis_mean']=temp[good].mean(0).tolist()
        D=_factor_basis(x,f)
        if c['model'].get('state_dependent',False):D=np.column_stack([D,D*z[:,None]])
        # Remove intercept/state-only components so nuisance and RBP effects do not duplicate.
        proj=np.linalg.lstsq(B[good],D[good],rcond=1e-8)[0];D=D-B@proj
        f['state_projection']=proj.tolist();f['start']=offset;f['width']=D.shape[1]
        parts.append(D);penalties.append(np.full(D.shape[1],c['model']['ridge']));spec['factors'].append(f)
        offset+=D.shape[1]
    for a,b in pair_indices:
        # Residualize against state and the two local linear main effects.
        fa=next(f for f in spec['factors'] if f['index']==a);fb=next(f for f in spec['factors'] if f['index']==b)
        u=(X[:,a]-fa['mean'])/fa['scale'];v=(X[:,b]-fb['mean'])/fb['scale']
        D=(u*v)[:,None]
        base=np.column_stack([B,u,v])
        proj=_stable_projection(base[good],D[good]);D=D-base@proj
        if np.std(D[good])<1e-6:continue
        spec['interactions'].append({'a':int(a),'b':int(b),'projection':proj.tolist(),
                                     'projection_mode':'local_linear','column':offset})
        parts.append(D);penalties.append(np.array([c['model']['interaction_ridge']]));offset+=1
    spec['time_reference']={'context_time':state['context_time']}
    # Missing source observations are explicit nuisance gates, not filled fractions.
    spec['composition_masks']=[];spec['compositions']=[]
    if composition_features:
        for j in sorted({f['mask_index'] for f in composition_features}):
            parts.append(X[:,j,None]);penalties.append(np.array([c['model']['state_ridge']]))
            spec['composition_masks'].append(j)
            offset+=1
        cfg=composition_settings(c)
        minimum=(max(1,min(cfg['min_joint_bins'],int(frame.valid.sum()))) if 'bin_id' in frame
                 else cfg['min_joint_cells'])
        for feature in composition_features:
            observed=(X[:,feature['mask_index']]>0)&good
            if observed.sum()<minimum:continue
            value=X[:,feature['index']];scale=max(float(np.std(value[observed])),1e-8)
            D=value[:,None]/scale
            if c['model'].get('state_dependent',False):D=np.column_stack([D,D*z[:,None]])
            proj=_stable_projection(B[observed],D[observed])
            residual=D-B@proj
            if np.max(np.std(residual[observed],axis=0))<1e-6:continue
            residual[X[:,feature['mask_index']]<=0]=0
            support_key='joint_observed_bins' if 'bin_id' in frame else 'joint_observed_cells'
            saved={**feature,'scale':scale,'projection':proj.tolist(),'projection_mode':'state_only','start':offset,
                   'width':D.shape[1],support_key:int(observed.sum()),
                   'WT_range':[float(value[observed].min()),float(value[observed].max())]}
            spec['compositions'].append(saved)
            parts.append(residual);penalties.append(np.full(D.shape[1],composition_settings(c)['ridge']))
            offset+=D.shape[1]
    return np.column_stack(parts),np.concatenate(penalties),spec


def design(X,frame,t,spec):
    B,_=state_basis(frame,t,spec['state']);parts=[B]
    z=standardized_time(frame,t,spec['state'])
    flookup={}
    for f in spec['factors']:
        D=_factor_basis(X[:,f['index']],f)
        if spec['state_dependent']:D=np.column_stack([D,D*z[:,None]])
        parts.append(D-B@np.asarray(f['state_projection']));flookup[f['index']]=f
    interactions=spec['interactions']
    main=np.column_stack(parts) if any(term.get('projection_mode')!='local_linear' for term in interactions) else None
    for term in interactions:
        a,b=term['a'],term['b'];fa,fb=flookup[a],flookup[b]
        u=(X[:,a]-fa['mean'])/fa['scale'];v=(X[:,b]-fb['mean'])/fb['scale']
        base=np.column_stack([B,u,v]) if term.get('projection_mode')=='local_linear' else main
        parts.append((u*v)[:,None]-base@np.asarray(term['projection']))
    for j in spec.get('composition_masks',[]):parts.append(X[:,j,None])
    compositions=spec.get('compositions',[])
    main_composition=(np.column_stack(parts)
                      if any(feature.get('projection_mode')!='state_only' for feature in compositions) else None)
    for feature in compositions:
        D=X[:,feature['index'],None]/feature['scale']
        if spec['state_dependent']:D=np.column_stack([D,D*z[:,None]])
        base=B if feature.get('projection_mode')=='state_only' else main_composition
        residual=D-base@np.asarray(feature['projection'])
        residual[X[:,feature['mask_index']]<=0]=0
        parts.append(residual)
    return np.column_stack(parts)


def dm_objective(flat,D,C,penalty,kappa,offset=None):
    n,k=C.shape;coef=flat.reshape(D.shape[1],k-1)
    eta=D@coef
    if offset is not None:
        eta+=offset
    logits=np.column_stack([eta,np.zeros(n)])
    p=softmax(logits,axis=1);concentration=kappa*p;alpha=np.maximum(concentration,1e-12);N=C.sum(1)
    # Parameter-independent multinomial constants omitted; same rows in all comparisons.
    ll=gammaln(kappa)-gammaln(kappa+N)+np.sum(gammaln(alpha+C)-gammaln(alpha),axis=1)
    reg=.5*np.sum(penalty[:,None]*coef**2)
    loss=-ll.mean()+reg/n
    active=concentration>1e-12
    gp=np.where(active,-kappa*(digamma(alpha+C)-digamma(alpha))/n,0.)
    gl=p*(gp-np.sum(gp*p,axis=1,keepdims=True))
    grad=D.T@gl[:,:-1]+penalty[:,None]*coef/n
    return float(loss),grad.ravel()


def fit_count_design(D,C,pen,c,offset=None):
    """Fit one fixed WT design; used for full fits and cell bootstrap refits."""
    k=C.shape[1];theta=np.zeros((D.shape[1],k-1))
    if offset is None:
        mean=(C.sum(0)+.5)/(C.sum()+k*.5);theta[0]=np.log(mean[:-1]/mean[-1])
    fit=minimize(dm_objective,theta.ravel(),args=(D,C,pen,c['model']['dm_concentration'],offset),
                 method='L-BFGS-B',jac=True,options={'maxiter':c['model']['max_iterations'],'ftol':1e-10,'gtol':1e-6})
    if not fit.success:return {'status':'optimizer_not_converged','message':str(fit.message),'iterations':int(fit.nit)}
    if not np.isfinite(fit.fun):return {'status':'numerical_failure'}
    return {'status':'ok','coefficients':fit.x.reshape(D.shape[1],k-1).tolist(),
            'fit_objective':float(fit.fun),'iterations':int(fit.nit)}


def fit_event(C,X,frame,t,indices,pairs,c,composition_features=()):
    good=observed_rows(C,frame)
    if C.shape[1]>c['model'].get('max_event_branches',8):return {'status':'too_many_branches_unmodeled'}
    minimum=c['priors']['min_event_cells']
    if 'bin_id' in frame:
        minimum=max(1,min(binning_settings(c)['min_oob_bins'],int(frame.valid.sum())))
    if good.sum()<minimum:return {'status':'insufficient_event_counts'}
    D,pen,spec=fit_design(X,frame,t,indices,pairs,c,composition_features,good)
    if len(indices)==1 and spec['factors'][0]['kind']=='constant':
        return {'status':'constant_input_on_observed_cells','weak_identifiability':True}
    fit=fit_count_design(D[good],C[good],pen,c)
    if fit['status']!='ok':return fit
    coeff=np.asarray(fit['coefficients']);pred=softmax(np.column_stack([D@coeff,np.zeros(len(D))]),axis=1)
    B,_=state_basis(frame,t,spec['state']);diagnostics={}
    for j in indices:
        other=[l for l in indices if l!=j]
        controls=np.column_stack([B,X[:,other]]) if other else B
        resid=X[:,j]-controls@np.linalg.lstsq(controls[good],X[good,j],rcond=1e-8)[0]
        residual_sd=float(np.std(resid[good]));corr=[]
        for other_j in other:
            if np.std(X[good,j])>1e-8 and np.std(X[good,other_j])>1e-8:
                corr.append(float(abs(np.corrcoef(X[good,j],X[good,other_j])[0,1])))
        diagnostics[str(j)]={'conditional_residual_sd':residual_sd,'max_abs_input_correlation':max(corr,default=0),
            'weak_identifiability':residual_sd<c['model']['min_conditional_residual_sd'] or max(corr,default=0)>c['model']['collinearity_threshold']}
    support_indices=list(indices)+[f['index'] for f in spec.get('compositions',[])]
    support_values=np.column_stack([X[:,support_indices],standardized_time(frame,t,spec['state'])])
    ref=support_values[good];std=np.maximum(ref.std(0),1e-6);mean=ref.mean(0);u=(ref-mean)/std
    k=min(c['model'].get('support_neighbors',10),max(len(u)-1,0))
    # Match densities within independent trajectories, never borrow another branch.
    radii={};branch_neighbors={}
    ref_branches=frame.branch.to_numpy()[good]
    for branch in sorted(set(ref_branches)):
        ub=u[ref_branches==branch];kb=min(k,len(ub)-1)
        branch_neighbors[str(branch)]=kb
        radii[str(branch)]=(max(float(np.quantile(NearestNeighbors(n_neighbors=kb+1).fit(ub).kneighbors(ub)[0][:,kb],
            c['model']['support_distance_quantile'])),1e-6) if kb else 0.)
    radius=max(radii.values(),default=0.)
    return {'status':'ok','design':spec,'coefficients':coeff.tolist(),'kappa':c['model']['dm_concentration'],
        'input_indices':list(indices),'pairs':[list(p) for p in pairs],'fit_objective':fit['fit_objective'],'iterations':fit['iterations'],
        'counted_cells':int(good.sum()),'molecules':float(C.sum()),'training_only':True,'input_diagnostics':diagnostics,
        'support':{'mean':mean.tolist(),'scale':std.tolist(),'reference_inputs':ref,
                   'feature_indices':support_indices,'includes_time':True,'neighbors':k,'branch_radii':radii,'branch_neighbors':branch_neighbors,
                   'reference_branches':frame.branch[good].tolist(),'radius':max(radius,1e-6)},
        'WT_prediction':pred,
        'likelihood':'gamma_extended_DM_composite_pseudolikelihood_for_author_fractional_UMI_estimates' if c['model'].get('count_unit')=='author_fractional_umi_estimates' else 'Dirichlet_multinomial_composite_counts_fixed_shared_concentration',
        'composition_effects':'WT_learned_signed_coefficients_per_downstream_event; expression_weighted_local_isoform_contrasts',
        'state_input_indices':list(indices)+[f['index'] for f in spec.get('compositions',[])]}


def logits(model,X,frame,t):
    D=design(X,frame,t,model['design']);a=D@np.asarray(model['coefficients'])
    return np.column_stack([a,np.zeros(len(X))])


def predict_change(model,Xold,Xnew,frame,t,with_flags=True):
    if model.get('model_kind')==HYBRID:
        return _predict_hybrid(model,Xold,Xnew,frame,t,with_flags)
    old=np.asarray(model['WT_prediction'])
    indices=model.get('state_input_indices',model['input_indices'])
    changed_rows=np.any(Xold[:,indices]!=Xnew[:,indices],axis=1) if indices else np.zeros(len(Xold),bool)
    if not changed_rows.any():return old.copy(),np.zeros_like(old),[[] for _ in range(len(old))]
    eta_old=logits(model,Xold,frame,t);eta_new=logits(model,Xnew,frame,t)
    # Condition on the frozen WT model baseline and re-evaluate all joint terms once.
    new=softmax(np.log(np.maximum(old,1e-300))+(eta_new-eta_old),axis=1)
    new[~changed_rows]=old[~changed_rows]
    return new,new-old,prediction_flags(model,Xold,Xnew,frame,t,changed_rows) if with_flags else None


def prediction_flags(model,Xold,Xnew,frame,t,changed_rows,reference_rows=None):
    flags=[[] for _ in range(len(Xold))];changed=[j for j in model['input_indices'] if np.any(Xnew[:,j]!=Xold[:,j])]
    support=model['support'];support_indices=support.get('feature_indices',model['input_indices'])
    scale=np.asarray(support['scale']);mean=np.asarray(support['mean'])
    values=Xnew[:,support_indices]
    if support.get('includes_time'):
        time=standardized_time(frame,t,model['design']['state'])
        values=np.column_stack([values,time])
    if reference_rows is None:
        ref=np.asarray(support['reference_inputs'])
        reference_branches=np.asarray(support['reference_branches'])
    else:
        ref=Xold[reference_rows][:,support_indices]
        if support.get('includes_time'):ref=np.column_stack([ref,time[reference_rows]])
        reference_branches=frame.branch.to_numpy()[reference_rows]
    refu=(ref-mean)/scale;query=(values-mean)/scale
    for branch,rows in frame.groupby('branch').indices.items():
        ix=np.asarray(rows)[changed_rows[np.asarray(rows)]]
        if not len(ix):continue
        same=reference_branches==branch
        if not same.any():
            for i in ix:flags[i].append('branch_unsupported')
        else:
            kb=support.get('branch_neighbors',{}).get(str(branch),support.get('neighbors',1))
            kb=max(1,min(kb,int(same.sum())))
            distances=NearestNeighbors(n_neighbors=kb).fit(refu[same]).kneighbors(query[ix],return_distance=True)[0][:,kb-1]
            radius=support.get('branch_radii',{}).get(str(branch),support['radius'])
            for i in ix[distances>radius]:flags[i].append('conditional_extrapolation')
    for f in model['design']['factors']:
        j=f['index']
        if j not in changed:continue
        low,high=f['range']
        for i in np.flatnonzero(changed_rows&((Xnew[:,j]<low-1e-10)|(Xnew[:,j]>high+1e-10))):flags[i].append('input_outside_WT_range')
        for i in np.flatnonzero(changed_rows&(Xnew[:,j]<f.get('p05',low)-1e-10)):flags[i].append('input_below_observed_WT_p05')
        if model['input_diagnostics'][str(j)]['weak_identifiability']:
            for i in np.flatnonzero(changed_rows):flags[i].append('individual_effect_unresolved')
    for f in model['design'].get('compositions',[]):
        j=f['index'];mask=Xnew[:,f['mask_index']]>0;low,high=f['WT_range']
        out=mask&changed_rows&((Xnew[:,j]<low-1e-10)|(Xnew[:,j]>high+1e-10))
        for i in np.flatnonzero(out):flags[i].append('isoform_input_outside_WT_range')
    return [sorted(set(f)) for f in flags]


def _bounded_probability(p, wt_min, wt_max, cap, wt_row=None):
    span=wt_max-wt_min
    lower=np.maximum(1e-6,wt_min-cap*span)
    upper=np.minimum(1-1e-6,wt_max+cap*span)
    # A WT fit at numerical extremes can make the probability floor infeasible.
    lower=np.minimum(lower,wt_max)
    upper=np.maximum(upper,wt_min)
    if lower.sum()>1:
        lower=np.maximum(0.,wt_min-cap*span)
    if upper.sum()<1:
        upper=np.minimum(1.,wt_max+cap*span)
    if wt_row is not None:
        lower=np.minimum(wt_row,lower)
        upper=np.maximum(wt_row,upper)
    clipped=np.clip(p,lower,upper)
    rest=1.-clipped.sum(1,keepdims=True)
    room=np.where(rest>=0,upper-clipped,clipped-lower)
    available=room.sum(1,keepdims=True)
    fraction=np.divide(np.abs(rest),available,out=np.zeros_like(rest),where=available>0)
    return np.clip(clipped+np.sign(rest)*np.minimum(fraction,1.)*room,lower,upper)


def _group_term(group,X,frame,t):
    D,_=_group_design(X,frame,t,group['input_indices'],group)
    return np.column_stack([D[:,-1,None]*np.asarray(group['coefficient']),np.zeros(len(X))])


def _predict_hybrid(model,Xold,Xnew,frame,t,with_flags):
    old=np.asarray(model['WT_prediction'])
    indices=model['state_input_indices']
    if not np.any(Xold[:,indices]!=Xnew[:,indices]):
        return old.copy(),np.zeros_like(old),[[] for _ in old] if with_flags else None
    cap=model['output_cap_lambda']
    effect=np.zeros_like(old)
    flags=[set() for _ in old] if with_flags else None
    observed=np.asarray(model['wt_observed_rows'],dtype=int)
    for single in model['singles']:
        indices=single['state_input_indices']
        changed=np.any(Xold[:,indices]!=Xnew[:,indices],axis=1)
        if not changed.any():
            continue
        before=logits(single,Xold,frame,t)
        change=logits(single,Xnew,frame,t)[changed]-before[changed]
        responsive=np.any(change!=0,axis=1)
        moving=np.flatnonzero(changed)[responsive]
        wt=softmax(before[moving],axis=1)
        log_wt=np.log(np.maximum(wt,1e-300))
        raw=softmax(log_wt+change[responsive],axis=1)
        capped=_bounded_probability(raw,np.asarray(single['wt_prob_min']),np.asarray(single['wt_prob_max']),cap,wt_row=wt)
        effect[moving]+=single['weight']*(np.log(np.maximum(capped,1e-300))-log_wt)
        if flags is not None:
            for i,values in enumerate(prediction_flags(single,Xold,Xnew,frame,t,changed,observed)):
                flags[i].update(values)
    for group in model['groups']:
        indices=group['input_indices']
        changed=np.any(Xold[:,indices]!=Xnew[:,indices],axis=1)
        if not changed.any():
            continue
        before=_group_term(group,Xold,frame,t)[changed]
        after=_group_term(group,Xnew,frame,t)[changed]
        lower=np.asarray(group['wt_min']);upper=np.asarray(group['wt_max'])
        padding=cap*(upper-lower)
        bounded=np.clip(after,lower-padding,upper+padding)-np.clip(before,lower-padding,upper+padding)
        effect[changed]+=bounded
    new=old.copy()
    moving=np.any(effect!=0,axis=1)
    new[moving]=softmax(np.log(np.maximum(old[moving],1e-300))+effect[moving],axis=1)
    return new,new-old,[sorted(f) for f in flags] if flags is not None else None


def _group_design(X,frame,t,indices,spec=None):
    fitting=spec is None
    values=X[:,indices]
    if fitting:
        B,state=state_basis(frame,t)
        mean=values.mean(0)
        scale=np.maximum(values.std(0),1e-6)
        u=(values-mean)/scale
        local=np.column_stack([B,u])
        projection=_stable_projection(local,u[:,0]*u[:,1])
        spec=dict(input_indices=list(indices),state=state,mean=mean.tolist(),scale=scale.tolist(),
                  projection=projection.tolist())
    else:
        B,_=state_basis(frame,t,spec['state'])
        u=(values-spec['mean'])/spec['scale']
        local=np.column_stack([B,u])
        projection=np.asarray(spec['projection'])
    residual=u[:,0]*u[:,1]-local@projection
    return np.column_stack([local,residual]),spec


def _log_ratio(probability):
    log_probability=np.log(np.maximum(probability,1e-300))
    return log_probability[:,:-1]-log_probability[:,-1:]


def _fit_group_counts(C,X,frame,t,indices,c,baseline):
    D,group=_group_design(X,frame,t,indices)
    residual=D[:,-1:]
    if np.std(residual)<1e-6:
        return None,{'status':'constant_interaction'}
    penalty=np.array([c['model']['interaction_ridge']])
    fit=fit_count_design(residual,C,penalty,c,offset=_log_ratio(baseline))
    if fit['status']!='ok':
        return None,fit
    group['coefficient']=fit['coefficients'][0]
    terms=np.column_stack([residual*group['coefficient'],np.zeros(len(D))])
    group.update(wt_min=terms.min(0).tolist(),wt_max=terms.max(0).tolist())
    return group,fit


def refit_singles(singles,C,X,frame,times,c):
    fitted=[]
    for single in singles:
        m=fit_event(C,X,frame,times,single['input_indices'],[],c,single['design'].get('compositions',[]))
        if m['status']!='ok':
            return None
        fitted.append(dict(m,weight=single['weight']))
    return fitted


def _group_draws(X,C,frame,times,c,draws,observed,seed):
    from .network import bootstrap_samples, selection_settings
    if 'bin_id' in frame:
        minimum=binning_settings(c)['min_oob_bins']
        for X_tr,C_tr,valid_tr,X_oo,C_oo,valid_oo,oob_cells in draws:
            if oob_cells<minimum:
                yield None
                continue
            train=np.flatnonzero(valid_tr & observed_rows(C_tr,frame))
            oob=np.flatnonzero(valid_oo & observed_rows(C_oo,frame))
            yield _group_split(X_tr,C_tr,X_oo,C_oo,frame,times,train,oob,minimum)
    else:
        cfg=selection_settings(c)
        rows=np.flatnonzero(observed)
        for sample in bootstrap_samples(frame,observed,cfg['bootstrap_repeats'],seed):
            oob=np.ones(len(rows),bool);oob[sample]=False
            yield _group_split(X,C,X,C,frame,times,rows[sample],rows[oob],1)


def _group_split(X_tr,C_tr,X_oo,C_oo,frame,times,train,oob,minimum):
    if len(train)<4 or len(oob)<minimum:
        return None
    return (X_tr[train],C_tr[train],frame.iloc[train],times[train],
            X_oo[oob],C_oo[oob],frame.iloc[oob],times[oob])


def fit_supported_groups(X,C,frame,times,pairs,c,*,singles,bin_plan=None,C_cell_N=None,draws=(),seed=42):
    from .network import count_log_likelihood
    if not pairs:
        return [],[]
    cfg=group_settings(c)
    observed=observed_rows(C,frame)
    n=int(observed.sum())
    support=CellSupport(float(n),float(n))
    if bin_plan is not None:
        counted=C_cell_N>0
        weights=bin_plan.K_raw[observed]*counted[None,:]
        mass=weights.sum(1,keepdims=True)
        usable=mass[:,0]>0
        support=CellSupport.from_rows(weights[usable]/mass[usable])
        support=CellSupport(min(support.pooled_unique_cell_ess,float(counted.sum())),support.kernel_effective_rows)
    diagnostic=dict(observed_rows=n,**vars(support))
    if bin_plan is not None:
        diagnostic['expression_cell_support']=vars(CellSupport.from_rows(bin_plan.W_agg[observed]))
    if n<cfg['min_group_rows']:
        diagnostic['reason']='insufficient_group_rows'
    elif support.pooled_unique_cell_ess<cfg['min_group_cell_ess']:
        diagnostic['reason']='training_cell_ess'
    elif support.kernel_effective_rows<cfg['min_group_effective_rows']:
        diagnostic['reason']='kernel_effective_rows'
    if 'reason' in diagnostic:
        return [],[dict(diagnostic,input_indices=list(pair),accepted=False) for pair in pairs]

    groups=[];audit=[]
    X_wt,C_wt,fr_wt,t_wt=X[observed],C[observed],frame.iloc[observed],times[observed]
    baseline_wt=None
    splits=None
    baselines={}
    kappa=c['model']['dm_concentration']
    for pair in pairs:
        record=dict(diagnostic,input_indices=list(pair),accepted=False)
        audit.append(record)
        x=X_wt[:,pair]
        u=(x-x.mean(0))/np.maximum(x.std(0),1e-6)
        condition=float(np.linalg.cond(u)) if np.linalg.matrix_rank(u)==2 else np.inf
        record['geometry_condition']=condition if np.isfinite(condition) else None
        if condition>cfg['max_geometry_condition']:
            record['reason']='input_geometry'
            continue
        if baseline_wt is None:
            baseline_wt=single_sum_probability(singles,X_wt,fr_wt,t_wt)
        group,fit=_fit_group_counts(C_wt,X_wt,fr_wt,t_wt,pair,c,baseline_wt)
        if group is None:
            record['reason']=fit['status']
            continue
        gains=[]
        if splits is None:
            splits=list(_group_draws(X,C,frame,times,c,draws,observed,seed))
        for draw,split in enumerate(splits):
            if split is None:
                continue
            X_tr,C_tr,fr_tr,t_tr,X_oo,C_oo,fr_oo,t_oo=split
            if draw not in baselines:
                refits=refit_singles(singles,C_tr,X_tr,fr_tr,t_tr,c)
                baselines[draw]=None if refits is None else (
                    single_sum_probability(refits,X_tr,fr_tr,t_tr),
                    single_sum_probability(refits,X_oo,fr_oo,t_oo))
            if baselines[draw] is None:
                continue
            p_tr,p_oo=baselines[draw]
            fitted,_=_fit_group_counts(C_tr,X_tr,fr_tr,t_tr,pair,c,p_tr)
            if fitted is None:
                continue
            logits_oo=np.column_stack([_log_ratio(p_oo),np.zeros(len(p_oo))])
            p_interaction=softmax(logits_oo+_group_term(fitted,X_oo,fr_oo,t_oo),axis=1)
            gain=count_log_likelihood(C_oo,p_interaction,kappa)
            baseline=count_log_likelihood(C_oo,p_oo,kappa)
            gains.append(float(np.mean(gain-baseline)))
        mean=float(np.mean(gains)) if gains else None
        fraction=sum(g>0 for g in gains)/len(splits) if splits else 0.
        record.update(oob_fit_gain=mean,oob_positive_gain_fraction=fraction,oob_successes=len(gains))
        if mean is None or mean<cfg['min_group_oob_gain'] or fraction<cfg['min_group_positive_gain_fraction']:
            record['reason']='no_valid_oob_fits' if mean is None else 'insufficient_group_oob_gain'
            continue
        record['accepted']=True
        groups.append(group)
    return groups,audit


def single_sum_probability(singles,X,frame,times):
    p=sum(single['weight']*softmax(logits(single,X,frame,times),axis=1) for single in singles)
    return p/p.sum(1,keepdims=True)


def assemble_hybrid(singles,groups,group_diagnostics,X,frame,times,observed,c):
    wt=single_sum_probability(singles,X,frame,times)
    return dict(status='ok',model_kind=HYBRID,singles=singles,groups=groups,group_diagnostics=group_diagnostics,
                output_cap_lambda=c['model'].get('output_cap_lambda',.25),WT_prediction=wt,
                wt_observed_rows=np.flatnonzero(observed).tolist(),
                input_indices=[j for single in singles for j in single['input_indices']],
                state_input_indices=sorted({j for single in singles for j in single['state_input_indices']}),
                pairs=[group['input_indices'] for group in groups],likelihood=singles[0]['likelihood'],
                design=dict(compositions=[f for single in singles for f in single['design'].get('compositions',[])]))


def fit_models(data,frame,times,prior,X,c,event_limit=None,network_path=None,model_root=None,
               *,bin_plan=None,X_cell_factors=None):
    from .network import learn_network,network_summary,binned_event_draws,event_counts,selection_settings
    from .storage import ModelWriter
    writer=None
    models={};audit=[];findex={f:j for j,f in enumerate(prior['factor_ids'])}
    pairs={tuple(sorted(p['members'])) for p in prior['pairs']}
    hybrid=hybrid_selection(c)
    if not binning_settings(c)['enabled']:
        bin_plan=None
    network,views=learn_network(data,frame,times,prior,X,c,fit_event,logits,network_path,event_limit,
                                bin_plan=bin_plan,X_cell_factors=X_cell_factors,return_bin_views=True)
    if bin_plan is not None:
        frame,times=bin_plan.bin_frame,bin_plan.centers
    try:
        writer=ModelWriter(model_root) if model_root is not None else None
        summary=network_summary(network)
        event_ids=[r[0] for r in network.execute('SELECT DISTINCT event_index FROM single_relations ORDER BY event_index').fetchall()]
        for e in event_ids:
            event=data['events'][e]
            rows=network.execute('''SELECT candidate_json,response_score,selection_rank,WT_weight,selection_score,fit_gain,stability,oob_fit_gain,model_json
                FROM single_relations WHERE event_index=? AND selected ORDER BY selection_rank''',[e]).fetchall()
            candidates=[{**json.loads(raw),'WT_response_score':score,'selection_rank':rank,'WT_network_weight':weight,
                         'selection_score':quality,'fit_gain':gain,'oob_fit_gain':oob_gain,'stability':stability}
                        for raw,score,rank,weight,quality,gain,stability,oob_gain,_ in rows]
            indices=[findex[r['factor_id']] for r in candidates]
            scores={r['factor_id']:r['selection_score'] for r in candidates}
            accepted=set(scores)
            ranked_pairs=sorted((p for p in pairs if set(p)<=accepted),
                                key=lambda p:(-min(scores[f] for f in p),-sum(scores[f] for f in p),p))
            ipairs=[(findex[a],findex[b]) for a,b in ranked_pairs[:c['model']['max_interactions']]]
            C_cell=event_counts(data,e)
            C_cell_N=C_cell.sum(1)
            C=C_cell if bin_plan is None else bin_plan.K_eff@C_cell
            if C.shape[1]>c['model'].get('max_event_branches',8):
                model={'status':'too_many_branches_unmodeled'}
            elif not indices:
                model={'status':'no_supported_factor_response'}
            elif hybrid:
                singles=[dict(json.loads(row[-1]),weight=row[3]) for row in rows]
                draws=binned_event_draws(views,C_cell) if views is not None and ipairs else ()
                groups,diagnostics=fit_supported_groups(X,C,frame,times,ipairs,c,
                    singles=singles,bin_plan=bin_plan,C_cell_N=C_cell_N,draws=draws,
                    seed=selection_settings(c)['bootstrap_seed']+e)
                model=assemble_hybrid(singles,groups,diagnostics,X,frame,times,observed_rows(C,frame),c)
            else:
                composition=[f for f in prior.get('composition_features',[]) if f['factor_index'] in indices
                             and f['gene_id']!=event['gene']]
                model=(fit_event(C,X,frame,times,indices,ipairs,c,composition) if composition
                       else fit_event(C,X,frame,times,indices,ipairs,c))
            model.update(event_index=e,event_id=event['event'],gene_id=event['gene'],gene_name=event['gene_name'],
                         branch_count=C.shape[1],candidates=candidates,
                         network_candidate_count=len(prior['events'][str(e)]),
                         selection_policy=summary['selection_policy'])
            if writer is not None:writer.append(e,model)
            else:models[str(e)]=model
            audit.append({'event':event['event'],'gene':event['gene'],'status':model['status'],
                'network_candidates':len(prior['events'][str(e)]),'factors':len(indices),'interactions':len(model.get('pairs',[])),
                'isoform_inputs':len(model.get('design',{}).get('compositions',[])),
                'counted_cells':int((C_cell_N>0).sum()),
                'molecules':float(C_cell.sum()),'message':model.get('message','')})
        if network_path is not None:
            from .io import dump
            from pathlib import Path
            dump(Path(network_path).with_suffix('.summary.json'),summary)
            network.execute('CHECKPOINT')
    except BaseException:
        if writer is not None:writer.close()
        raise
    finally:
        network.close()
    if writer is not None:models=writer.finish()
    return models,pd.DataFrame(audit)
