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


def fit_design(X,frame,t,factor_indices,pair_indices,c,composition_features=(),observed_target=None):
    good=np.ones(len(frame),bool) if observed_target is None else np.asarray(observed_target,bool)
    B,state=state_basis(frame,t,observed=good);parts=[B];penalties=[np.r_[0,np.repeat(c['model']['state_ridge'],B.shape[1]-1)]]
    spec={'state':state,'factors':[],'interactions':[],'state_width':B.shape[1],'state_dependent':c['model']['state_dependent']}
    z=standardized_time(frame,t,state)
    rawbases={};feature_slices={};offset=B.shape[1]
    for j in factor_indices:
        x=X[:,j];mean=float(np.mean(x[good]));scale=max(float(np.std(x[good])),1e-6);u=(x-mean)/scale;ug=u[good]
        f={'index':int(j),'mean':mean,'scale':scale,'range':[float(x[good].min()),float(x[good].max())],'p05':float(np.quantile(x[good],.05)),
           'kind':c['model']['response']}
        if np.ptp(x[good])<1e-6:f['kind']='constant';f['weak_identifiability']=True
        if f['kind']=='spline':
            internal=np.unique(np.quantile(ug,np.linspace(0,1,c['model']['spline_knots']+2)[1:-1]))
            internal=internal[(internal>ug.min()+1e-8)&(internal<ug.max()-1e-8)]
            if np.ptp(ug)<1e-6:f['kind']='linear'
            else:
                f['knots']=np.r_[np.repeat(ug.min(),4),internal,np.repeat(ug.max(),4)].tolist()
                temp=_spline_basis(u,f['knots']);f['basis_mean']=temp[good].mean(0).tolist()
        D=_factor_basis(x,f);rawbases[j]=D
        if c['model']['state_dependent']:D=np.column_stack([D,D*z[:,None]])
        # Remove intercept/state-only components so nuisance and RBP effects do not duplicate.
        proj=np.linalg.lstsq(B[good],D[good],rcond=1e-8)[0];D=D-B@proj
        f['state_projection']=proj.tolist();f['start']=offset;f['width']=D.shape[1]
        parts.append(D);penalties.append(np.full(D.shape[1],c['model']['ridge']));spec['factors'].append(f)
        feature_slices[j]=(offset,offset+D.shape[1]);offset+=D.shape[1]
    main=np.column_stack(parts)
    for a,b in pair_indices:
        # Low-dimensional conditional interaction; residualize against all main effects.
        fa=next(f for f in spec['factors'] if f['index']==a);fb=next(f for f in spec['factors'] if f['index']==b)
        u=(X[:,a]-fa['mean'])/fa['scale'];v=(X[:,b]-fb['mean'])/fb['scale']
        D=(u*v)[:,None]
        proj=np.linalg.lstsq(main[good],D[good],rcond=1e-8)[0];D=D-main@proj
        if np.std(D[good])<1e-6:continue
        spec['interactions'].append({'a':int(a),'b':int(b),'projection':proj.tolist(),'main_width':main.shape[1],'column':offset})
        parts.append(D);penalties.append(np.array([c['model']['interaction_ridge']]));offset+=1
    spec['time_reference']={'context_time':state['context_time']}
    # Missing source observations are explicit nuisance gates, not filled fractions.
    spec['composition_masks']=[];spec['compositions']=[]
    if composition_features:
        for j in sorted({f['mask_index'] for f in composition_features}):
            parts.append(X[:,j,None]);penalties.append(np.array([c['model']['state_ridge']]))
            spec['composition_masks'].append(j)
            offset+=1
        main_composition=np.column_stack(parts)
        for feature in composition_features:
            observed=(X[:,feature['mask_index']]>0)&good
            if observed.sum()<composition_settings(c)['min_joint_cells']:continue
            value=X[:,feature['index']];scale=max(float(np.std(value[observed])),1e-8)
            D=value[:,None]/scale
            if c['model']['state_dependent']:D=np.column_stack([D,D*z[:,None]])
            proj=np.linalg.lstsq(main_composition[observed],D[observed],rcond=1e-8)[0]
            residual=D-main_composition@proj
            if np.max(np.std(residual[observed],axis=0))<1e-6:continue
            residual[X[:,feature['mask_index']]<=0]=0
            saved={**feature,'scale':scale,'projection':proj.tolist(),'start':offset,
                   'width':D.shape[1],'joint_observed_cells':int(observed.sum()),
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
    main=np.column_stack(parts)
    for term in spec['interactions']:
        a,b=term['a'],term['b'];fa,fb=flookup[a],flookup[b]
        u=(X[:,a]-fa['mean'])/fa['scale'];v=(X[:,b]-fb['mean'])/fb['scale']
        parts.append((u*v)[:,None]-main@np.asarray(term['projection']))
    for j in spec.get('composition_masks',[]):parts.append(X[:,j,None])
    main_composition=np.column_stack(parts)
    for feature in spec.get('compositions',[]):
        D=X[:,feature['index'],None]/feature['scale']
        if spec['state_dependent']:D=np.column_stack([D,D*z[:,None]])
        residual=D-main_composition@np.asarray(feature['projection'])
        residual[X[:,feature['mask_index']]<=0]=0
        parts.append(residual)
    return np.column_stack(parts)


def dm_objective(flat,D,C,penalty,kappa):
    n,k=C.shape;coef=flat.reshape(D.shape[1],k-1)
    logits=np.column_stack([D@coef,np.zeros(n)])
    p=softmax(logits,axis=1);alpha=np.maximum(kappa*p,1e-12);N=C.sum(1)
    # Parameter-independent multinomial constants omitted; same rows in all comparisons.
    ll=gammaln(kappa)-gammaln(kappa+N)+np.sum(gammaln(alpha+C)-gammaln(alpha),axis=1)
    reg=.5*np.sum(penalty[:,None]*coef**2)
    loss=-ll.mean()+reg/n
    gp=-kappa*(digamma(alpha+C)-digamma(alpha))/n
    gl=p*(gp-np.sum(gp*p,axis=1,keepdims=True))
    grad=D.T@gl[:,:-1]+penalty[:,None]*coef/n
    return float(loss),grad.ravel()


def fit_count_design(D,C,pen,c):
    """Fit one fixed WT design; used for full fits and cell bootstrap refits."""
    k=C.shape[1];theta=np.zeros((D.shape[1],k-1))
    mean=(C.sum(0)+.5)/(C.sum()+k*.5);theta[0]=np.log(mean[:-1]/mean[-1])
    fit=minimize(dm_objective,theta.ravel(),args=(D,C,pen,c['model']['dm_concentration']),
                 method='L-BFGS-B',jac=True,options={'maxiter':c['model']['max_iterations'],'ftol':1e-10,'gtol':1e-6})
    if not fit.success:return {'status':'optimizer_not_converged','message':str(fit.message),'iterations':int(fit.nit)}
    if not np.isfinite(fit.fun):return {'status':'numerical_failure'}
    return {'status':'ok','coefficients':fit.x.reshape(D.shape[1],k-1).tolist(),
            'fit_objective':float(fit.fun),'iterations':int(fit.nit)}


def fit_event(C,X,frame,t,indices,pairs,c,composition_features=()):
    good=C.sum(1)>0
    if C.shape[1]>c['model'].get('max_event_branches',8):return {'status':'too_many_branches_unmodeled'}
    if good.sum()<c['priors']['min_event_cells']:return {'status':'insufficient_event_counts'}
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
    radii={}
    ref_branches=frame.branch.to_numpy()[good]
    for branch in sorted(set(ref_branches)):
        ub=u[ref_branches==branch];kb=min(k,len(ub)-1)
        radii[str(branch)]=(max(float(np.quantile(NearestNeighbors(n_neighbors=kb+1).fit(ub).kneighbors(ub)[0][:,kb],
            c['model']['support_distance_quantile'])),1e-6) if kb else 0.)
    radius=max(radii.values(),default=0.)
    return {'status':'ok','design':spec,'coefficients':coeff.tolist(),'kappa':c['model']['dm_concentration'],
        'input_indices':list(indices),'pairs':[list(p) for p in pairs],'fit_objective':fit['fit_objective'],'iterations':fit['iterations'],
        'counted_cells':int(good.sum()),'molecules':float(C.sum()),'training_only':True,'input_diagnostics':diagnostics,
        'support':{'mean':mean.tolist(),'scale':std.tolist(),'reference_inputs':ref,
                   'feature_indices':support_indices,'includes_time':True,'neighbors':k,'branch_radii':radii,
                   'reference_branches':frame.branch[good].tolist(),'radius':max(radius,1e-6)},
        'WT_prediction':pred,
        'likelihood':'gamma_extended_DM_composite_pseudolikelihood_for_author_fractional_UMI_estimates' if c['model'].get('count_unit')=='author_fractional_umi_estimates' else 'Dirichlet_multinomial_composite_counts_fixed_shared_concentration',
        'composition_effects':'WT_learned_signed_coefficients_per_downstream_event; expression_weighted_local_isoform_contrasts',
        'state_input_indices':list(indices)+[f['index'] for f in spec.get('compositions',[])]}


def logits(model,X,frame,t):
    D=design(X,frame,t,model['design']);a=D@np.asarray(model['coefficients'])
    return np.column_stack([a,np.zeros(len(X))])


def predict_change(model,Xold,Xnew,frame,t,with_flags=True):
    old=np.asarray(model['WT_prediction'])
    indices=model.get('state_input_indices',model['input_indices'])
    changed_rows=np.any(Xold[:,indices]!=Xnew[:,indices],axis=1) if indices else np.zeros(len(Xold),bool)
    if not changed_rows.any():return old.copy(),np.zeros_like(old),[[] for _ in range(len(old))]
    eta_old=logits(model,Xold,frame,t);eta_new=logits(model,Xnew,frame,t)
    # Condition on the frozen WT model baseline and re-evaluate all joint terms once.
    new=softmax(np.log(np.maximum(old,1e-300))+(eta_new-eta_old),axis=1)
    new[~changed_rows]=old[~changed_rows]
    if not with_flags:return new,new-old,None
    flags=[[] for _ in range(len(Xold))];changed=[j for j in model['input_indices'] if np.any(Xnew[:,j]!=Xold[:,j])]
    support=model['support'];ref=np.asarray(support['reference_inputs']);support_indices=support.get('feature_indices',model['input_indices'])
    scale=np.asarray(support['scale']);mean=np.asarray(support['mean'])
    values=Xnew[:,support_indices]
    if support.get('includes_time'):values=np.column_stack([values,standardized_time(frame,t,model['design']['state'])])
    refu=(ref-mean)/scale;query=(values-mean)/scale
    for branch,rows in frame.groupby('branch').indices.items():
        ix=np.asarray(rows)[changed_rows[np.asarray(rows)]]
        if not len(ix):continue
        same=np.asarray(support['reference_branches'])==branch
        if not same.any():
            for i in ix:flags[i].append('branch_unsupported')
        else:
            distances=NearestNeighbors(n_neighbors=1).fit(refu[same]).kneighbors(query[ix],return_distance=True)[0][:,0]
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
    return new,new-old,[sorted(set(f)) for f in flags]


def fit_models(data,frame,times,prior,X,c,event_limit=None,network_path=None,model_root=None):
    from .network import learn_network,network_summary
    from .storage import ModelWriter
    writer=None
    models={};audit=[];findex={f:j for j,f in enumerate(prior['factor_ids'])}
    pairs={tuple(sorted(p['members'])) for p in prior['pairs']}
    network=learn_network(data,frame,times,prior,X,c,fit_event,logits,network_path,event_limit)
    try:
        writer=ModelWriter(model_root) if model_root is not None else None
        summary=network_summary(network)
        event_ids=[r[0] for r in network.execute('SELECT DISTINCT event_index FROM single_relations ORDER BY event_index').fetchall()]
        for e in event_ids:
            event=data['events'][e]
            rows=network.execute('''SELECT candidate_json,response_score,selection_rank,WT_weight,selection_score,fit_gain,stability,oob_fit_gain
                FROM single_relations WHERE event_index=? AND selected ORDER BY selection_rank''',[e]).fetchall()
            candidates=[{**json.loads(raw),'WT_response_score':score,'selection_rank':rank,'WT_network_weight':weight,
                         'selection_score':quality,'fit_gain':gain,'oob_fit_gain':oob_gain,'stability':stability}
                        for raw,score,rank,weight,quality,gain,stability,oob_gain in rows]
            indices=[findex[r['factor_id']] for r in candidates]
            scores={r['factor_id']:r['selection_score'] for r in candidates}
            accepted=set(scores)
            ranked_pairs=sorted((p for p in pairs if set(p)<=accepted),
                                key=lambda p:(-min(scores[f] for f in p),-sum(scores[f] for f in p),p))
            ipairs=[(findex[a],findex[b]) for a,b in ranked_pairs[:c['model']['max_interactions']]]
            C=data['C'][:,data['offsets'][e]:data['offsets'][e+1]].toarray()
            composition=[f for f in prior.get('composition_features',[]) if f['factor_index'] in indices
                         and f['gene_id']!=event['gene']]
            model=(fit_event(C,X,frame,times,indices,ipairs,c,composition) if composition
                   else fit_event(C,X,frame,times,indices,ipairs,c)) if indices else {'status':'no_supported_factor_response'}
            if C.shape[1]>c['model'].get('max_event_branches',8):model={'status':'too_many_branches_unmodeled'}
            model.update(event_index=e,event_id=event['event'],gene_id=event['gene'],gene_name=event['gene_name'],
                         branch_count=C.shape[1],candidates=candidates,
                         network_candidate_count=len(prior['events'][str(e)]),
                         selection_policy=summary['selection_policy'])
            if writer is not None:writer.append(e,model)
            else:models[str(e)]=model
            audit.append({'event':event['event'],'gene':event['gene'],'status':model['status'],
                'network_candidates':len(prior['events'][str(e)]),'factors':len(indices),'interactions':len(ipairs),
                'isoform_inputs':len(model.get('design',{}).get('compositions',[])),
                'counted_cells':int((C.sum(1)>0).sum()),
                'molecules':float(C.sum()),'message':model.get('message','')})
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
