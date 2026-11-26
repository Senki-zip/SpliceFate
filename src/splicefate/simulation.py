"""Synchronous frozen-WT composition propagation and bounded output writers."""
from __future__ import annotations
import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import log_softmax

from .model import HYBRID,predict_change,logits
from .structure import iter_structure_predictions
from .composition import apply_intervention,next_inputs,state_diagnostics
from .binning import settings as binning_settings

EVENT_COLUMNS=['event','gene_id','gene_name','cell','cell_stage','cluster','pseudotime','branch_index',
    'probability_WT','probability_perturbed','delta_probability','delta_expression_only','delta_direct_composition',
    'delta_cascade','observed_event_molecules','cell_event_observed','status','target_relation_status',
    'flags','strength','targets','bootstrap_low','bootstrap_high','bootstrap_successful_fits']
CHAIN_COLUMNS=['gene_id','cell','chain_id','usage_WT','usage_perturbed','delta_usage','structure_scope','retained_chain_molecule_fraction']
AUDIT_COLUMNS=['gene','cell_index','status','rank','conflict']


def apply_dose(X,indices,strength):
    if not np.isfinite(strength) or not 0<=strength<=1:raise ValueError('Dose must be within [0,1]')
    if len(set(indices))!=len(indices):raise ValueError('Duplicate intervention targets')
    out=X.copy()
    if strength:out[:,indices]=np.log1p((1-strength)*np.expm1(X[:,indices]))
    return out


def resolve_targets(prior,targets):
    result=[]
    for target in targets:
        hits=[j for j,row in enumerate(zip(prior['factor_ids'],prior['gene_ids'],prior['symbols'])) if target in row]
        if len(hits)!=1:raise ValueError(f'Target {target} unavailable or ambiguous in prepared model')
        result.append(hits[0])
    return result


def propagate(models,X,frame,times,prior,targets,strength,c,*,bin_plan=None):
    if bin_plan is not None and binning_settings(c)['enabled']:
        frame,times=bin_plan.bin_frame,bin_plan.centers
    indices=resolve_targets(prior,targets)
    direct=apply_intervention(X,indices,strength,prior);current=direct.copy()
    source_ids={str(s['source_event_index']) for s in prior.get('isoform_states',[])}
    source_models={e:models[e] for e in sorted(source_ids) if e in models and models[e]['status']=='ok'}
    wt={e:np.asarray(m['WT_prediction']) for e,m in source_models.items()}
    history=[];stop='depth_budget';previous_state=None
    has_composition=any(m['design'].get('compositions') for m in models.values() if m['status']=='ok')
    for r in range(c['simulation']['max_rounds']):
        new={};logit_change=0.;maximum=0.
        changed_cols=set(np.flatnonzero(np.any(current!=previous_state,axis=0))) if previous_state is not None else set()
        if changed_cols:
            for m in models.values():
                if m['status']!='ok':continue
                used=m.get('state_input_indices',m['input_indices'])
                if changed_cols.isdisjoint(used):continue
                if m.get('model_kind')==HYBRID:
                    p=predict_change(m,X,current,frame,times,with_flags=False)[0]
                    p_prev=predict_change(m,X,previous_state,frame,times,with_flags=False)[0]
                    lp=np.log(np.maximum(p,1e-300))
                    old_lp=np.log(np.maximum(p_prev,1e-300))
                else:
                    lp=log_softmax(logits(m,current,frame,times),axis=1)
                    old_lp=log_softmax(logits(m,previous_state,frame,times),axis=1)
                logit_change=max(logit_change,float(np.max(np.abs(lp-old_lp))))
        for e,m in source_models.items():
            p,_,_=predict_change(m,X,current,frame,times,with_flags=False);new[e]=p
            maximum=max(maximum,float(np.max(np.abs(p-wt[e]))))
        n=len(prior['factor_ids'])
        entry={'round':r+1,'max_source_event_change':maximum,'max_adjacent_round_logit_change':logit_change,
            'changed_expression_inputs':[prior['factor_ids'][j] for j in range(n) if np.any(current[:,j]!=X[:,j])],
            'changed_isoform_inputs':[f['index'] for f in prior.get('composition_features',[]) if np.any(current[:,f['index']]!=X[:,f['index']])],
            'RNA_expression_unchanged_except_direct_proxy':True}
        history.append(entry)
        if logit_change>c['simulation'].get('divergence_limit',10.):
            entry['discarded_round']=True;current=previous_state;stop='divergence_limit_exceeded';break
        if strength==0:stop='no_op';break
        if not has_composition:stop='no_supported_isoform_inputs';break
        next_state,updates=next_inputs(X,direct,wt,new,source_models,prior,c)
        entry['isoform_updates']=updates
        if not np.isfinite(next_state).all():stop='nonfinite_isoform_state';break
        difference=next_state-current;delta=float(np.max(np.abs(difference)));norm=float(np.linalg.norm(difference))
        nonzero=np.count_nonzero(difference)
        rms=norm/float(np.sqrt(nonzero)) if nonzero else 0.
        entry.update(next_input_max_change=delta,next_input_update_norm=norm,next_input_rms_change=rms)
        if max(logit_change,delta)>c['simulation'].get('divergence_limit',10.):
            stop='divergence_limit_exceeded';break
        if delta<=c['simulation']['convergence_tolerance']:stop='no_new_isoform_input_change';break
        if r+1<c['simulation']['max_rounds']:previous_state=current;current=next_state
    return current,direct,history,stop


def relation_status(model,indices,changed,network_rows):
    if model is None or model['status']!='ok':return 'event_unmodeled'
    if set(indices)&set(model['input_indices']):return 'modeled_direct'
    if changed:return 'modeled_cascade_only'
    return 'not_selected_by_quality' if network_rows else 'no_prior_path'


def summarize_event(event,model,wt,new,components,N,frame,flags,status,network,targets,strength,stage_groups=None,*,unfitted_status=None,covered=None):
    rows=[];k=len(event.get('branches',[])) or (model or {}).get('branch_count',0)
    fit_status=model['status'] if model is not None else unfitted_status
    for stage,ix in (stage_groups if stage_groups is not None else frame.groupby('stage',sort=True).indices).items():
        ix=np.asarray(ix)
        ix_cov,ix_uncov=(ix,ix[:0]) if covered is None else (ix[covered[ix]],ix[~covered[ix]])
        n=N[ix_cov];good=n>0;total=float(n.sum())
        uncovered_cells=int(len(ix_uncov));uncovered_molecules=float(N[ix_uncov].sum())
        for branch in range(k):
            row={'event':event['event'],'gene_id':event['gene'],'cell_stage':stage,'branch_index':branch,
                 'targets':';'.join(sorted(targets)),'strength':strength,'target_relation_status':status,'event_fit_status':fit_status,
                 'network_best_rank':network.get('selection_rank'),'network_best_selection_score':network.get('selection_score'),
                 'network_best_oob_gain':network.get('oob_fit_gain'),'cells':len(ix),'observed_cells':int(good.sum()),
                 'observed_molecules':total,'uncovered_cells':uncovered_cells,'uncovered_molecules':uncovered_molecules,
                 'aggregation_scope':'WT_event_molecule_weighted; observed_cells; all_cells'}
            if wt is not None:
                d=new[ix_cov,branch]-wt[ix_cov,branch]
                mw=float(n@d/total) if total else np.nan
                row.update(probability_WT=float(n@wt[ix_cov,branch]/total) if total else np.nan,
                    probability_perturbed=float(n@new[ix_cov,branch]/total) if total else np.nan,
                    delta_probability=mw,delta_probability_molecule_weighted=mw,
                    delta_probability_observed_cells=float(d[good].mean()) if good.any() else np.nan,
                    delta_probability_all_cells=float(d.mean()) if len(d) else np.nan,
                    within_support_cells=sum(not flags[i] for i in ix_cov),
                    outside_support_cells=sum(any('extrapolation' in f or 'outside_WT' in f or 'unsupported' in f or 'p05' in f for f in flags[i]) for i in ix_cov),
                    weak_identifiability_cells=sum('individual_effect_unresolved' in flags[i] for i in ix_cov))
                for name,values in zip(['delta_expression_only_mw','delta_direct_composition_mw','delta_cascade_mw'],components):
                    row[name]=float(n@values[ix_cov,branch]/total) if total else np.nan
            else:
                row.update({key:np.nan for key in ['probability_WT','probability_perturbed','delta_probability',
                    'delta_probability_molecule_weighted','delta_probability_observed_cells','delta_probability_all_cells',
                    'delta_expression_only_mw','delta_direct_composition_mw','delta_cascade_mw']})
                row.update(within_support_cells=0,outside_support_cells=0,weak_identifiability_cells=0)
            rows.append(row)
    return rows


def simulate(models,X,frame,times,prior,structures,targets,strength,c,output=None,data=None,network=None,intervals=None,
             *,bin_plan=None):
    """Production streams per event/gene; the in-memory adapter supports small tests."""
    model_frame,model_times=frame,times
    covered=None
    if not binning_settings(c)['enabled']:
        bin_plan=None
    if bin_plan is not None:
        model_frame,model_times=bin_plan.bin_frame,bin_plan.centers
        covered=bin_plan.M_cell.any(axis=1)
        dominant=bin_plan.M_cell.argmax(axis=1) if len(model_frame) else np.zeros(len(frame),int)

    def cell_probability(p):
        if bin_plan is None:
            return p
        result=bin_plan.M_cell@p
        # Unsupported partitions carry a neutral, unchanged probability and flag.
        result[~covered]=1./p.shape[1]
        return result

    current,direct,history,stop=propagate(models,X,frame,times,prior,targets,strength,c,bin_plan=bin_plan)
    indices=resolve_targets(prior,targets);expression=apply_dose(X,indices,strength)
    min_delta_threshold=c['simulation'].get('min_delta_threshold',1e-4)
    write_unchanged_cells=c['simulation'].get('write_unchanged_cells',False)
    targets_label=';'.join(sorted(targets))
    event_rows=[];chain_rows=[];audit_rows=[];summary_rows=[];handles=[]
    if output is not None:
        output=Path(output);output.mkdir(parents=True,exist_ok=True)
        def gzip_writer(name,columns):
            h=gzip.open(output/name,'wt',compresslevel=1,newline='')
            handles.append(h)
            w=csv.writer(h,delimiter='\t');w.writerow(columns)
            return w
        ew=gzip_writer('events.tsv.gz',EVENT_COLUMNS);cw=gzip_writer('chains.tsv.gz',CHAIN_COLUMNS)
        h=(output/'structure_diagnostics.tsv').open('w',newline='');handles.append(h)
        aw=csv.DictWriter(h,fieldnames=AUDIT_COLUMNS,delimiter='\t',extrasaction='ignore');aw.writeheader()
    affected_genes=set();counts_written=0
    eligible_events=set(prior.get('eligible_events',[]))
    candidate_events=set(prior.get('events',{}))
    stage_groups=frame.groupby('stage',sort=True).indices
    event_keys=[str(i) for i in range(len(data['events']))] if data is not None else list(models)
    try:
        for e in event_keys:
            m=models.get(e);event=data['events'][int(e)] if data is not None else {'event':m['event_id'],'gene':m['gene_id']}
            if data is not None:
                j=int(e);N=np.asarray(data['C'][:,data['offsets'][j]:data['offsets'][j+1]].sum(1)).ravel()
            else:N=np.asarray(m.get('observed_denominators',np.zeros(len(frame))))
            relation_rows=(network or {}).get(e,[])
            best=min(relation_rows,key=lambda r:(r.get('selection_rank') is None,r.get('selection_rank') or float('inf'),-float(r.get('selection_score') or 0)),default={})
            unfitted_status=None
            if m is not None and m['status']=='ok':
                wt=cell_probability(np.asarray(m['WT_prediction']))
                new,_,flags=predict_change(m,X,current,model_frame,model_times)
                new=cell_probability(new)
                pe=cell_probability(predict_change(m,X,expression,model_frame,model_times,with_flags=False)[0])
                pdirect=cell_probability(predict_change(m,X,direct,model_frame,model_times,with_flags=False)[0])
                if bin_plan is not None:
                    flags=[flags[b] if supported else ['branch_unsupported'] for b,supported in zip(dominant,covered)]
                delta=new-wt
                components=(pe-wt,pdirect-pe,new-pdirect);changed=bool(np.any(delta!=0))
                status=relation_status(m,indices,changed,relation_rows)
                material_change=bool(np.max(np.abs(delta),initial=0.)>=min_delta_threshold) if min_delta_threshold>0 else changed
                if (status=='modeled_direct' and changed) or material_change:
                    affected_genes.add(m['gene_id'])
                write=write_unchanged_cells or output is None or status=='modeled_direct' or material_change
                if write:
                    ci=(intervals or {}).get(e)
                    event_fields=(m['event_id'],m['gene_id'],m['gene_name'])
                    successful_fits=ci['successful_fits'] if ci else 0
                    stages,clusters=frame.stage.array,frame.leiden.array
                    for i,(cell,stage,cluster,time,n,cell_flags) in enumerate(zip(frame.cell,stages,clusters,times,N,flags)):
                        observed=n>0
                        prediction_status='conditional_prediction' if observed else 'unobserved_event_prediction_not_imputation'
                        row_flags=';'.join(cell_flags)
                        for k in range(m['branch_count']):
                            record=(*event_fields,cell,stage,cluster,float(time),k,wt[i,k],new[i,k],delta[i,k],
                                    components[0][i,k],components[1][i,k],components[2][i,k],n,observed,
                                    prediction_status,status,row_flags,strength,targets_label,
                                    ci['low'][i,k] if ci else np.nan,ci['high'][i,k] if ci else np.nan,successful_fits)
                            if output is not None:ew.writerow(record)
                            else:event_rows.append(dict(zip(EVENT_COLUMNS,record)))
                            counts_written+=1
            else:
                wt=new=None;flags=None;components=()
                if m is None:
                    if len(event.get('branches',[]))>c['model'].get('max_event_branches',8):
                        unfitted_status='too_many_branches_unmodeled'
                    elif e in candidate_events:unfitted_status='not_selected_smoke_limit'
                    elif int(e) in eligible_events:unfitted_status='no_prior_path'
                    else:unfitted_status='insufficient_event_counts'
                status='no_prior_path' if unfitted_status=='no_prior_path' else 'event_unmodeled'
            summary_rows.extend(summarize_event(event,m,wt,new,components,N,frame,flags,status,best,targets,strength,stage_groups,
                unfitted_status=unfitted_status,covered=covered))
        # Each worker gets only one gene's arrays. Limit submitted results in the iterator.
        def gene_jobs():
            for gene in sorted(affected_genes):
                if gene not in structures:continue
                structure=structures[gene];wt={};new={}
                for e in structure['event_branch_maps']:
                    if e not in models or models[e]['status']!='ok':continue
                    m=models[e];wt[e]=cell_probability(np.asarray(m['WT_prediction']))
                    new[e]=cell_probability(predict_change(m,X,current,model_frame,model_times,with_flags=False)[0])
                yield gene,structure,wt,new,c
        workers=c.get('_runtime_threads',c.get('runtime',{}).get('threads',1)) if output is not None else 1
        for gene,result,audit in iter_structure_predictions(gene_jobs(),workers):
            for row in audit:
                if output is not None:aw.writerow(row)
                else:audit_rows.append(row)
            for i in np.flatnonzero(result['valid']):
                cell=frame.cell.iloc[i]
                for k,chain in enumerate(result['chain_ids']):
                    reference,perturbed=result['reference'][i,k],result['perturbed'][i,k]
                    row=(gene,cell,chain,reference,perturbed,perturbed-reference,result['scope'],result['retained_fraction'])
                    if output is not None:cw.writerow(row)
                    else:chain_rows.append(dict(zip(CHAIN_COLUMNS,row)))
    finally:
        for h in handles:h.close()
    if output is not None:pd.DataFrame(summary_rows).to_csv(output/'event_stage_summary.tsv',sep='\t',index=False)
    metadata={'targets':targets,'strength':strength,'rounds_executed':len(history),'stop_reason':stop,'history':history,
        'effect_type':'WT_learned_expression_and_local_isoform_composition_response','event_rows_written':counts_written,
        'affected_genes':len(affected_genes),'event_summary_rows':len(summary_rows),
        'isoform_state_diagnostics':state_diagnostics(prior,models),
        'propagation':'synchronous_frozen_WT_relative_composition_updates; no_fixed_time_lag',
        'manual_functional_weights_required':False,'no_S_completion':True,'WT_state_frozen':True,
        'no_gene_expression_feedback_inferred':True}
    return pd.DataFrame(event_rows,columns=EVENT_COLUMNS),pd.DataFrame(chain_rows,columns=CHAIN_COLUMNS),pd.DataFrame(audit_rows,columns=AUDIT_COLUMNS),metadata
