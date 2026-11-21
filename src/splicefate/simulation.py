"""Synchronous conditional perturbation; no numerical transmission over naked PPI."""
from __future__ import annotations
import numpy as np
import pandas as pd

from .model import predict_change
from .structure import structure_predictions
from .composition import apply_intervention,next_inputs,state_diagnostics


def apply_dose(X,indices,strength):
    if not np.isfinite(strength) or not 0<=strength<=1:raise ValueError('Dose must be within [0,1]')
    if len(set(indices))!=len(indices):raise ValueError('Duplicate intervention targets')
    out=X.copy()
    if strength:
        out[:,indices]=np.log1p((1-strength)*np.expm1(X[:,indices]))
    return out


def simulate(models,X,frame,times,prior,structures,targets,strength,c):
    target_indices=[]
    for t in targets:
        hits=[j for j,(f,g,s) in enumerate(zip(prior['factor_ids'],prior['gene_ids'],prior['symbols'])) if t in (f,g,s)]
        if len(hits)!=1:raise ValueError(f'Target {t} unavailable or ambiguous in prepared model')
        target_indices.append(hits[0])
    direct=apply_intervention(X,target_indices,strength,prior);current=direct.copy()
    event_wt={e:np.asarray(m['WT_prediction']) for e,m in models.items() if m['status']=='ok'}
    history=[];stop='depth_budget';event_new={};flags={};structure_result={};structure_audit=[]
    for r in range(c['simulation']['max_rounds']):
        event_new={};flags={}
        for e,m in models.items():
            if m['status']!='ok':continue
            p,delta,ff=predict_change(m,X,current,frame,times)
            event_new[e]=p;flags[e]=ff
        n_factors=len(prior['factor_ids'])
        history.append({'round':r+1,'max_event_change':max((float(np.max(np.abs(event_new[e]-event_wt[e]))) for e in event_new),default=0.),
            'changed_expression_inputs':[prior['factor_ids'][j] for j in range(n_factors) if np.any(current[:,j]!=X[:,j])],
            'changed_isoform_inputs':[f['index'] for f in prior.get('composition_features',[]) if np.any(current[:,f['index']]!=X[:,f['index']])],
            'RNA_expression_unchanged_except_direct_proxy':True})
        if strength==0:stop='no_op';break
        if not any(m['design'].get('compositions') for m in models.values() if m['status']=='ok'):
            stop='no_supported_isoform_inputs';break
        next_state,updates=next_inputs(X,direct,event_wt,event_new,models,prior)
        history[-1]['isoform_updates']=updates
        if np.any(~np.isfinite(next_state)):
            stop='nonfinite_isoform_state';break
        delta=float(np.max(np.abs(next_state-current)));history[-1]['next_input_max_change']=delta
        if delta<=c['simulation']['convergence_tolerance']:
            stop='no_new_isoform_input_change';break
        if r+1<c['simulation']['max_rounds']:current=next_state
    structure_result,structure_audit=structure_predictions(structures,event_wt,event_new,c)
    records=[]
    for e,m in models.items():
        if m['status']!='ok':
            records.append({'event':m['event_id'],'gene_id':m['gene_id'],'cell':'','branch_index':None,
                'status':m['status'],'delta_probability':None});continue
        for i,cell in enumerate(frame.cell):
            observed=float(m.get('observed_denominators',[0]*len(frame))[i])
            for k in range(m['branch_count']):
                records.append({'event':m['event_id'],'gene_id':m['gene_id'],'gene_name':m['gene_name'],
                    'cell':cell,'cell_stage':frame.stage.iloc[i],'cluster':frame.leiden.iloc[i],
                    'pseudotime':float(times[i]),'branch_index':k,
                    'probability_WT':float(event_wt[e][i,k]),'probability_perturbed':float(event_new[e][i,k]),
                    'delta_probability':float(event_new[e][i,k]-event_wt[e][i,k]),
                    'observed_event_molecules':observed,'cell_event_observed':observed>0,
                    'status':'conditional_prediction' if observed>0 else 'unobserved_event_prediction_not_imputation',
                    'flags':';'.join(flags[e][i]),'strength':strength,'targets':';'.join(targets)})
    chain_rows=[]
    for gene,s in structure_result.items():
        for i,cell in enumerate(frame.cell):
            if not s['valid'][i]:continue
            for k,chain in enumerate(s['chain_ids']):
                chain_rows.append({'gene_id':gene,'cell':cell,'chain_id':chain,'usage_WT':s['reference'][i,k],
                    'usage_perturbed':s['perturbed'][i,k],'delta_usage':s['perturbed'][i,k]-s['reference'][i,k],
                    'structure_scope':s['scope'],'retained_chain_molecule_fraction':s['retained_fraction']})
    metadata={'targets':targets,'strength':strength,'rounds_executed':len(history),'stop_reason':stop,'history':history,
        'effect_type':'WT_learned_expression_and_local_isoform_composition_response',
        'isoform_state_diagnostics':state_diagnostics(prior,models),
        'propagation':'synchronous_frozen_WT_relative_composition_updates; no_fixed_time_lag',
        'manual_functional_weights_required':False,
        'no_S_completion':True,'WT_state_frozen':True,'no_gene_expression_feedback_inferred':True}
    return pd.DataFrame(records),pd.DataFrame(chain_rows),pd.DataFrame(structure_audit),metadata
