"""Prepare -> fit -> simulate; every stage binds inputs and source versions."""
from __future__ import annotations
import argparse
import json
import logging
from pathlib import Path
import time
from uuid import uuid4
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .io import configuration,dump,load_data,provenance,digest
from .binning import BinPlan,fit_bin_plan,settings as binning_settings

LOG=logging.getLogger('splicefate')


def input_files(c):
    gp=Path(c['input']['genes']);ep=Path(c['input']['events'])
    files=[gp/n for n in ['gene_counts.npz','genes.tsv','cells.tsv','isoforms.tsv','isoform_counts.npz']]
    files+=[ep/n for n in ['branch_counts.npz','branch_offsets.npy','groups.json','genes.tsv','normalization_factors.tsv','cells.tsv']]
    files += [Path(c['input']['database']),Path(c['input']['fasta']),Path(c['input']['fasta']+'.fai')]
    if (gp/'quantification.json').is_file():files.append(gp/'quantification.json')
    if c['trajectory'].get('root_method')=='cytotrace':
        from .cytotrace import source_files
        files+=source_files(c)
    if c['input'].get('velocity_times'):files.append(Path(c['input']['velocity_times']))
    if c['input'].get('intron_counts'):
        files += [Path(c['input']['intron_counts'])/n for n in ['exon_counts.npz','intron_counts.npz','genes.tsv','cells.tsv','quantification.json']]
    return files


def check_prepared(c):
    p=Path(c['output']);manifest=p/'provenance.json'
    if not manifest.is_file():raise ValueError('Prepared provenance missing; prepare the run first')
    saved=json.loads(manifest.read_text())
    old=dict(saved['config']);new=dict(c)
    old.pop('_config_path',None);new.pop('_config_path',None)
    old.pop('_runtime_threads',None);new.pop('_runtime_threads',None)
    if old!=new:raise ValueError('Configuration differs from prepared run; use a new output or prepare again')
    for file,sha in saved['inputs'].items():
        if digest(file)!=sha:raise ValueError('Prepared input changed: '+file)
    for file,sha in saved['source_sha256'].items():
        if digest(file)!=sha:raise ValueError('Source changed since prepare; prepare a new run: '+file)
    for file,sha in saved.get('prepared_artifacts',{}).items():
        if digest(p/file)!=sha:raise ValueError('Prepared artifact changed: '+file)
    return p


def prepare(c):
    from .trajectory import prepare_trajectory
    from .priors import build_priors
    p=Path(c['output']);p.mkdir(parents=True,exist_ok=True)
    LOG.info('Loading full gene counts and source event support')
    data=load_data(c)
    LOG.info('Reference Leiden and splice-calibrated WT trajectory, all %d cells',len(data['cells']))
    frame,times,pcs,meta=prepare_trajectory(data,c)
    LOG.info('Native %s motifs, CLIP, protein partners and complex candidates',c['assembly'])
    prior,X_cell_factors=build_priors(data,c,output=p)
    from .composition import build_states
    bin_plan=fit_bin_plan(frame,times,c)
    X=X_cell_factors if bin_plan is None else bin_plan.W_agg@X_cell_factors
    prior,X=build_states(data,prior,X,c,bin_plan)
    frame.to_csv(p/'cells.tsv',sep='\t',index=False)
    arrays=dict(times=times,pcs=pcs,X=X)
    if bin_plan is not None:
        arrays.update(bin_plan.arrays(),X_cell_factors=X_cell_factors)
        dump(p/'binning_diagnostics.json',bin_plan.diagnostics())
    np.savez_compressed(p/'prepared_arrays.npz',**arrays)
    dump(p/'priors_meta.json',{k:v for k,v in prior.items() if k not in ('events','pairs')})
    dump(p/'trajectory_diagnostics.json',meta)
    summary={'status':'prepared','cells':len(frame),'events_total':len(data['events']),
             'events_count_eligible':len(prior['eligible_events']),'events_with_candidate_path':len(prior['events']),
             'factors':len(prior['factor_ids']),'seed_factor_ids':prior['seed_factor_ids'],
             'isoform_states':len(prior.get('isoform_states',[])),
             'isoform_contrasts':len(prior.get('composition_features',[])),
             'S_count_unit':data.get('S_count_unit','integer_umi'),
             'no_test_cells':True,'no_S_completion':True}
    manifest=provenance(c,input_files(c))
    manifest['prepare_generation']=uuid4().hex
    artifacts=['cells.tsv','prepared_arrays.npz','priors_meta.json','priors.duckdb','trajectory_diagnostics.json']
    if bin_plan is not None:artifacts.append('binning_diagnostics.json')
    if c['trajectory'].get('root_method')=='cytotrace':artifacts.append('cytotrace/scores.tsv')
    manifest['prepared_artifacts']={n:digest(p/n) for n in artifacts}
    dump(p/'provenance.json',manifest)
    dump(p/'prepare_summary.json',summary)
    LOG.info('Prepared %s',summary)


def fitted_inputs(c):
    p=check_prepared(c);data=load_data(c)
    frame=pd.read_csv(p/'cells.tsv',sep='\t',dtype={'cell':str,'batch':str,'stage':str,'branch':str,'leiden':str})
    if frame.cell.tolist()!=data['cells'].cell.tolist():raise ValueError('Prepared cell axis differs')
    from .storage import load_priors
    prior=load_priors(p)
    with np.load(p/'prepared_arrays.npz') as arrays:
        X,times=arrays['X'],arrays['times']
        bin_plan=None
        if binning_settings(c)['enabled']:
            bin_plan=BinPlan.load(arrays,json.loads((p/'binning_diagnostics.json').read_text()))
            data['X_cell_factors']=arrays['X_cell_factors']
    return p,data,frame,times,prior,X,bin_plan


def fit(c,event_limit=None):
    from .model import fit_models,hybrid_selection
    from .structure import load_structures
    p,data,frame,times,prior,X,*plans=fitted_inputs(c)
    bin_plan=plans[0] if plans else None
    manifest=p/'provenance.json'
    if not plans and not manifest.exists():
        # Legacy in-memory providers supply and verify the fitted inputs themselves.
        dump(manifest,{**provenance(c,()),'prepare_generation':uuid4().hex,'input_scope':'in_memory_fitted_inputs'})
    prepared_provenance_sha256=digest(manifest)
    LOG.info('Learning and evaluating all candidate single relations; event_limit=%s',event_limit)
    network_path=p/'regulatory_network.duckdb'
    models,audit=fit_models(data,frame,times,prior,X,c,event_limit,network_path=network_path,model_root=p,bin_plan=bin_plan)
    import duckdb
    with duckdb.connect(str(network_path),read_only=True) as network:
        screening_metadata=dict(network.execute('''SELECT key,value FROM network_metadata
            WHERE key IN ('screening_mode','training_evidence','oob_scope')''').fetchall())
        network.execute('''COPY (SELECT event_id,gene_id,factor_id,gene_symbol,status,
            response_score,supported_adjacent_pairs,quality_status,fit_gain,partial_r2,penalized_fit_gain,stability,
            oob_fit_gain,oob_positive_gain_fraction,bootstrap_successes,bootstrap_positive_gain_fraction,bootstrap_effect_agreement,
            bootstrap_gain_mean,normalized_fit_gain,normalized_response,selection_score,
            selection_rank,selected,WT_weight
            FROM single_relations ORDER BY event_index,selection_rank NULLS LAST,factor_id)
            TO ? (HEADER,DELIMITER '\t',COMPRESSION GZIP)''',[str(p/'regulatory_network.tsv.gz')])
    extra=[{'event':event['event'],'gene':event['gene'],'status':'too_many_branches_unmodeled','network_candidates':0,'factors':0,'interactions':0,'isoform_inputs':0,'counted_cells':int((np.asarray(data['C'][:,data['offsets'][e]:data['offsets'][e+1]].sum(1)).ravel()>0).sum()),'molecules':float(data['C'][:,data['offsets'][e]:data['offsets'][e+1]].sum()),'message':''} for e,event in enumerate(data['events']) if str(e) not in models and len(event.get('branches',[]))>c['model']['max_event_branches']]
    if extra:audit=pd.concat([audit,pd.DataFrame(extra)],ignore_index=True)
    audit.to_csv(p/'fit_diagnostics.tsv',sep='\t',index=False)
    weights=[];single_counts=[];events_with_groups=0;events_with_isoform_inputs=0
    def observed_cells(e):
        return np.asarray(data['C'][:,data['offsets'][e]:data['offsets'][e+1]].sum(1)).ravel()>0
    source_observed={s['source_event_index']:observed_cells(s['source_event_index'])
                     for s in prior.get('isoform_states',[])} if bin_plan is not None else {}
    for m in models.values():
        if m['status']!='ok':continue
        singles=m.get('singles')
        if singles is not None:
            single_counts.append(len(singles))
            events_with_groups+=bool(m['groups'])
        has_composition=bool(m['design'].get('compositions'))
        events_with_isoform_inputs+=has_composition
        target_observed=observed_cells(m['event_index']) if bin_plan is not None and has_composition else None
        for part in singles if singles is not None else [m]:
            coef=np.asarray(part['coefficients'])
            for f in part['design'].get('compositions',[]):
                joint=f.get('joint_observed_cells') if target_observed is None else int(
                    np.count_nonzero(target_observed & source_observed[f['source_event_index']]))
                for row in range(f['width']):
                    for k in range(m['branch_count']-1):
                        weights.append({'target_event':m['event_id'],'target_gene':m['gene_id'],
                            'factor_id':f['factor_id'],'source_event':f['source_event'],
                            'source_branch':f['branch_index'],'source_reference_branch':f['reference_branch'],
                            'target_branch':k,'target_reference_branch':m['branch_count']-1,
                            'term':'main' if row==0 else 'time_interaction',
                            'coefficient':float(coef[f['start']+row,k]),
                            'input_scale':f['scale'],'joint_observed_cells':joint})
    pd.DataFrame(weights,columns=['target_event','target_gene','factor_id','source_event','source_branch',
                 'source_reference_branch','target_branch','target_reference_branch','term','coefficient',
                 'input_scale','joint_observed_cells']).to_csv(p/'isoform_response_weights.tsv',sep='\t',index=False)
    LOG.info('Mapping observed chains to event branches without completion')
    structures=load_structures(data,models,c)
    from .storage import save_structures
    save_structures(p,structures,len(frame))
    modeled={int(e) for e in models}
    status=[]
    for e,event in enumerate(data['events']):
        reason=models[str(e)]['status'] if e in modeled else ('not_selected_smoke_limit' if str(e) in prior['events'] else 'no_candidate_anchor_or_insufficient_counts')
        if len(event.get('branches',[]))>c['model']['max_event_branches']:reason='too_many_branches_unmodeled'
        status.append({'event':event['event'],'gene':event['gene'],'status':reason})
    pd.DataFrame(status).to_csv(p/'all_event_status.tsv',sep='\t',index=False)
    summary={'status':'fitted','events_attempted':len(models),'events_ok':sum(m['status']=='ok' for m in models.values()),
        'events_total':len(data['events']),'event_limit':event_limit,'genes_with_structure':len(structures),
        'training_regime':'strict_WT_only_all_cells','fit_status_counts':audit.status.value_counts().to_dict() if len(audit) else {},
        'predictive_benchmark':None,'intervals':'not_estimated_unless_bootstrap_requested',
        'events_with_isoform_inputs':events_with_isoform_inputs,
        'binning':bin_plan.summary() if bin_plan is not None else {'enabled':False},
        'selection_mode':c['model'].get('selection_mode','all_passing_singles_hybrid'),
        'events_with_hybrid':len(single_counts),'events_with_groups':events_with_groups,
        'avg_retained_singles':float(np.mean(single_counts)) if single_counts else 0.,
        'manual_functional_weights_required':False,'S_count_unit':data.get('S_count_unit','integer_umi'),
        'artifacts':{n:digest(p/n) for n in ['models.duckdb','model_arrays.npz','structures_meta.json','structure_counts.npz']},
        'prepared_provenance_sha256':prepared_provenance_sha256,**screening_metadata,
        'network':json.loads(network_path.with_suffix('.summary.json').read_text()),
        'network_sha256':digest(network_path),
        'max_factors_after_network_learning':c['model'].get('max_single_factors') if hybrid_selection(c) else c['model']['max_factors']}
    dump(p/'fit_summary.json',summary);LOG.info('%s',summary)


class IntervalStore:
    def __init__(self):
        import tempfile
        self.temp=tempfile.TemporaryDirectory(prefix='splicefate_intervals_');self.root=Path(self.temp.name)
    def get(self,e,default=None):
        path=self.root/(str(e)+'.npz')
        if not path.exists():return default
        with np.load(path) as z:return {'low':z['low'],'high':z['high'],'successful_fits':int(z['success'])}
    def close(self):self.temp.cleanup()


def bootstrap_intervals(c,data,frame,times,prior,X,models,targets,q,repeats,bin_plan=None):
    from .composition import build_composition_view
    from .model import (HYBRID,fit_event,logits,predict_change,assemble_hybrid,
                        observed_rows,_fit_group_counts,refit_singles,single_sum_probability)
    from .network import prepare_bin_bootstrap,event_counts,compact_model
    from .simulation import propagate,resolve_targets
    from scipy.special import softmax
    import tempfile
    import shutil
    from .storage import ModelWriter
    rng=np.random.default_rng(42);store=IntervalStore()
    states=prior.get('isoform_states',[])
    source_ids={str(state['source_event_index']) for state in states}
    source_ids={e for e in source_ids if e in models and models[e]['status']=='ok'}
    resolve_targets(prior,targets)
    dependents={}
    for e,wt_model in models.items():
        for feature in wt_model.get('design',{}).get('compositions',[]):
            source=str(feature['source_event_index'])
            dependents.setdefault(source,set()).add(str(e))
    if bin_plan is None:
        strata=list(frame.groupby(['branch','stage','batch'],sort=True).indices.values())
        views=[np.concatenate([rng.choice(ix,len(ix),replace=True) for ix in strata]) for _ in range(repeats)]
    else:
        views=prepare_bin_bootstrap(bin_plan,frame,prior,data,data['X_cell_factors'],c,repeats=repeats,seed=42)
        frame,times=bin_plan.bin_frame,bin_plan.centers
        eval_counts={state['source_event_index']:np.asarray(state['counts']) for state in states}
        eval_factors=X[:,:len(prior['factor_ids'])]
    with tempfile.TemporaryDirectory(prefix='splicefate_draws_') as tmp:
        root=Path(tmp)
        for repeat,view in enumerate(views):
            X_eval,prior_draw=X,prior
            if bin_plan is None:
                Xt,ft,tt=X[view],frame.iloc[view].reset_index(drop=True),times[view]
            else:
                Xt,K_tr=view.X_tr,view.K_tr
                ft,tt=frame.assign(valid=view.valid_tr),times
                if states:
                    X_eval=build_composition_view(eval_factors,eval_counts,prior,c,pooled=view.pooled)
                    prior_draw=dict(prior,isoform_states=[
                        dict(state,pooled_composition=view.pooled[state['source_event_index']].tolist()) for state in states])
            repeat_root=root/f'refit_{repeat}';writer=ModelWriter(repeat_root)
            try:
                for e,m in models.items():
                    if m['status']!='ok':continue
                    C=event_counts(data,int(e))
                    Ct=C[view] if bin_plan is None else K_tr@C
                    observed=observed_rows(C,frame) if bin_plan is None else observed_rows(bin_plan.K_eff@C,frame)
                    if m.get('model_kind')==HYBRID:
                        singles=refit_singles(m['singles'],Ct,Xt,ft,tt,c)
                        if singles is None:continue
                        good=observed_rows(Ct,ft)
                        groups=[]
                        if m['groups']:
                            X_good,C_good,fr_good,t_good=Xt[good],Ct[good],ft.iloc[good],tt[good]
                            baseline=single_sum_probability(singles,X_good,fr_good,t_good)
                            for group in m['groups']:
                                g,_=_fit_group_counts(C_good,X_good,fr_good,t_good,group['input_indices'],c,baseline)
                                if g is None:break
                                groups.append(g)
                            if len(groups)!=len(m['groups']):continue
                        compact=[dict(compact_model(single,good),weight=single['weight']) for single in singles]
                        bm=assemble_hybrid(compact,groups,[],X_eval,frame,times,observed,c)
                    else:
                        bm=fit_event(Ct,Xt,ft,tt,m['input_indices'],[tuple(p) for p in m['pairs']],c,m['design'].get('compositions',[]))
                        if bm['status']!='ok':continue
                        bm['WT_prediction']=softmax(logits(bm,X_eval,frame,times),axis=1)
                    bm.update({k:m[k] for k in ['event_id','event_index','gene_id','gene_name','branch_count']})
                    writer.append(e,bm)
            except BaseException:
                writer.close()
                raise
            fitted=writer.finish()
            failed_sources=source_ids-fitted.keys()
            affected_by_failure=set()
            pending=list(failed_sources)
            while pending:
                source=pending.pop()
                for dependent in dependents.get(source,()):
                    if dependent not in affected_by_failure:
                        affected_by_failure.add(dependent)
                        pending.append(dependent)
            current,_,_,_=propagate(fitted,X_eval,frame,times,prior_draw,targets,q,c)
            for e,m in fitted.items():
                if e in affected_by_failure:continue
                delta=predict_change(m,X_eval,current,frame,times,with_flags=False)[1]
                if bin_plan is not None:delta=bin_plan.M_cell@delta
                np.save(root/f'{e}_{repeat}.npy',delta)
            del fitted
            shutil.rmtree(repeat_root)
        for e,m in models.items():
            files=list(root.glob(str(e)+'_*.npy'))
            if not files:continue
            values=np.stack([np.load(file) for file in files])
            np.savez(store.root/(str(e)+'.npz'),low=np.quantile(values,.025,axis=0),high=np.quantile(values,.975,axis=0),success=len(files))
    return store


def simulate_run(c,targets,q,bootstrap=0):
    from .simulation import simulate
    from .storage import ModelStore,StructureStore
    import duckdb
    p,data,frame,times,prior,X,*plans=fitted_inputs(c)
    bin_plan=plans[0] if plans else None
    summary=json.loads((p/'fit_summary.json').read_text())
    manifest=p/'provenance.json'
    if not manifest.is_file() or summary.get('prepared_provenance_sha256')!=digest(manifest):
        raise ValueError('Prepared provenance differs from fit; fit the current prepared run')
    network_path=p/'regulatory_network.duckdb'
    if not network_path.is_file() or summary.get('network_sha256')!=digest(network_path):
        raise ValueError('Fit artifact changed: regulatory_network.duckdb')
    for file,sha in summary['artifacts'].items():
        if digest(p/file)!=sha:raise ValueError('Fit artifact changed: '+file)
    models=ModelStore(p);structures=StructureStore(p)
    label='+'.join(sorted(targets)).replace('/','_')+f'_q{q:g}'
    out=p/'simulations'/label;out.mkdir(parents=True,exist_ok=True)
    wanted={f for target in targets for f,g,s in zip(prior['factor_ids'],prior['gene_ids'],prior['symbols']) if target in (f,g,s)}
    with duckdb.connect(str(network_path),read_only=True) as db:
        rows=db.execute('SELECT event_index,factor_id,selection_rank,selection_score,oob_fit_gain FROM single_relations WHERE factor_id IN (SELECT unnest(?))',[sorted(wanted)]).fetchall()
    network={}
    for e,f,rank,score,gain in rows:
        network.setdefault(str(e),[]).append({'factor_id':f,'selection_rank':rank,'selection_score':score,'oob_fit_gain':gain})
    ci=bootstrap_intervals(c,data,frame,times,prior,X,models,targets,q,bootstrap,bin_plan) if bootstrap else None
    try:
        _,_,_,meta=simulate(models,X,frame,times,prior,structures,targets,q,c,output=out,data=data,network=network,intervals=ci,bin_plan=bin_plan)
    finally:
        if ci:ci.close()
    meta['uncertainty']='WT stratified cell bootstrap; event-specific successful refits' if bootstrap else 'not_estimated'
    dump(out/'simulation.json',meta);LOG.info('Saved %s; %s',out,meta['stop_reason'])


def evaluate(prediction,truth,output,metric='delta_probability_molecule_weighted'):
    """External fixed truth, never used for fitting/selection."""
    p=pd.read_csv(prediction,sep='\t');t=pd.read_csv(truth,sep='\t')
    keys=['targets','strength','event','cell_stage','branch_index']
    if not set(keys)<=set(p) or not set(keys)<=set(t):raise ValueError('Benchmark inputs require target/dose/event/stage/branch keys: '+str(keys))
    if p.duplicated(keys).any() or t.duplicated(keys).any():raise ValueError('Benchmark keys must be unique')
    truth_column='delta_probability_truth' if 'delta_probability_truth' in t else 'delta_probability'
    if truth_column not in t:raise ValueError('Truth needs delta_probability_truth')
    # Legacy prediction adapters can explicitly request delta_probability.
    if metric not in p:raise ValueError('Prediction missing requested metric: '+metric)
    keep=keys+[metric]+(['target_relation_status'] if 'target_relation_status' in p else [])
    left=p[keep].rename(columns={metric:'prediction_value'})
    right=t[keys+[truth_column]].rename(columns={truth_column:'truth_value'})
    m=left.merge(right,on=keys,validate='one_to_one').dropna(subset=['prediction_value','truth_value'])
    if not len(m):raise ValueError('No valid matched benchmark records')
    def metrics(rows):
        x=rows.prediction_value.to_numpy();y=rows.truth_value.to_numpy()
        rho=spearmanr(x,y).statistic if len(x)>1 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
        return {'matched':len(rows),'MAE':float(np.mean(np.abs(x-y))),'zero_effect_MAE':float(np.mean(np.abs(y))),
                'spearman':float(rho) if np.isfinite(rho) else None}
    report={**metrics(m),'metric':metric,'truth_rows':len(t),'prediction_rows':len(p),
        'truth_file_sha256':digest(truth),'prediction_file_sha256':digest(prediction),'tuning_performed':False}
    if 'target_relation_status' in m:
        report['by_target_relation_status']={str(k):metrics(group) for k,group in m.groupby('target_relation_status')}
        report['prediction_status_counts']=p.target_relation_status.value_counts().to_dict()
    report['by_target']={str(k):metrics(group) for k,group in m.groupby('targets')}
    dump(output,report)


def main():
    parser=argparse.ArgumentParser(description='WT-conditioned splicing perturbation with RNA and protein priors')
    sub=parser.add_subparsers(dest='command',required=True)
    for name in ['prepare','fit','simulate','run']:
        p=sub.add_parser(name);p.add_argument('--config',required=True);p.add_argument('--output')
        if name in ('fit','run'):p.add_argument('--event-limit',type=int,help='Explicit smoke run scope; never reported as a full fit')
        if name in ('simulate','run'):
            p.add_argument('--targets',nargs='+',required=True);p.add_argument('--strength',type=float,default=1.)
            p.add_argument('--bootstrap',type=int,default=0)
    p=sub.add_parser('evaluate');p.add_argument('--prediction',required=True);p.add_argument('--truth',required=True);p.add_argument('--output',required=True)
    p.add_argument('--metric',choices=['delta_probability_molecule_weighted','delta_probability_observed_cells','delta_probability_all_cells','delta_probability'],default='delta_probability_molecule_weighted')
    args=parser.parse_args();logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s')
    if getattr(args,'event_limit',None) is not None and args.event_limit<1:parser.error('--event-limit must be positive')
    if getattr(args,'bootstrap',0)<0:parser.error('--bootstrap must be nonnegative')
    if args.command=='evaluate':evaluate(args.prediction,args.truth,args.output,args.metric);return
    c=configuration(args.config)
    c['_runtime_threads']=c.get('runtime',{}).get('threads',1)
    if args.output:c['output']=str(Path(args.output).resolve())
    out=Path(c['output']);out.mkdir(parents=True,exist_ok=True)
    status={'status':'running','command':args.command,'started_at':time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    dump(out/'status.json',status)
    try:
        if args.command in ('prepare','run'):prepare(c)
        if args.command in ('fit','run'):fit(c,args.event_limit)
        if args.command in ('simulate','run'):simulate_run(c,args.targets,args.strength,args.bootstrap)
        status.update(status='complete',finished_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
    except BaseException as exc:
        status.update(status='failed',error=str(exc));dump(out/'status.json',status);raise
    dump(out/'status.json',status)


if __name__=='__main__':main()
