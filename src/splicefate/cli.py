"""Prepare -> fit -> simulate; every stage binds inputs and source versions."""
from __future__ import annotations
import argparse
import json
import logging
from pathlib import Path
import time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .io import configuration,dump,load_data,provenance,digest

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
    p=Path(c['output']);saved=json.loads((p/'provenance.json').read_text())
    old=dict(saved['config']);new=dict(c)
    old.pop('_config_path',None);new.pop('_config_path',None)
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
    prior,X=build_priors(data,c)
    from .composition import build_states
    prior,X=build_states(data,prior,X,c)
    frame.to_csv(p/'cells.tsv',sep='\t',index=False)
    np.savez_compressed(p/'prepared_arrays.npz',times=times,pcs=pcs,X=X)
    dump(p/'priors.json',prior);dump(p/'trajectory_diagnostics.json',meta)
    summary={'status':'prepared','cells':len(frame),'events_total':len(data['events']),
             'events_count_eligible':len(prior['eligible_events']),'events_with_candidate_path':len(prior['events']),
             'factors':len(prior['factor_ids']),'seed_factor_ids':prior['seed_factor_ids'],
             'isoform_states':len(prior.get('isoform_states',[])),
             'isoform_contrasts':len(prior.get('composition_features',[])),
             'S_count_unit':data.get('S_count_unit','integer_umi'),
             'no_test_cells':True,'no_S_completion':True}
    manifest=provenance(c,input_files(c))
    artifacts=['cells.tsv','prepared_arrays.npz','priors.json','trajectory_diagnostics.json']
    if c['trajectory'].get('root_method')=='cytotrace':artifacts.append('cytotrace/scores.tsv')
    manifest['prepared_artifacts']={n:digest(p/n) for n in artifacts}
    dump(p/'provenance.json',manifest)
    dump(p/'prepare_summary.json',summary)
    LOG.info('Prepared %s',summary)


def fitted_inputs(c):
    p=check_prepared(c);data=load_data(c)
    frame=pd.read_csv(p/'cells.tsv',sep='\t',dtype={'cell':str,'batch':str,'stage':str,'branch':str,'leiden':str})
    if frame.cell.tolist()!=data['cells'].cell.tolist():raise ValueError('Prepared cell axis differs')
    arrays=np.load(p/'prepared_arrays.npz');prior=json.loads((p/'priors.json').read_text())
    return p,data,frame,arrays['times'],prior,arrays['X']


def fit(c,event_limit=None):
    from .model import fit_models
    from .structure import load_structures
    p,data,frame,times,prior,X=fitted_inputs(c)
    LOG.info('Fitting joint conditional count models; event_limit=%s',event_limit)
    models,audit=fit_models(data,frame,times,prior,X,c,event_limit)
    for e,m in models.items():
        m['observed_denominators']=np.asarray(data['C'][:,data['offsets'][int(e)]:data['offsets'][int(e)+1]].sum(1)).ravel().tolist()
    dump(p/'models.json',models);audit.to_csv(p/'fit_diagnostics.tsv',sep='\t',index=False)
    weights=[]
    for e,m in models.items():
        if m['status']!='ok':continue
        coef=np.asarray(m['coefficients'])
        for f in m['design'].get('compositions',[]):
            for row in range(f['width']):
                for k in range(m['branch_count']-1):
                    weights.append({'target_event':m['event_id'],'target_gene':m['gene_id'],
                        'factor_id':f['factor_id'],'source_event':f['source_event'],
                        'source_branch':f['branch_index'],'source_reference_branch':f['reference_branch'],
                        'target_branch':k,'target_reference_branch':m['branch_count']-1,
                        'term':'main' if row==0 else 'time_interaction',
                        'coefficient':float(coef[f['start']+row,k]),
                        'input_scale':f['scale'],'joint_observed_cells':f['joint_observed_cells']})
    pd.DataFrame(weights,columns=['target_event','target_gene','factor_id','source_event','source_branch',
                 'source_reference_branch','target_branch','target_reference_branch','term','coefficient',
                 'input_scale','joint_observed_cells']).to_csv(p/'isoform_response_weights.tsv',sep='\t',index=False)
    LOG.info('Mapping observed chains to event branches without completion')
    structures=load_structures(data,models,c);dump(p/'structures.json',structures)
    modeled={int(e) for e in models}
    status=[]
    for e,event in enumerate(data['events']):
        reason=models[str(e)]['status'] if e in modeled else ('not_selected_smoke_limit' if str(e) in prior['events'] else 'no_candidate_anchor_or_insufficient_counts')
        status.append({'event':event['event'],'gene':event['gene'],'status':reason})
    pd.DataFrame(status).to_csv(p/'all_event_status.tsv',sep='\t',index=False)
    summary={'status':'fitted','events_attempted':len(models),'events_ok':sum(m['status']=='ok' for m in models.values()),
        'events_total':len(data['events']),'event_limit':event_limit,'genes_with_structure':len(structures),
        'training_regime':'strict_WT_only_all_cells','fit_status_counts':audit.status.value_counts().to_dict() if len(audit) else {},
        'predictive_benchmark':None,'intervals':'not_estimated_unless_bootstrap_requested',
        'events_with_isoform_inputs':sum(bool(m.get('design',{}).get('compositions')) for m in models.values()),
        'manual_functional_weights_required':False,'S_count_unit':data.get('S_count_unit','integer_umi'),
        'models_sha256':digest(p/'models.json'),'structures_sha256':digest(p/'structures.json')}
    dump(p/'fit_summary.json',summary);LOG.info('%s',summary)


def bootstrap_intervals(c,data,frame,times,prior,X,models,targets,q,repeats):
    from .model import fit_event,logits
    from .simulation import simulate
    from scipy.special import softmax
    indices=[next(j for j,(f,g,s) in enumerate(zip(prior['factor_ids'],prior['gene_ids'],prior['symbols'])) if t in (f,g,s)) for t in targets]
    rng=np.random.default_rng(42);draws={}
    # Cells resampled within stage x batch, not falsely treated as embryo-level replication.
    strata=list(frame.groupby(['stage','batch'],sort=True).indices.values())
    samples=[np.concatenate([rng.choice(ix,len(ix),replace=True) for ix in strata]) for _ in range(repeats)]
    raw={e:data['C'][:,data['offsets'][int(e)]:data['offsets'][int(e)+1]].toarray() for e,m in models.items() if m['status']=='ok'}
    values={e:[] for e in raw}
    for sample in samples:
        fitted={}
        for e,C in raw.items():
            m=models[e]
            bm=fit_event(C[sample],X[sample],frame.iloc[sample].reset_index(drop=True),times[sample],
                         m['input_indices'],[tuple(p) for p in m['pairs']],c,m['design'].get('compositions',[]))
            if bm['status']!='ok':continue
            bm.update({k:m[k] for k in ['event_id','event_index','gene_id','gene_name','branch_count','observed_denominators']})
            bm['WT_prediction']=softmax(logits(bm,X,frame,times),axis=1).tolist();fitted[e]=bm
        if len(fitted)!=len(raw):continue
        prediction,_,_,_=simulate(fitted,X,frame,times,prior,{},targets,q,c)
        for e,m in fitted.items():
            v=prediction[prediction.event==m['event_id']].delta_probability.to_numpy().reshape(len(frame),m['branch_count'])
            values[e].append(v)
    for e,v in values.items():
        if v:draws[e]={'low':np.quantile(v,.025,axis=0),'high':np.quantile(v,.975,axis=0),'successful_fits':len(v)}
    return draws


def simulate_run(c,targets,q,bootstrap=0):
    from .simulation import simulate
    p,data,frame,times,prior,X=fitted_inputs(c)
    summary=json.loads((p/'fit_summary.json').read_text())
    for file,key in [('models.json','models_sha256'),('structures.json','structures_sha256')]:
        if digest(p/file)!=summary[key]:raise ValueError('Fit artifact changed: '+file)
    models=json.loads((p/'models.json').read_text());structures=json.loads((p/'structures.json').read_text())
    for target in targets:
        if not any(target in row for row in zip(prior['factor_ids'],prior['gene_ids'],prior['symbols'])):raise ValueError('Target not prepared: '+target)
    events,chains,audit,meta=simulate(models,X,frame,times,prior,structures,targets,q,c)
    if bootstrap:
        ci=bootstrap_intervals(c,data,frame,times,prior,X,models,targets,q,bootstrap)
        rowids={cell:i for i,cell in enumerate(frame.cell)};eventids={m['event_id']:e for e,m in models.items()}
        lo=[];hi=[];success=[]
        for r in events.itertuples():
            e=eventids[r.event]
            if e not in ci or r.cell not in rowids or pd.isna(r.branch_index):lo.append(np.nan);hi.append(np.nan);success.append(0);continue
            i=rowids[r.cell];k=int(r.branch_index);lo.append(ci[e]['low'][i,k]);hi.append(ci[e]['high'][i,k]);success.append(ci[e]['successful_fits'])
        events['bootstrap_low']=lo;events['bootstrap_high']=hi;events['bootstrap_successful_fits']=success
        meta['uncertainty']='WT cell bootstrap within stage x batch; complete learned composition cascade re-fitted'
    else:
        meta['uncertainty']='not_estimated; no fabricated confidence intervals'
    label='+'.join(targets).replace('/','_')+f'_q{q:g}';out=p/'simulations'/label;out.mkdir(parents=True,exist_ok=True)
    events.to_csv(out/'events.tsv.gz',sep='\t',index=False,compression='gzip')
    chains.to_csv(out/'chains.tsv.gz',sep='\t',index=False,compression='gzip')
    audit.to_csv(out/'structure_diagnostics.tsv',sep='\t',index=False)
    if 'probability_WT' in events:
        valid=events.dropna(subset=['delta_probability']).copy()
        valid['within_support']=valid['flags'].fillna('').eq('')
        valid['outside_support']=valid['flags'].fillna('').str.contains('extrapolation|outside_WT|unsupported')
        valid['weak_identifiability']=valid['flags'].fillna('').str.contains('unresolved')
        grouped=valid.groupby(['event','gene_id','cell_stage','branch_index'],as_index=False).agg(
            probability_WT=('probability_WT','mean'),probability_perturbed=('probability_perturbed','mean'),
            delta_probability=('delta_probability','mean'),cells=('cell','size'),observed_cells=('cell_event_observed','sum'),
            within_support_cells=('within_support','sum'),outside_support_cells=('outside_support','sum'),
            weak_identifiability_cells=('weak_identifiability','sum'))
        grouped['targets']=';'.join(sorted(targets));grouped['strength']=q
        grouped['aggregation_scope']='same_finite_cells_WT_and_perturbed; includes_flagged_conditional_predictions'
        grouped.to_csv(out/'event_stage_summary.tsv',sep='\t',index=False)
    dump(out/'simulation.json',meta);LOG.info('Saved %s; %s',out,meta['stop_reason'])


def evaluate(prediction,truth,output):
    """External fixed truth only; no fitting or hyperparameter selection."""
    p=pd.read_csv(prediction,sep='\t');t=pd.read_csv(truth,sep='\t')
    keys=['targets','strength','event','cell_stage','branch_index']
    if not set(keys)<=set(p) or not set(keys)<=set(t):raise ValueError('Benchmark inputs require target/dose/event/stage/branch keys: '+str(keys))
    if p.duplicated(keys).any() or t.duplicated(keys).any():raise ValueError('Benchmark keys must be unique')
    if 'delta_probability' not in t:raise ValueError('Truth needs delta_probability')
    m=p.merge(t,on=keys,suffixes=('_prediction','_truth'),validate='one_to_one').dropna(subset=['delta_probability_prediction','delta_probability_truth'])
    if not len(m):raise ValueError('No valid matched benchmark records')
    x=m.delta_probability_prediction.to_numpy();y=m.delta_probability_truth.to_numpy();rho=spearmanr(x,y).statistic
    report={'matched':len(m),'truth_rows':len(t),'prediction_rows':len(p),'MAE':float(np.mean(np.abs(x-y))),
        'zero_effect_MAE':float(np.mean(np.abs(y))),'spearman':float(rho) if np.isfinite(rho) else None,
        'truth_file_sha256':digest(truth),'prediction_file_sha256':digest(prediction),'tuning_performed':False}
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
    args=parser.parse_args();logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s')
    if getattr(args,'event_limit',None) is not None and args.event_limit<1:parser.error('--event-limit must be positive')
    if getattr(args,'bootstrap',0)<0:parser.error('--bootstrap must be nonnegative')
    if args.command=='evaluate':evaluate(args.prediction,args.truth,args.output);return
    c=configuration(args.config)
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
