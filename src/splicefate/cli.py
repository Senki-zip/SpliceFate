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
    prior,X=build_priors(data,c,output=p)
    from .composition import build_states
    prior,X=build_states(data,prior,X,c)
    frame.to_csv(p/'cells.tsv',sep='\t',index=False)
    np.savez_compressed(p/'prepared_arrays.npz',times=times,pcs=pcs,X=X)
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
    artifacts=['cells.tsv','prepared_arrays.npz','priors_meta.json','priors.duckdb','trajectory_diagnostics.json']
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
    arrays=np.load(p/'prepared_arrays.npz');prior=load_priors(p)
    return p,data,frame,arrays['times'],prior,arrays['X']


def fit(c,event_limit=None):
    from .model import fit_models
    from .structure import load_structures
    p,data,frame,times,prior,X=fitted_inputs(c)
    LOG.info('Learning all candidate single relations, then selecting joint inputs; event_limit=%s',event_limit)
    network_path=p/'regulatory_network.duckdb'
    models,audit=fit_models(data,frame,times,prior,X,c,event_limit,network_path=network_path,model_root=p)
    import duckdb
    with duckdb.connect(str(network_path),read_only=True) as network:
        network.execute('''COPY (SELECT event_id,gene_id,factor_id,gene_symbol,status,
            response_score,supported_adjacent_pairs,quality_status,fit_gain,stability,
            oob_fit_gain,oob_positive_gain_fraction,bootstrap_successes,bootstrap_positive_gain_fraction,bootstrap_effect_agreement,
            bootstrap_gain_mean,normalized_fit_gain,normalized_response,selection_score,
            selection_rank,selected,WT_weight
            FROM single_relations ORDER BY event_index,selection_rank NULLS LAST,factor_id)
            TO ? (HEADER,DELIMITER '\t',COMPRESSION GZIP)''',[str(p/'regulatory_network.tsv.gz')])
    extra=[{'event':event['event'],'gene':event['gene'],'status':'too_many_branches_unmodeled','network_candidates':0,'factors':0,'interactions':0,'isoform_inputs':0,'counted_cells':int((np.asarray(data['C'][:,data['offsets'][e]:data['offsets'][e+1]].sum(1)).ravel()>0).sum()),'molecules':float(data['C'][:,data['offsets'][e]:data['offsets'][e+1]].sum()),'message':''} for e,event in enumerate(data['events']) if str(e) not in models and len(event.get('branches',[]))>c['model']['max_event_branches']]
    if extra:audit=pd.concat([audit,pd.DataFrame(extra)],ignore_index=True)
    audit.to_csv(p/'fit_diagnostics.tsv',sep='\t',index=False)
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
        'events_with_isoform_inputs':sum(bool(m.get('design',{}).get('compositions')) for m in models.values()),
        'manual_functional_weights_required':False,'S_count_unit':data.get('S_count_unit','integer_umi'),
        'artifacts':{n:digest(p/n) for n in ['models.duckdb','model_arrays.npz','structures_meta.json','structure_counts.npz']},
        'network':json.loads(network_path.with_suffix('.summary.json').read_text()),
        'network_sha256':digest(network_path),
        'max_factors_after_network_learning':c['model']['max_factors']}
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


def bootstrap_intervals(c,data,frame,times,prior,X,models,targets,q,repeats):
    from .model import fit_event,logits,predict_change
    from .simulation import propagate
    from scipy.special import softmax
    import tempfile
    import shutil
    from .storage import ModelWriter
    rng=np.random.default_rng(42);store=IntervalStore()
    strata=list(frame.groupby(['branch','stage','batch'],sort=True).indices.values())
    with tempfile.TemporaryDirectory(prefix='splicefate_draws_') as tmp:
        root=Path(tmp)
        for repeat in range(repeats):
            sample=np.concatenate([rng.choice(ix,len(ix),replace=True) for ix in strata])
            repeat_root=root/f'refit_{repeat}';writer=ModelWriter(repeat_root)
            for e,m in models.items():
                if m['status']!='ok':continue
                j=int(e);C=data['C'][:,data['offsets'][j]:data['offsets'][j+1]].toarray()
                bm=fit_event(C[sample],X[sample],frame.iloc[sample].reset_index(drop=True),times[sample],
                    m['input_indices'],[tuple(p) for p in m['pairs']],c,m['design'].get('compositions',[]))
                if bm['status']!='ok':continue
                bm.update({k:m[k] for k in ['event_id','event_index','gene_id','gene_name','branch_count']})
                bm['WT_prediction']=softmax(logits(bm,X,frame,times),axis=1);writer.append(e,bm)
            fitted=writer.finish()
            current,_,_,_=propagate(fitted,X,frame,times,prior,targets,q,c)
            for e,m in fitted.items():
                np.save(root/f'{e}_{repeat}.npy',predict_change(m,X,current,frame,times,with_flags=False)[1])
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
    p,data,frame,times,prior,X=fitted_inputs(c)
    summary=json.loads((p/'fit_summary.json').read_text())
    for file,sha in summary['artifacts'].items():
        if digest(p/file)!=sha:raise ValueError('Fit artifact changed: '+file)
    models=ModelStore(p,data);structures=StructureStore(p)
    label='+'.join(sorted(targets)).replace('/','_')+f'_q{q:g}'
    out=p/'simulations'/label;out.mkdir(parents=True,exist_ok=True)
    wanted={f for target in targets for f,g,s in zip(prior['factor_ids'],prior['gene_ids'],prior['symbols']) if target in (f,g,s)}
    with duckdb.connect(str(p/'regulatory_network.duckdb'),read_only=True) as db:
        rows=db.execute('SELECT event_index,factor_id,selection_rank,selection_score,oob_fit_gain FROM single_relations WHERE factor_id IN (SELECT unnest(?))',[sorted(wanted)]).fetchall()
    network={}
    for e,f,rank,score,gain in rows:
        network.setdefault(str(e),[]).append({'factor_id':f,'selection_rank':rank,'selection_score':score,'oob_fit_gain':gain})
    ci=bootstrap_intervals(c,data,frame,times,prior,X,models,targets,q,bootstrap) if bootstrap else None
    try:
        _,_,_,meta=simulate(models,X,frame,times,prior,structures,targets,q,c,output=out,data=data,network=network,intervals=ci)
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
