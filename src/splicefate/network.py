"""Learn all prior-supported WT single relations before selecting joint inputs."""
from __future__ import annotations

import json
import hashlib
import logging
import os
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import numpy as np
from scipy.special import softmax, gammaln


LOG = logging.getLogger('splicefate')
POLICY = 'all_single_relations_then_weighted_quality'
SCORE = 'mean_adjacent_WT_factor_response_total_variation'


def selection_settings(c):
    return {'method':POLICY, 'min_response_score':0., 'min_oob_fit_gain':0.,
            'min_positive_gain_fraction':.6,'min_stability':.3,
            'bootstrap_repeats':10, 'bootstrap_seed':42,
            'weights':{'fit_gain':.45, 'stability':.45, 'response':.10},
            **c['model'].get('network_selection', {})}


def count_log_likelihood(counts, probabilities, concentration):
    """Per-cell DM log likelihood with parameter-independent terms omitted."""
    alpha = np.maximum(concentration * probabilities, 1e-12)
    return (gammaln(concentration)-gammaln(concentration+counts.sum(1))
            + (gammaln(alpha+counts)-gammaln(alpha)).sum(1))


def bootstrap_samples(frame, observed, repeats, seed):
    strata = frame.iloc[np.flatnonzero(observed)].groupby(
        ['branch','stage','batch'], sort=True, observed=True).indices.values()
    strata = [np.asarray(ix) for ix in strata]
    rng = np.random.default_rng(seed)
    # Positions refer to the observed-row design. No cells are held out.
    return [np.concatenate([rng.choice(ix,len(ix),replace=True) for ix in strata])
            for _ in range(repeats)]


def quality_design(model,counts,X,frame,times,c):
    from .model import design
    from .composition import settings
    spec=model['design'];observed=counts.sum(1)>0
    D=design(X,frame,times,spec)[observed];width=spec['state_width']
    mask_start=width+sum(f['width'] for f in spec['factors'])+len(spec['interactions'])
    null_columns=list(range(width))+list(range(mask_start,mask_start+len(spec.get('composition_masks',[]))))
    penalty=np.full(D.shape[1],c['model']['ridge'],dtype=float)
    penalty[:width]=c['model']['state_ridge'];penalty[0]=0.
    for column in range(mask_start,mask_start+len(spec.get('composition_masks',[]))):penalty[column]=c['model']['state_ridge']
    for feature in spec.get('compositions',[]):
        penalty[feature['start']:feature['start']+feature['width']]=settings(c)['ridge']
    return D,counts[observed],penalty,null_columns


def _probability(D,fit):
    values=D@np.asarray(fit['coefficients'])
    return softmax(np.column_stack([values,np.zeros(len(D))]),axis=1)


def baseline_fits(D,C,penalty,samples,c):
    from .model import fit_count_design
    return (fit_count_design(D,C,penalty,c),
            [fit_count_design(D[sample],C[sample],penalty,c) for sample in samples])


def relation_quality(model,counts,X,frame,times,c,samples,baseline_cache=None):
    """Paired full/reduced fits scored only on each draw's out-of-bag rows."""
    from .model import fit_count_design
    D,C,penalty,columns=quality_design(model,counts,X,frame,times,c)
    null=D[:,columns];null_penalty=penalty[columns]
    # Mask-free reduced designs are exactly the shared event state basis.
    key=(null.shape,null.tobytes(),null_penalty.tobytes()) if baseline_cache is not None else None
    if key is not None and key in baseline_cache:baseline,reduced_draws=baseline_cache[key]
    else:
        baseline,reduced_draws=baseline_fits(null,C,null_penalty,samples,c)
        if key is not None:baseline_cache[key]=(baseline,reduced_draws)
    empty={'quality_status':'baseline_'+baseline['status'],'fit_gain':None,'oob_fit_gain':None,
           'stability':0.,'bootstrap_successes':0,'oob_positive_gain_fraction':0.,
           'bootstrap_positive_gain_fraction':0.,'bootstrap_effect_agreement':0.,'bootstrap_gain_mean':None}
    if baseline['status']!='ok':return empty
    kappa=c['model']['dm_concentration'];cfg=selection_settings(c)
    p=np.asarray(model['WT_prediction'])[counts.sum(1)>0];p0=_probability(null,baseline)
    gain=float(np.mean(count_log_likelihood(C,p,kappa)-count_log_likelihood(C,p0,kappa)))
    stability=[];agreements=[];gains=[];audit=[];successes=0
    for sample,reduced in zip(samples,reduced_draws):
        oob=np.ones(len(C),bool);oob[np.unique(sample)]=False
        if not oob.any():
            stability.append(0.);audit.append({'status':'empty_oob','oob_cells':0});continue
        full=fit_count_design(D[sample],C[sample],penalty,c)
        if full['status']!='ok' or reduced['status']!='ok':
            stability.append(0.);audit.append({'status':'fit_failed','full':full['status'],'reduced':reduced['status'],'oob_cells':int(oob.sum())});continue
        pb=_probability(D[oob],full);p0b=_probability(null[oob],reduced)
        bgain=float(np.mean(count_log_likelihood(C[oob],pb,kappa)-count_log_likelihood(C[oob],p0b,kappa)))
        effect=(p[oob]-p0[oob]).ravel();beffect=(pb-p0b).ravel()
        denominator=np.linalg.norm(effect)*np.linalg.norm(beffect)
        agreement=float(np.clip(effect@beffect/denominator,0,1)) if denominator>1e-16 else 0.
        successes+=1;gains.append(bgain);agreements.append(agreement)
        stability.append(agreement if bgain>cfg['min_oob_fit_gain'] else 0.)
        audit.append({'status':'ok','oob_cells':int(oob.sum()),'gain':bgain,'agreement':agreement})
    positive=sum(g>cfg['min_oob_fit_gain'] for g in gains)/len(samples)
    return {'quality_status':'ok' if successes else 'no_valid_oob_fits','fit_gain':gain,
            'oob_fit_gain':float(np.mean(gains)) if gains else None,
            'oob_positive_gain_fraction':positive,'stability':float(np.mean(stability)),
            'baseline_mean_log_likelihood':float(count_log_likelihood(C,p0,kappa).mean()),
            'bootstrap_successes':successes,'bootstrap_positive_gain_fraction':positive,
            'bootstrap_effect_agreement':float(np.mean(agreements)) if agreements else 0.,
            'bootstrap_gain_mean':float(np.mean(gains)) if gains else None,'oob_draws':audit}


def rank_relations(db,c):
    """Normalize gain/amplitude within each event; bootstrap stability is [0,1]."""
    db.execute('UPDATE single_relations SET selection_rank=NULL,selected=false,WT_weight=0,selection_score=NULL,normalized_fit_gain=NULL,normalized_response=NULL')
    cfg = selection_settings(c); weights = cfg['weights']; total = sum(weights.values())
    for event, in db.execute('SELECT DISTINCT event_index FROM single_relations').fetchall():
        rows = db.execute('''SELECT factor_id,oob_fit_gain,stability,response_score FROM single_relations
            WHERE event_index=? AND status='ok' AND quality_status='ok'
            AND oob_fit_gain>? AND response_score>? AND stability>=? AND oob_positive_gain_fraction>=?''',
                          [event,cfg['min_oob_fit_gain'],cfg['min_response_score'],cfg['min_stability'],cfg['min_positive_gain_fraction']]).fetchall()
        if not rows:
            continue
        gain_max = max(r[1] for r in rows); response_max = max(r[3] for r in rows)
        values = []
        for fid,gain,stability,response in rows:
            g = gain/gain_max; r = response/response_max
            score = (weights['fit_gain']*g+weights['stability']*stability+weights['response']*r)/total
            values.append((g,r,score,event,fid))
        db.executemany('''UPDATE single_relations SET normalized_fit_gain=?,normalized_response=?,
            selection_score=? WHERE event_index=? AND factor_id=?''',values)
    db.execute('''UPDATE single_relations AS edge SET selection_rank=ranked.r,
        selected=ranked.r<=?, WT_weight=ranked.w FROM (
        SELECT event_index,factor_id,
            row_number() OVER (PARTITION BY event_index ORDER BY selection_score DESC,factor_id) r,
            selection_score/sum(selection_score) OVER (PARTITION BY event_index) w
        FROM single_relations WHERE selection_score>0) ranked
        WHERE edge.event_index=ranked.event_index AND edge.factor_id=ranked.factor_id''',
               [c['model']['max_factors']])


def adjacent_pairs(frame, times):
    """Consecutive cells within context, stage and batch; preserve missing slots."""
    left, right = [], []
    for indices in frame.groupby(['branch', 'stage', 'batch'], sort=True, observed=True).indices.values():
        indices = np.asarray(indices)
        order = np.lexsort((frame.cell.iloc[indices].astype(str).to_numpy(), times[indices]))
        ordered = indices[order]
        a, b = ordered[:-1], ordered[1:]
        valid = np.isfinite(times[a]) & np.isfinite(times[b]) & (times[b] > times[a])
        left.extend(a[valid]); right.extend(b[valid])
    return np.asarray(left, dtype=int), np.asarray(right, dtype=int)


def response_score(model, X, reference, frame, times, observed, pairs, predict_logits):
    """Adapt CausalBridge WT response changes to conditional branch probabilities.

    Hold the source factor at its observed-WT mean and its centered composition
    at zero to remove state-only response. Unsupported adjacent pairs contribute
    zero to a common pair denominator across all candidates for the event.
    """
    j = model['input_indices'][0]
    reference[:, j] = np.mean(X[observed, j])
    for feature in model['design'].get('compositions', []):
        reference[:, feature['index']] = 0.
    baseline = softmax(predict_logits(model, reference, frame, times), axis=1)
    contribution = np.asarray(model['WT_prediction']) - baseline
    left, right = pairs
    supported = observed[left] & observed[right]
    variation = .5 * np.abs(contribution[right[supported]] - contribution[left[supported]]).sum(1)
    score = float(variation.sum() / len(left)) if len(left) else 0.
    # Floating-point differences around a constant fitted response are zero.
    if score <= 32 * np.finfo(float).eps:
        score = 0.
    reference[:, j] = X[:, j]
    for feature in model['design'].get('compositions', []):
        reference[:, feature['index']] = X[:, feature['index']]
    return score, int(supported.sum())


def compact_model(model):
    return {key:model[key] for key in ('status','coefficients','fit_objective','iterations','counted_cells','molecules','message') if key in model}


_WORKER={}


def _init_worker(payload):
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
        os.environ[name]='1'
    from threadpoolctl import threadpool_limits
    from numba import set_num_threads
    threadpool_limits(1);set_num_threads(1)
    _WORKER.update(payload)


def _event_worker(e):
    from .model import fit_event,logits
    return event_relations(e,**_WORKER,fit_single=fit_event,predict_logits=logits)


def event_relations(e,data,frame,times,prior,X,c,fit_single,predict_logits):
    event=data['events'][e];counts=data['C'][:,data['offsets'][e]:data['offsets'][e+1]].toarray()
    observed=counts.sum(1)>0;cfg=selection_settings(c)
    seed=(cfg['bootstrap_seed']+int.from_bytes(hashlib.sha256(event['event'].encode()).digest()[:4],'little'))%(2**32)
    samples=bootstrap_samples(frame,observed,cfg['bootstrap_repeats'],seed)
    pairs=adjacent_pairs(frame,times);reference=X.copy();cache={};records=[]
    findex={f:j for j,f in enumerate(prior['factor_ids'])}
    for candidate in sorted(prior['events'][str(e)],key=lambda r:r['factor_id']):
        fid=candidate['factor_id'];j=findex[fid]
        composition=[f for f in prior.get('composition_features',[]) if f['factor_index']==j and f['gene_id']!=event['gene']]
        model=fit_single(counts,X,frame,times,[j],[],c,composition)
        score,supported=(response_score(model,X,reference,frame,times,observed,pairs,predict_logits)
                         if model['status']=='ok' else (0. if model['status']=='constant_input_on_observed_cells' else None,0))
        quality=(relation_quality(model,counts,X,frame,times,c,samples,cache) if model['status']=='ok' else
                 {'quality_status':'single_fit_failed','fit_gain':None,'oob_fit_gain':None,'oob_positive_gain_fraction':0.,
                  'stability':0.,'bootstrap_successes':0,'bootstrap_positive_gain_fraction':0.,
                  'bootstrap_effect_agreement':0.,'bootstrap_gain_mean':None})
        from .storage import compact_candidate
        compact=compact_candidate(candidate)
        compact.setdefault('anchor_count',len(candidate.get('anchors',[])))
        compact.setdefault('path_count',len(candidate.get('paths',[])))
        records.append((e,event['event'],event['gene'],fid,candidate.get('symbol',fid),model['status'],score,supported,
            json.dumps(compact,allow_nan=False),json.dumps(compact_model(model),allow_nan=False),
            *[quality.get(key) for key in ['quality_status','fit_gain','oob_fit_gain','oob_positive_gain_fraction','stability',
                'bootstrap_successes','bootstrap_positive_gain_fraction','bootstrap_effect_agreement','bootstrap_gain_mean']],
            json.dumps(quality,allow_nan=False)))
    return e,records


def network_fingerprint(data,frame,times,prior,X,c,event_ids):
    from .io import digest
    config={k:v for k,v in c.items() if k not in ('runtime','execution','_config_path','_runtime_threads')}
    sources={p.name:digest(p) for p in Path(__file__).parent.glob('*.py')}
    h=hashlib.sha256()
    for a in (X,times,data['C'].data,data['C'].indices,data['C'].indptr,data['offsets']):
        v=np.ascontiguousarray(a);h.update(str(v.shape).encode());h.update(str(v.dtype).encode());h.update(v.tobytes())
    h.update(frame.to_csv(index=False).encode());h.update(json.dumps(data['events'],sort_keys=True).encode())
    meta={k:v for k,v in prior.items() if k!='events'}
    h.update(json.dumps(meta,sort_keys=True,default=lambda v:v.tolist() if isinstance(v,np.ndarray) else str(v)).encode())
    if hasattr(prior['events'],'path'):h.update(digest(prior['events'].path).encode())
    else:h.update(json.dumps(prior['events'],sort_keys=True).encode())
    return hashlib.sha256(json.dumps({'config':config,'input':h.hexdigest(),'sources':sources,'events':event_ids},sort_keys=True).encode()).hexdigest()


def learn_network(data,frame,times,prior,X,c,fit_single,predict_logits,path=None,event_limit=None):
    """Main process owns the DB; each committed event is an atomic resume unit."""
    if path is not None:Path(path).parent.mkdir(parents=True,exist_ok=True)
    db=duckdb.connect(str(path) if path is not None else ':memory:')
    try:
        event_ids=sorted(map(int,prior['events']))
        if event_limit is not None:event_ids=event_ids[:event_limit]
        fingerprint=network_fingerprint(data,frame,times,prior,X,c,event_ids)
        tables={r[0] for r in db.execute('SHOW TABLES').fetchall()}
        if 'network_metadata' in tables:
            previous=dict(db.execute('SELECT key,value FROM network_metadata').fetchall())
            if previous.get('fingerprint')!=fingerprint:
                raise ValueError('Network resume fingerprint differs (config/input/code/scope); use a new output directory')
        else:
            db.execute("""CREATE TABLE single_relations(event_index INTEGER,event_id VARCHAR,gene_id VARCHAR,
                factor_id VARCHAR,gene_symbol VARCHAR,status VARCHAR,response_score DOUBLE,supported_adjacent_pairs INTEGER,
                candidate_json VARCHAR,model_json VARCHAR,quality_status VARCHAR,fit_gain DOUBLE,oob_fit_gain DOUBLE,
                oob_positive_gain_fraction DOUBLE,stability DOUBLE,bootstrap_successes INTEGER,
                bootstrap_positive_gain_fraction DOUBLE,bootstrap_effect_agreement DOUBLE,bootstrap_gain_mean DOUBLE,
                quality_json VARCHAR,normalized_fit_gain DOUBLE,normalized_response DOUBLE,selection_score DOUBLE,
                selection_rank INTEGER,selected BOOLEAN DEFAULT false,WT_weight DOUBLE DEFAULT 0,
                PRIMARY KEY(event_index,factor_id))""")
            db.execute('CREATE TABLE completed_events(event_index INTEGER PRIMARY KEY,finished_at VARCHAR,candidate_count INTEGER)')
            db.execute('CREATE TABLE network_metadata(key VARCHAR PRIMARY KEY,value VARCHAR)')
            db.executemany('INSERT INTO network_metadata VALUES (?,?)',[
                ('fingerprint',fingerprint),('status','learning_single_relations'),('selection_policy',POLICY),
                ('response_score_definition',SCORE),('quality_settings',json.dumps(selection_settings(c))),
                ('max_factors_after_learning',str(c['model']['max_factors'])),
                ('event_scope','all_prior_supported_events' if event_limit is None else 'limited_event_trial')])
        completed={r[0] for r in db.execute('SELECT event_index FROM completed_events').fetchall()}
        pending=[e for e in event_ids if e not in completed]
        LOG.info('Network resume: %d complete, %d pending events',len(completed),len(pending))
        payload=dict(data=data,frame=frame,times=times,prior=prior,X=X,c=c)
        workers=c.get('_runtime_threads',1)
        executor=None
        if workers>1 and pending:
            # Fork shares read-only count/design arrays; each worker limits BLAS to one thread.
            import multiprocessing as mp
            executor=ProcessPoolExecutor(max_workers=min(workers,len(pending)),mp_context=mp.get_context('fork'),initializer=_init_worker,initargs=(payload,))
            # Bound pending result buffers; a slow early event cannot retain the whole network.
            def results():
                from concurrent.futures import wait,FIRST_COMPLETED
                iterator=iter(pending);active={}
                for e in iterator:
                    active[executor.submit(_event_worker,e)]=e
                    if len(active)>=2*workers:break
                while active:
                    done,_=wait(active,return_when=FIRST_COMPLETED)
                    for future in done:
                        active.pop(future);yield future.result()
                        e=next(iterator,None)
                        if e is not None:active[executor.submit(_event_worker,e)]=e
            result_iterator=results()
        else:
            result_iterator=(event_relations(e,**payload,fit_single=fit_single,predict_logits=predict_logits) for e in pending)
        try:
            for e,records in result_iterator:
                db.execute('BEGIN TRANSACTION')
                try:
                    if records:
                        import pandas as pd
                        columns=['event_index','event_id','gene_id','factor_id','gene_symbol','status','response_score','supported_adjacent_pairs',
                            'candidate_json','model_json','quality_status','fit_gain','oob_fit_gain','oob_positive_gain_fraction','stability',
                            'bootstrap_successes','bootstrap_positive_gain_fraction','bootstrap_effect_agreement','bootstrap_gain_mean','quality_json']
                        batch=pd.DataFrame(records,columns=columns);db.register('event_batch',batch)
                        db.execute('INSERT INTO single_relations ('+','.join(columns)+') SELECT * FROM event_batch');db.unregister('event_batch')
                    db.execute('INSERT INTO completed_events VALUES (?,?,?)',[e,datetime.now(timezone.utc).isoformat(),len(records)])
                    db.execute('COMMIT')
                except BaseException:db.execute('ROLLBACK');raise
                LOG.info('Single-relation network: event %s committed (%d candidates)',data['events'][e]['event'],len(records))
        finally:
            if executor:executor.shutdown(wait=True,cancel_futures=True)
        rank_relations(db,c)
        db.execute("UPDATE network_metadata SET value='selected' WHERE key='status'")
        db.execute('CHECKPOINT')
        return db
    except BaseException:db.close();raise


def network_summary(db):
    counts = db.execute('''SELECT count(*),count(*) FILTER (WHERE status='ok'),
        count(*) FILTER (WHERE response_score>0), count(*) FILTER (WHERE selected),
        count(DISTINCT event_index),count(DISTINCT factor_id) FROM single_relations''').fetchone()
    return dict(zip(['candidate_relations', 'fitted_single_relations', 'positive_score_relations',
                     'selected_relations', 'events', 'factors'], counts),
                selection_policy=POLICY, response_score_definition=SCORE,
                quality_eligible_relations=db.execute('SELECT count(*) FROM single_relations WHERE selection_rank IS NOT NULL').fetchone()[0],
                quality_settings=json.loads(db.execute("SELECT value FROM network_metadata WHERE key='quality_settings'").fetchone()[0]),
                fit_status_counts=dict(db.execute(
                    'SELECT status,count(*) FROM single_relations GROUP BY status').fetchall()))
