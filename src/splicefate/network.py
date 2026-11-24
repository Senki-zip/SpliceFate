"""Learn all prior-supported WT single relations before selecting joint inputs."""
from __future__ import annotations

import json
import hashlib
import logging
from pathlib import Path

import duckdb
import numpy as np
from scipy.special import softmax, gammaln


LOG = logging.getLogger('splicefate')
POLICY = 'all_single_relations_then_weighted_quality'
SCORE = 'mean_adjacent_WT_factor_response_total_variation'


def selection_settings(c):
    return {'method':POLICY, 'min_response_score':0., 'min_fit_gain':0.,
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


def relation_quality(model, counts, X, frame, times, c, samples):
    """Matched nested gain and fixed-WT-basis bootstrap response reproducibility."""
    from .model import design, fit_count_design
    spec = model['design']; observed = counts.sum(1)>0
    D = design(X,frame,times,spec)[observed]; C = counts[observed]
    width = spec['state_width']
    mask_start = width+sum(f['width'] for f in spec['factors'])+len(spec['interactions'])
    null_columns = list(range(width))+list(range(mask_start,mask_start+len(spec.get('composition_masks',[]))))
    null = D[:,null_columns]
    penalty = np.full(D.shape[1],c['model']['ridge'],dtype=float)
    penalty[:width] = c['model']['state_ridge']; penalty[0] = 0.
    for column in range(mask_start,mask_start+len(spec.get('composition_masks',[]))):
        penalty[column] = c['model']['state_ridge']
    from .composition import settings
    for feature in spec.get('compositions',[]):
        penalty[feature['start']:feature['start']+feature['width']] = settings(c)['ridge']
    null_penalty = penalty[null_columns]
    def predict(matrix, fit):
        values = matrix@np.asarray(fit['coefficients'])
        return softmax(np.column_stack([values,np.zeros(len(matrix))]),axis=1)
    baseline = fit_count_design(null,C,null_penalty,c)
    if baseline['status']!='ok':
        return {'quality_status':'baseline_'+baseline['status'], 'fit_gain':None,
                'stability':0., 'bootstrap_successes':0, 'bootstrap_positive_gain_fraction':0.,
                'bootstrap_effect_agreement':0., 'bootstrap_gain_mean':None}
    kappa = c['model']['dm_concentration']
    p = np.asarray(model['WT_prediction'])[observed]; p0 = predict(null,baseline)
    gain = float(np.mean(count_log_likelihood(C,p,kappa)-count_log_likelihood(C,p0,kappa)))
    effect = (p-p0).ravel(); norm = np.linalg.norm(effect)
    stability = []; agreements = []; gains = []; successes = 0
    for sample in samples:
        full = fit_count_design(D[sample],C[sample],penalty,c)
        reduced = fit_count_design(null[sample],C[sample],null_penalty,c)
        if full['status']!='ok' or reduced['status']!='ok':
            stability.append(0.); continue
        successes += 1
        pb = predict(D,full); p0b = predict(null,reduced)
        bgain = float(np.mean(count_log_likelihood(C,pb,kappa)-count_log_likelihood(C,p0b,kappa)))
        beffect = (pb-p0b).ravel(); denominator = norm*np.linalg.norm(beffect)
        agreement = float(np.clip(effect@beffect/denominator,0,1)) if denominator>1e-16 else 0.
        gains.append(bgain); agreements.append(agreement)
        stability.append(agreement if bgain>selection_settings(c)['min_fit_gain'] else 0.)
    return {'quality_status':'ok', 'fit_gain':gain, 'baseline_mean_log_likelihood':float(count_log_likelihood(C,p0,kappa).mean()),
            'stability':float(np.mean(stability)), 'bootstrap_successes':successes,
            'bootstrap_positive_gain_fraction':sum(g>selection_settings(c)['min_fit_gain'] for g in gains)/len(samples),
            'bootstrap_effect_agreement':float(np.mean(agreements)) if agreements else 0.,
            'bootstrap_gain_mean':float(np.mean(gains)) if gains else None}


def rank_relations(db,c):
    """Normalize gain/amplitude within each event; bootstrap stability is [0,1]."""
    cfg = selection_settings(c); weights = cfg['weights']; total = sum(weights.values())
    for event, in db.execute('SELECT DISTINCT event_index FROM single_relations').fetchall():
        rows = db.execute('''SELECT factor_id,fit_gain,stability,response_score FROM single_relations
            WHERE event_index=? AND status='ok' AND quality_status='ok'
            AND fit_gain>? AND response_score>? AND stability>0''',
                          [event,cfg['min_fit_gain'],cfg['min_response_score']]).fetchall()
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
    """Store fitted curves; cell-level arrays remain in the prepared WT inputs."""
    return {key: value for key, value in model.items()
            if key not in ('WT_prediction', 'support', 'observed_denominators')}


def learn_network(data, frame, times, prior, X, c, fit_single, predict_logits, path=None,
                  event_limit=None):
    """Complete the single-relation pass, then rank and select at most K inputs."""
    if path is not None:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    db = duckdb.connect(str(path) if path is not None else ':memory:')
    try:
        db.execute('''CREATE OR REPLACE TABLE single_relations(
            event_index INTEGER, event_id VARCHAR, gene_id VARCHAR,
            factor_id VARCHAR, gene_symbol VARCHAR, status VARCHAR,
            response_score DOUBLE, supported_adjacent_pairs INTEGER,
            candidate_json VARCHAR, model_json VARCHAR,
            quality_status VARCHAR, fit_gain DOUBLE, stability DOUBLE,
            bootstrap_successes INTEGER, bootstrap_positive_gain_fraction DOUBLE,
            bootstrap_effect_agreement DOUBLE, bootstrap_gain_mean DOUBLE,
            quality_json VARCHAR, normalized_fit_gain DOUBLE, normalized_response DOUBLE,
            selection_score DOUBLE,
            selection_rank INTEGER, selected BOOLEAN DEFAULT false,
            WT_weight DOUBLE DEFAULT 0,
            PRIMARY KEY(event_index, factor_id))''')
        db.execute('CREATE OR REPLACE TABLE network_metadata(key VARCHAR, value VARCHAR)')
        db.executemany('INSERT INTO network_metadata VALUES (?, ?)', [
            ('status', 'learning_single_relations'), ('selection_policy', POLICY),
            ('response_score_definition', SCORE),
            ('quality_settings',json.dumps(selection_settings(c))),
            ('reference_implementation', 'CausalBridge_en:saturating_hybrid.fit_target/_wt_scores'),
            ('max_factors_after_learning', str(c['model']['max_factors'])),
            ('event_scope', 'all_prior_supported_events' if event_limit is None else 'limited_event_trial'),
        ])
        findex = {f: j for j, f in enumerate(prior['factor_ids'])}
        pairs = adjacent_pairs(frame, times)
        reference = X.copy()
        event_ids = sorted(map(int, prior['events']))
        if event_limit is not None:
            event_ids = event_ids[:event_limit]
        for number, e in enumerate(event_ids):
            event = data['events'][e]
            counts = data['C'][:, data['offsets'][e]:data['offsets'][e + 1]].toarray()
            observed = counts.sum(1) > 0
            cfg = selection_settings(c)
            event_seed = (cfg['bootstrap_seed']+int.from_bytes(
                hashlib.sha256(event['event'].encode()).digest()[:4],'little')) % (2**32)
            samples = bootstrap_samples(frame,observed,cfg['bootstrap_repeats'],event_seed)
            records = []
            for candidate in sorted(prior['events'][str(e)], key=lambda r: r['factor_id']):
                fid = candidate['factor_id']; j = findex[fid]
                composition = [f for f in prior.get('composition_features', [])
                               if f['factor_index'] == j and f['gene_id'] != event['gene']]
                model = fit_single(counts, X, frame, times, [j], [], c, composition)
                score, supported = (response_score(model, X, reference, frame, times,
                                                  observed, pairs, predict_logits)
                                    if model['status'] == 'ok' else (None, 0))
                quality = (relation_quality(model,counts,X,frame,times,c,samples)
                           if model['status']=='ok' else {'quality_status':'single_fit_failed',
                               'fit_gain':None,'stability':0.,'bootstrap_successes':0,
                               'bootstrap_positive_gain_fraction':0.,'bootstrap_effect_agreement':0.,
                               'bootstrap_gain_mean':None})
                records.append((e, event['event'], event['gene'], fid,
                                candidate.get('symbol', fid), model['status'], score, supported,
                                json.dumps(candidate, allow_nan=False),
                                json.dumps(compact_model(model), allow_nan=False),
                                *[quality[key] for key in ['quality_status','fit_gain','stability',
                                    'bootstrap_successes','bootstrap_positive_gain_fraction',
                                    'bootstrap_effect_agreement','bootstrap_gain_mean']],
                                json.dumps(quality,allow_nan=False)))
            if records:
                db.executemany('''INSERT INTO single_relations(event_index,event_id,gene_id,
                    factor_id,gene_symbol,status,response_score,supported_adjacent_pairs,
                    candidate_json,model_json,quality_status,fit_gain,stability,bootstrap_successes,
                    bootstrap_positive_gain_fraction,bootstrap_effect_agreement,bootstrap_gain_mean,
                    quality_json) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''', records)
            LOG.info('Single-relation network: %d/%d events, %d candidates fitted for %s',
                     number + 1, len(event_ids), len(records), event['event'])
        db.execute("UPDATE network_metadata SET value='single_relations_learned' WHERE key='status'")
        rank_relations(db,c)
        db.execute("UPDATE network_metadata SET value='selected' WHERE key='status'")
        return db
    except BaseException:
        db.close()
        raise


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
