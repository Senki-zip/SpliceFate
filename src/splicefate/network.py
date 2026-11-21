"""Learn all prior-supported WT single relations before selecting joint inputs."""
from __future__ import annotations

import json
import logging
from pathlib import Path

import duckdb
import numpy as np
from scipy.special import softmax


LOG = logging.getLogger('splicefate')
POLICY = 'all_single_relations_then_top_response'
SCORE = 'mean_adjacent_WT_factor_response_total_variation'


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
            selection_rank INTEGER, selected BOOLEAN DEFAULT false,
            WT_weight DOUBLE DEFAULT 0,
            PRIMARY KEY(event_index, factor_id))''')
        db.execute('CREATE OR REPLACE TABLE network_metadata(key VARCHAR, value VARCHAR)')
        db.executemany('INSERT INTO network_metadata VALUES (?, ?)', [
            ('status', 'learning_single_relations'), ('selection_policy', POLICY),
            ('response_score_definition', SCORE),
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
            records = []
            for candidate in sorted(prior['events'][str(e)], key=lambda r: r['factor_id']):
                fid = candidate['factor_id']; j = findex[fid]
                composition = [f for f in prior.get('composition_features', [])
                               if f['factor_index'] == j and f['gene_id'] != event['gene']]
                model = fit_single(counts, X, frame, times, [j], [], c, composition)
                score, supported = (response_score(model, X, reference, frame, times,
                                                  observed, pairs, predict_logits)
                                    if model['status'] == 'ok' else (None, 0))
                records.append((e, event['event'], event['gene'], fid,
                                candidate.get('symbol', fid), model['status'], score, supported,
                                json.dumps(candidate, allow_nan=False),
                                json.dumps(compact_model(model), allow_nan=False)))
            if records:
                db.executemany('''INSERT INTO single_relations(event_index,event_id,gene_id,
                    factor_id,gene_symbol,status,response_score,supported_adjacent_pairs,
                    candidate_json,model_json) VALUES (?,?,?,?,?,?,?,?,?,?)''', records)
            LOG.info('Single-relation network: %d/%d events, %d candidates fitted for %s',
                     number + 1, len(event_ids), len(records), event['event'])
        db.execute("UPDATE network_metadata SET value='single_relations_learned' WHERE key='status'")
        minimum = float(c['model']['network_selection'].get('min_response_score', 0.))
        db.execute('''UPDATE single_relations AS edge SET selection_rank=ranked.r,
            selected=ranked.r<=?, WT_weight=ranked.w
            FROM (SELECT event_index,factor_id,
                row_number() OVER (PARTITION BY event_index ORDER BY response_score DESC,factor_id) AS r,
                response_score/sum(response_score) OVER (PARTITION BY event_index) AS w
                FROM single_relations WHERE status='ok' AND response_score>?) ranked
            WHERE edge.event_index=ranked.event_index AND edge.factor_id=ranked.factor_id''',
                   [c['model']['max_factors'], minimum])
        db.execute("UPDATE network_metadata SET value='selected' WHERE key='status'")
        return db
    except BaseException:
        db.close()
        raise


def network_summary(db):
    counts = db.execute('''SELECT count(*),count(*) FILTER (WHERE status='ok'),
        count(*) FILTER (WHERE selection_rank IS NOT NULL), count(*) FILTER (WHERE selected),
        count(DISTINCT event_index),count(DISTINCT factor_id) FROM single_relations''').fetchone()
    return dict(zip(['candidate_relations', 'fitted_single_relations', 'positive_score_relations',
                     'selected_relations', 'events', 'factors'], counts),
                selection_policy=POLICY, response_score_definition=SCORE,
                fit_status_counts=dict(db.execute(
                    'SELECT status,count(*) FROM single_relations GROUP BY status').fetchall()))
