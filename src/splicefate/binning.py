"""Frozen, partition-local Gaussian kernels and cell-level bootstrap views."""
from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple

import numpy as np
import pandas as pd


DEFAULTS = dict(enabled=True, mode='fixed_gaussian', n_bins=50, min_cells_per_bin=10,
                overlap_factor=2., min_support_cells=5, min_oob_bins=4)
PARTITION = ['branch', 'stage', 'batch']
BIN_COLUMNS = ['bin_id', *PARTITION, 'pseudotime', 'ess', 'support_cells', 'valid']
ARRAY_NAMES = ('K_raw', 'W_agg', 'K_eff', 'M_cell', 'centers', 'sigmas', 'valid', 'ess', 'support_cells')


def settings(c):
    return {**DEFAULTS, **c.get('binning', {})}


@dataclass(frozen=True)
class BinPlan:
    K_raw: np.ndarray
    W_agg: np.ndarray
    K_eff: np.ndarray
    M_cell: np.ndarray
    centers: np.ndarray
    sigmas: np.ndarray
    valid: np.ndarray
    ess: np.ndarray
    support_cells: np.ndarray
    bin_frame: pd.DataFrame

    def __post_init__(self):
        for name in ARRAY_NAMES:
            getattr(self, name).setflags(write=False)

    def arrays(self):
        return {'bin_' + name: getattr(self, name) for name in ARRAY_NAMES}

    def summary(self):
        support = CellSupport.from_rows(self.W_agg[self.valid])
        return {'enabled': True, 'mode': 'fixed_gaussian', 'bins': len(self.valid),
                'valid_bins': int(self.valid.sum()), 'cells': self.M_cell.shape[0],
                'covered_cells': int(np.count_nonzero(self.M_cell.any(axis=1))), **vars(support)}

    def diagnostics(self):
        return {**self.summary(), 'bin_frame': self.bin_frame.to_dict('records')}

    @classmethod
    def load(cls, arrays, diagnostics):
        frame = pd.DataFrame(diagnostics['bin_frame'], columns=BIN_COLUMNS)
        return cls(*(arrays['bin_' + name] for name in ARRAY_NAMES), frame)


@dataclass(frozen=True)
class CellSupport:
    pooled_unique_cell_ess: float
    kernel_effective_rows: float

    @classmethod
    def from_rows(cls, weights):
        return cls(cell_ess(weights), effective_rows(weights))


class BootstrapBinView(NamedTuple):
    W_tr: np.ndarray
    K_tr: np.ndarray
    valid_tr: np.ndarray
    W_oo: np.ndarray
    K_oo: np.ndarray
    valid_oo: np.ndarray
    oob_cells: np.ndarray


def cell_ess(weights):
    if not len(weights):
        return 0.
    pooled = weights.mean(axis=0)
    mass = np.dot(pooled, pooled)
    return float(1. / mass) if mass else 0.


def effective_rows(weights):
    gram = weights @ weights.T
    mass = np.sum(gram ** 2)
    return float(np.trace(gram) ** 2 / mass) if mass else 0.


def fit_bin_plan(frame, times, c):
    cfg = settings(c)
    if not cfg['enabled']:
        return None
    times = np.asarray(times, dtype=float)
    if times.shape != (len(frame),) or not np.isfinite(times).all() or frame[PARTITION].isna().any().any():
        raise ValueError('Binning requires finite pseudotime and complete cell partitions')

    # Partition-local requests keep a common cells-per-bin scale. In aggregate,
    # the overlap budget grows with cell count plus the per-partition minimum.
    blocks = []
    strata = frame.groupby(PARTITION, sort=True, observed=True).indices
    for labels, ix in strata.items():
        local = times[ix]
        n = len(local)
        low, high = local.min(), local.max()
        if low == high:
            k = 1
        elif n < cfg['min_cells_per_bin']:
            k = min(cfg['n_bins'], max(1, n // cfg['min_support_cells']))
        else:
            k = min(cfg['n_bins'], max(5, int(cfg['overlap_factor'] * n / cfg['min_cells_per_bin'])))
        width = max(high - low, 1e-8) / k
        centers = np.linspace(low + width / 2, high - width / 2, k) if high > low else np.array([low])
        sigma = width * max(1., cfg['min_cells_per_bin'] / max(n / k, 1.))
        kernel = np.exp(-.5 * ((centers[:, None] - local) / sigma) ** 2)
        blocks.append((labels, ix, centers, sigma, kernel))

    raw = np.zeros((sum(len(block[2]) for block in blocks), len(frame)))
    records = []
    for labels, ix, centers, sigma, kernel in blocks:
        start = len(records)
        raw[start:start + len(centers), ix] = kernel
        supports = np.count_nonzero(kernel > 0, axis=1)
        sums = kernel.sum(1)
        sum_squares = np.square(kernel).sum(1)
        ess = np.divide(sums ** 2, sum_squares, out=np.zeros_like(sums), where=sum_squares > 0)
        for j, center in enumerate(centers):
            valid = sums[j] > 0 and supports[j] >= cfg['min_support_cells']
            records.append((start + j, *labels, center, ess[j], supports[j], valid, sigma))

    bins = pd.DataFrame(records, columns=[*BIN_COLUMNS, 'sigma'])
    sigmas = bins.pop('sigma').to_numpy(dtype=float)
    valid = bins.valid.to_numpy(dtype=bool)
    W = raw / (raw.sum(1, keepdims=True) + 1e-10)
    mass = raw[valid].sum(0)
    effective = np.divide(raw, mass, out=np.zeros_like(raw), where=valid[:, None] & (mass > 0))
    return BinPlan(raw, W, effective, effective.T.copy(), bins.pseudotime.to_numpy(dtype=float),
                   sigmas, valid, bins.ess.to_numpy(dtype=float), bins.support_cells.to_numpy(dtype=int), bins)


def cell_bootstrap_bin_views(plan, frame, repeats, seed, min_oob_bins=4):
    """Project disjoint stratified draws onto the unchanged WT kernels."""
    rng = np.random.default_rng(seed)
    strata = list(frame.groupby(PARTITION, sort=True, observed=True).indices.values())

    def project(multiplicity):
        K = plan.K_raw * multiplicity
        mass = K.sum(1, keepdims=True)
        K /= mass + 1e-10
        return K, plan.K_eff * multiplicity, plan.valid & (mass[:, 0] > 0)

    views = []
    for _ in range(repeats):
        sample = np.concatenate([rng.choice(ix, len(ix), replace=True) for ix in strata]) if strata else []
        train = np.bincount(np.asarray(sample, dtype=int), minlength=len(frame))
        oob_cells = train == 0
        W_tr, K_tr, valid_tr = project(train)
        W_oo, K_oo, valid_oo = project(oob_cells)
        if valid_oo.sum() < min_oob_bins:
            # Preserve draw alignment and cell identities for audit and refit.
            valid_oo[:] = False
        views.append(BootstrapBinView(W_tr, K_tr, valid_tr, W_oo, K_oo, valid_oo, oob_cells))
    return views
