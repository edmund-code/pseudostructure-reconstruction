"""Cross-species conservation of gene gradients by gradient strength, and resampling of the index (notebook 49).

Notebook 43's conservation index divides the noise-corrected human–mouse correlation of gene
gradients by same-species, cross-dataset ceilings (`zonation_reliability`). This module asks
*where* along the strength axis that agreement lives, and how uncertain the index is.

Strength bins must not be defined from the data being correlated. Binning on one dataset's
|gradient| selects its noise extremes (regression to the mean) and inflates that dataset's
apparent strength in the top bins. Here strength comes from one donor fold of the external
references and every correlation uses the other fold (cross-fitting). Each bin gets its own
ceilings, so a near-zero cross-species correlation is read against the agreement that two
datasets of one species reach for the same genes.

The estimator matches notebook 43 when ``match_reliability=False`` and ``min_reliability=0``.
With ``match_reliability=True``, reliability is computed on the same winsorised values as the
correlation. The sampling variance of a clipped gene is kept as the variance of its unclipped mean;
that overstates its noise (clipping is a contraction), so matched reliabilities are slightly low.
Only a fraction 2q of genes is affected. The upgrade path is to winsorise replicates before
averaging.
"""
from __future__ import annotations

import warnings
from hashlib import sha256

import numpy as np
import pandas as pd

STATS = ('genes', 'r_cross', 'rstar_cross', 'r_ceiling_human', 'rstar_ceiling_human', 'r_ceiling_mouse',
         'rstar_ceiling_mouse', 'index', 'r_external', 'rstar_external', 'index_external',
         'reliability_human', 'reliability_mouse', 'reliability_human_ref', 'reliability_mouse_ref',
         'sign_agreement')
DATASETS = ('human', 'mouse', 'human_ref', 'mouse_ref')
DEFAULT_EDGES = (0., .5, .8, .9, .95, .99, 1.)


def winsorize_array(values, q):
    """Clip at the q and 1 - q quantiles (linear interpolation, as pandas); q = 0 returns the input."""
    values = np.asarray(values, float)
    if q <= 0:
        return values
    low, high = np.quantile(values, [q, 1 - q])
    return np.clip(values, low, high)


def reliability_array(mean, noise):
    """1 - mean(noise) / var(mean) across genes, floored at 0 (notebook 43's definition)."""
    total = float(np.var(mean, ddof=1)) if len(mean) > 1 else 0.
    return max(1 - float(np.mean(noise)) / total, 0.) if total > 0 else np.nan


def _correlation(x, y):
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def subset_statistics(data, idx, *, q=.01, match_reliability=True, min_reliability=0.):
    """Raw and noise-corrected correlations, ceilings and both indices on the genes ``idx``.

    ``data`` maps 'human', 'mouse' (the cross-species pair), 'human_ref' and 'mouse_ref' (the
    same-species references) to (mean, noise) arrays aligned on one gene list. ``idx`` holds gene
    positions and may repeat (bootstrap). Winsorisation is applied within ``idx``. A corrected
    correlation is NaN unless both reliabilities are > 0 and >= ``min_reliability``; an index is
    NaN unless both ceilings are > 0. Returns a dict keyed by :data:`STATS`.
    """
    idx = np.asarray(idx, int)
    if len(idx) < 3:
        return {k: (len(idx) if k == 'genes' else np.nan) for k in STATS}
    values, rel = {}, {}
    for key in DATASETS:
        mean, noise = data[key]
        raw = np.asarray(mean, float)[idx]
        values[key] = winsorize_array(raw, q)
        rel[key] = reliability_array(values[key] if match_reliability else raw, np.asarray(noise, float)[idx])

    def corrected(a, b):
        r = _correlation(values[a], values[b])
        ok = all(np.isfinite(rel[k]) and rel[k] > 0 and rel[k] >= min_reliability for k in (a, b))
        return r, (r / np.sqrt(rel[a] * rel[b]) if ok else np.nan)

    def index(cross, ceiling_h, ceiling_m):
        return cross / np.sqrt(ceiling_h * ceiling_m) if ceiling_h > 0 and ceiling_m > 0 else np.nan

    r_cross, rs_cross = corrected('human', 'mouse')
    r_ch, rs_ch = corrected('human', 'human_ref')
    r_cm, rs_cm = corrected('mouse', 'mouse_ref')
    r_ext, rs_ext = corrected('human_ref', 'mouse_ref')
    h, m = np.asarray(data['human'][0], float)[idx], np.asarray(data['mouse'][0], float)[idx]
    return {'genes': len(idx), 'r_cross': r_cross, 'rstar_cross': rs_cross, 'r_ceiling_human': r_ch,
            'rstar_ceiling_human': rs_ch, 'r_ceiling_mouse': r_cm, 'rstar_ceiling_mouse': rs_cm,
            'index': index(rs_cross, rs_ch, rs_cm), 'r_external': r_ext, 'rstar_external': rs_ext,
            'index_external': index(rs_ext, rs_ch, rs_cm), 'reliability_human': rel['human'],
            'reliability_mouse': rel['mouse'], 'reliability_human_ref': rel['human_ref'],
            'reliability_mouse_ref': rel['mouse_ref'], 'sign_agreement': float(np.mean(np.sign(h) == np.sign(m)))}


def percentile_rank(values):
    """Percentile rank of |values| in (0, 1], ties averaged."""
    return pd.Series(np.abs(np.asarray(values, float))).rank(method='average', pct=True).to_numpy()


def strength_score(ranks, how='max'):
    """Combine per-species percentile ranks (list of arrays) into one strength score."""
    stacked = np.vstack([np.asarray(r, float) for r in ranks])
    if how == 'max':
        return stacked.max(axis=0)
    if how == 'mean':
        return stacked.mean(axis=0)
    raise ValueError(f'Unknown combination: {how}')


def strength_bins(score, edges=DEFAULT_EDGES):
    """Integer bin per gene from quantiles of ``score`` (ties broken by order, so bin sizes are exact)."""
    ranks = pd.Series(np.asarray(score, float)).rank(method='first')
    return pd.qcut(ranks, list(edges), labels=False).to_numpy(int)


def bin_labels(edges=DEFAULT_EDGES):
    return [f'{100 * a:g}–{100 * b:g}%' for a, b in zip(edges[:-1], edges[1:])]


def donor_folds(donors, n_folds=2):
    """Deterministic, balanced donor folds: sort by sha256(donor id), then deal round-robin."""
    donors = [str(d) for d in donors]
    order = sorted(donors, key=lambda d: sha256(d.encode()).hexdigest())
    fold = {d: i % n_folds for i, d in enumerate(order)}
    return np.array([fold[d] for d in donors], dtype=int)


def binned_statistics(directions, idx, subsets, **kwargs):
    """Per-subset statistics averaged over cross-fitting directions.

    ``directions``: list of (data, bins) pairs, one per direction; ``bins`` gives each gene's bin
    in that direction. ``subsets``: dict name -> list of bin ids. ``idx``: gene positions (a
    bootstrap draw, or all genes). Returns an array (subsets x STATS). A statistic that is NaN in
    any direction is NaN in the average.
    """
    idx = np.asarray(idx, int)
    out = np.zeros((len(subsets), len(STATS)))
    for data, bins in directions:
        drawn_bins = np.asarray(bins)[idx]
        for i, members in enumerate(subsets.values()):
            stats_ = subset_statistics(data, idx[np.isin(drawn_bins, members)], **kwargs)
            out[i] += [stats_[k] for k in STATS]
    return out / len(directions)


def gene_draws(n_genes, n_boot, seed):
    """Bootstrap draws of gene positions; the same sequence as notebook 43's gene bootstrap."""
    rng = np.random.default_rng(seed)
    for _ in range(n_boot):
        yield rng.integers(0, n_genes, n_genes)


def block_draws(blocks, n_boot, seed):
    """Block bootstrap: resample blocks with replacement and take all genes of each drawn block."""
    blocks = np.asarray(blocks)
    labels = np.unique(blocks)
    members = [np.flatnonzero(blocks == b) for b in labels]
    rng = np.random.default_rng(seed)
    for _ in range(n_boot):
        yield np.concatenate([members[c] for c in rng.integers(0, len(labels), len(labels))])


def percentile_interval(draws, *, level=.95, max_missing=.1):
    """Percentile interval over non-missing draws (axis 0); NaN where more than ``max_missing`` are missing."""
    draws = np.asarray(draws, float)
    missing = np.isnan(draws).mean(axis=0)
    tail = (1 - level) / 2 * 100
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        low, high = np.nanpercentile(draws, [tail, 100 - tail], axis=0)
    low, high = np.where(missing > max_missing, np.nan, low), np.where(missing > max_missing, np.nan, high)
    return low, high, missing


def resampled_gradient(per_donor, rng):
    """Mean gradient and sampling variance (var / n) after resampling donors (rows) with replacement."""
    per_donor = np.asarray(per_donor, float)
    rows = rng.integers(0, len(per_donor), len(per_donor))
    sample = per_donor[rows]
    return sample.mean(axis=0), sample.var(axis=0, ddof=1) / len(sample)


def coexpression_modules(matrix, specimen, species, *, n_components=50, n_modules=100, seed=0):
    """Gene modules from co-variation across structures (for a module-block bootstrap).

    ``matrix``: structures x genes (dense or sparse) log expression. Each gene is centred within
    specimen and scaled to unit SD within species; each species block is weighted by
    1/sqrt(structures) so both species count equally. Genes are embedded by a truncated SVD
    (loadings scaled by singular values, L2-normalised) and clustered by k-means.
    """
    from scipy import sparse
    from sklearn.cluster import KMeans
    from sklearn.utils.extmath import randomized_svd

    X = matrix.toarray() if sparse.issparse(matrix) else np.array(matrix, dtype=float)
    X = X.astype(np.float32)
    specimen, species = np.asarray(specimen), np.asarray(species)
    for name in np.unique(specimen):
        rows = specimen == name
        X[rows] -= X[rows].mean(axis=0)
    for name in np.unique(species):
        rows = species == name
        sd = X[rows].std(axis=0)
        X[rows] /= np.where(sd > 0, sd, 1)[None, :] * np.sqrt(rows.sum())
    _, singular, vt = randomized_svd(X, n_components=n_components, random_state=seed)
    loadings = vt.T * singular[None, :]
    loadings /= np.clip(np.linalg.norm(loadings, axis=1, keepdims=True), 1e-12, None)
    return KMeans(n_clusters=n_modules, n_init=10, random_state=seed).fit_predict(loadings)
