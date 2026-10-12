"""AKI along the PT on the package pseudospace (notebook 76): segment units, program profiles and their summaries,
empirical calibration against relabelings, and the placement of units that have no PT position.

Protocol: ``docs/results/pt-package-aki.md``. Everything here is descriptive arithmetic on fitted curves; the gene
model itself is ``pathway_pipelines.standard_gene_model``.

- ``segment_cuts`` / ``segment_units``: the S1->S2 and S2->S3 cuts that maximise balanced accuracy on reference
  units, and the piecewise-linear map of a position to segment units (S1 0-1, S2 1-2, S3 2-3).
- ``centre_rows`` / ``program_profile``: Delta(x) minus its median over reference genes at each x (notebook 68's
  centring), then the median over a program's members, in log2.
- ``profile_summary``: level (grid mean), drift (last minus first ``window`` points), centre of mass of the positive
  part in segment units, peak and the change across each segment's stretch.
- ``empirical_p``: a p-value calibrated against a pooled null sample (the relabelings' p).
- ``equal_weight_quantiles``: quantiles that give each group (specimen) the same total weight.
- ``knn_place``: the median value and label shares of each query's k nearest reference units.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

LN2 = np.log(2.)


def balanced_cut(position, late):
    """The cut c maximising balanced accuracy of "late iff position > c" (midpoints of sorted distinct positions)."""
    x = np.asarray(position, float)
    y = np.asarray(late, bool)
    if not (y.any() and (~y).any()):
        raise ValueError('Need units on both sides of the cut.')
    order = np.argsort(x, kind='stable')
    xs, ys = x[order], y[order]
    early_below = np.cumsum(~ys)                  # early units at or below each sorted position
    late_above = y.sum() - np.cumsum(ys)          # late units above each sorted position
    score = early_below / (~y).sum() + late_above / y.sum()
    distinct = np.r_[xs[1:] != xs[:-1], False]    # a cut can only fall between two distinct positions
    if not distinct.any():
        raise ValueError('All positions are equal.')
    score = np.where(distinct, score, -np.inf)
    i = int(np.argmax(score))
    return float((xs[i] + xs[i + 1]) / 2)


def segment_cuts(position, segment, labels=('PT-S1', 'PT-S2', 'PT-S3')):
    """(t12, t23): balanced-accuracy cuts between adjacent segment labels; must increase."""
    x = np.asarray(position, float)
    seg = np.asarray(segment).astype(str)
    s1, s2, s3 = labels
    a, b = np.isin(seg, [s1, s2]), np.isin(seg, [s2, s3])
    t12 = balanced_cut(x[a], seg[a] == s2)
    t23 = balanced_cut(x[b], seg[b] == s3)
    if not t12 < t23:
        raise ValueError(f'Segment cuts are not ordered: {t12} >= {t23}')
    return t12, t23


def segment_units(position, t12, t23, lo=0., hi=1.):
    """Piecewise-linear map [lo, t12] -> [0, 1], [t12, t23] -> [1, 2], [t23, hi] -> [2, 3], extended linearly."""
    x = np.asarray(position, float)
    if not lo < t12 < t23 < hi:
        raise ValueError('Need lo < t12 < t23 < hi.')
    return np.where(x < t12, (x - lo) / (t12 - lo),
                    np.where(x < t23, 1 + (x - t12) / (t23 - t12), 2 + (x - t23) / (hi - t23)))


def centre_rows(delta, reference):
    """Delta (genes x grid) minus, at each grid point, its median over the ``reference`` genes; also returns m(x)."""
    d = np.asarray(delta, float)
    ref = np.asarray(reference, bool)
    if d.ndim != 2 or ref.shape != (d.shape[0],) or not ref.any():
        raise ValueError('Need a genes x grid matrix and a nonempty reference mask over its rows.')
    m = np.nanmedian(d[ref], axis=0)
    return d - m, m


def program_profile(centred, genes, members, usable=None):
    """Median over the program's usable members of the centred natural-log Delta, in log2; (profile, n members)."""
    index = pd.Index(genes)
    rows = index.get_indexer(pd.Index(members).unique())
    rows = rows[rows >= 0]
    if usable is not None:
        rows = rows[np.asarray(usable, bool)[rows]]
    if len(rows) == 0:
        return np.full(np.asarray(centred).shape[1], np.nan), 0
    return np.nanmedian(np.asarray(centred, float)[rows], axis=0) / LN2, int(len(rows))


def centre_of_mass(profile, s):
    """Sum s * P+ / sum P+ (P+ = max(P, 0)); NaN when the profile has no positive part."""
    p = np.clip(np.asarray(profile, float), 0, None)
    total = p.sum()
    return float((np.asarray(s, float) * p).sum() / total) if total > 0 else np.nan


def profile_summary(profile, s, window=21):
    """Level, drift, centre of mass, peak (segment units) and the change across each segment's stretch of the grid."""
    p, s = np.asarray(profile, float), np.asarray(s, float)
    if p.shape != s.shape or len(p) < 2 * window:
        raise ValueError('Need a profile and segment units of the same length, at least two windows long.')
    out = {'level': float(p.mean()), 'drift': float(p[-window:].mean() - p[:window].mean()),
           'centre_of_mass': centre_of_mass(p, s), 'peak': float(s[int(np.argmax(p))]),
           'peak_value': float(p.max())}
    for k, name in enumerate(('S1', 'S2', 'S3')):
        inside = np.flatnonzero((s >= k) & (s <= k + 1))
        out[f'within_{name}'] = float(p[inside[-1]] - p[inside[0]]) if len(inside) >= 2 else np.nan
    return out


def empirical_p(p, null):
    """(1 + #{null <= p}) / (1 + N) for each p (a calibrated p against a pooled null sample)."""
    null = np.sort(np.asarray(null, float))
    null = null[np.isfinite(null)]
    if not len(null):
        raise ValueError('Empty null sample.')
    p = np.asarray(p, float)
    return (1 + np.searchsorted(null, p, side='right')) / (1 + len(null))


def equal_weight_quantiles(values, groups, q):
    """Quantiles of ``values`` where every group carries the same total weight (each unit 1 / its group size)."""
    v = np.asarray(values, float)
    g = pd.Series(np.asarray(groups).astype(str))
    w = (1 / g.map(g.value_counts())).to_numpy()
    order = np.argsort(v, kind='stable')
    cum = np.cumsum(w[order]) - .5 * w[order]
    cum /= w.sum()
    return np.interp(np.asarray(q, float), cum, v[order])


def knn_place(query, reference, values, labels, k=30):
    """For each query row, its k nearest reference rows (Euclidean): the median of their ``values`` and the share of
    each label. Returns a DataFrame with 'median_value', 'share_<label>' columns and 'majority' (the most frequent
    label; ties go to the label first in sorted order)."""
    from scipy.spatial import cKDTree

    ref = np.asarray(reference, float)
    qry = np.asarray(query, float)
    vals = np.asarray(values, float)
    labs = np.asarray(labels).astype(str)
    if ref.ndim != 2 or qry.ndim != 2 or ref.shape[1] != qry.shape[1] or len(vals) != len(ref) or len(labs) != len(ref):
        raise ValueError('Reference, values and labels must align, in the same space as the query.')
    if not 0 < k <= len(ref):
        raise ValueError('k must be between 1 and the number of reference rows.')
    _, idx = cKDTree(ref).query(qry, k=k)
    idx = np.asarray(idx).reshape(len(qry), k)
    out = pd.DataFrame({'median_value': np.median(vals[idx], axis=1)})
    names = np.unique(labs)
    shares = np.stack([(labs[idx] == name).mean(axis=1) for name in names], axis=1)
    for j, name in enumerate(names):
        out[f'share_{name}'] = shares[:, j]
    out['majority'] = names[np.argmax(shares, axis=1)]
    return out
