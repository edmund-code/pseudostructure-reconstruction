"""Positional programs along the PT against a cluster-level view (notebook 89).

Protocol: ``docs/results/pt-positional-programs.md``.

- ``fixed_bins``, ``bin_profiles``, ``reference_profile``, ``program_profiles``: binned gene and program profiles in
  mouse segment units (no spline).
- ``boundary_excess``: how far a profile's window at a segment boundary rises above (peak) or falls below (trough)
  both flanks; a monotone profile gives at most noise.
- ``strata_permutations``: H1's gene-label permutation within expression strata.
- ``redundancy_groups``, ``cluster_profiles``: E1's redundancy groups and shape clusters.
- ``bump``, ``ramp``, ``matched_step``, ``simulate_counts``: T's semi-synthetic programs.
- ``segment_call``: what the cluster view calls from two mice's segment profiles.
- ``centre_by_median``, ``pair_r2``, ``thin_to``: E4's Δ profiles and coordination along the PT.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from .transition_sharpness import binned_log_rates

EFFECT = .25              # log2; the project's effect threshold (76's O5)
WIDTH, LO, HI = .1, 0., 2.5
MIN_UNITS = 20
HALF, FLANK = .25, .5     # boundary window half-width and flank length (segment units)
MIN_BINS = 2
SEGMENTS = ('S1', 'S2', 'S3')


def fixed_bins(units, *, width=WIDTH, lo=LO, hi=HI):
    """Bin index of each unit on fixed ``width`` bins over [lo, hi) (-1 outside), and the bin centres."""
    units = np.asarray(units, float)
    edges = np.arange(lo, hi + width / 2, width)
    centres = (edges[:-1] + edges[1:]) / 2
    finite = np.isfinite(units)
    idx = np.full(len(units), -1)
    idx[finite] = np.floor((units[finite] - lo) / width + 1e-9).astype(int)
    idx[(idx < 0) | (idx >= len(centres))] = -1
    return idx, centres


def bin_profiles(counts, library, bins, n_bins, specimen, names, *, min_units=MIN_UNITS, pseudo=.5, centre=True):
    """Per specimen: genes x bins log2 rates (85's ``binned_log_rates``), minus each gene's mean over valid bins when
    ``centre``. A bin with fewer than ``min_units`` units of that specimen is NaN. Returns {name: array}."""
    counts = sparse.csr_matrix(counts)
    out = {}
    for name in names:
        rows = np.flatnonzero((np.asarray(specimen) == name) & (np.asarray(bins) >= 0))
        b = np.asarray(bins)[rows]
        y = binned_log_rates(counts[rows], np.asarray(library, float)[rows], b, n_bins, pseudo=pseudo)
        valid = np.bincount(b, minlength=n_bins) >= min_units
        y[:, ~valid] = np.nan
        out[name] = y - np.nanmean(y, axis=1, keepdims=True) if centre else y
    return out


def reference_profile(profile, usable):
    """The typical gene's profile: the median over ``usable`` rows, per bin."""
    return np.nanmedian(np.asarray(profile, float)[np.asarray(usable, bool)], axis=0)


def program_profiles(deviation, members, *, min_members=10):
    """Median over member rows of ``deviation`` (gene profile minus reference), per bin.

    ``members`` is a list of row-index arrays; programs with fewer than ``min_members`` rows are NaN. Returns
    programs x bins."""
    D = np.asarray(deviation, float)
    out = np.full((len(members), D.shape[1]), np.nan)
    for i, rows in enumerate(members):
        rows = np.asarray(rows, int)
        if len(rows) >= min_members:
            out[i] = np.nanmedian(D[rows], axis=0)
    return out


def _part_mean(P, mask, min_bins):
    sub = P[:, mask]
    n = np.isfinite(sub).sum(axis=1)
    with np.errstate(invalid='ignore'):
        m = np.nanmean(sub, axis=1) if sub.shape[1] else np.full(len(P), np.nan)
    m[n < min_bins] = np.nan
    return m


def boundary_excess(profiles, centres, k, *, half=HALF, flank=FLANK, hi=HI, min_bins=MIN_BINS):
    """Peak excess E = mean(window) - max(mean(left), mean(right)) and trough excess T = min(flanks) - mean(window) at
    boundary ``k`` (segment units), per profile row. Window [k - half, k + half); flanks [k - half - flank, k - half)
    and [k + half, k + half + flank) clipped at ``hi``. NaN when a part has fewer than ``min_bins`` valid bins."""
    P = np.atleast_2d(np.asarray(profiles, float))
    c = np.asarray(centres, float)
    w = _part_mean(P, (c >= k - half) & (c < k + half), min_bins)
    left = _part_mean(P, (c >= k - half - flank) & (c < k - half), min_bins)
    right = _part_mean(P, (c >= k + half) & (c < min(k + half + flank, hi)), min_bins)
    return w - np.maximum(left, right), np.minimum(left, right) - w


def strata_permutations(strata, n_perm, rng):
    """``n_perm`` x len(strata) arrays: each row maps a gene index to a gene of the same stratum (a permutation within
    every stratum)."""
    strata = np.asarray(strata)
    out = np.tile(np.arange(len(strata)), (n_perm, 1))
    for s in np.unique(strata):
        idx = np.flatnonzero(strata == s)
        for p in range(n_perm):
            out[p, idx] = rng.permutation(idx)
    return out


def redundancy_groups(sets, *, threshold=.5):
    """Connected components of sets linked by Jaccard >= ``threshold``. ``sets`` is {name: iterable}; returns a
    Series name -> group number (groups numbered by first appearance)."""
    names = list(sets)
    members = [set(sets[n]) for n in names]
    parent = list(range(len(names)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            union = len(members[i] | members[j])
            if union and len(members[i] & members[j]) / union >= threshold:
                parent[find(i)] = find(j)
    roots = [find(i) for i in range(len(names))]
    order = {r: k for k, r in enumerate(dict.fromkeys(roots))}
    return pd.Series([order[r] for r in roots], index=names, name='group')


def correlation_distance(profiles):
    """1 - Pearson r between rows (finite columns only)."""
    X = np.asarray(profiles, float)
    X = X[:, np.isfinite(X).all(axis=0)]
    r = np.corrcoef(X)
    d = np.clip(1 - r, 0, 2)
    np.fill_diagonal(d, 0)
    return (d + d.T) / 2


def cluster_profiles(profiles, *, k_range=(3, 10), k=None):
    """Average-linkage clusters of rows on 1 - r. With ``k`` None, K in ``k_range`` maximising the mean silhouette.
    Returns (labels 1..K, K, {K: silhouette})."""
    from sklearn.metrics import silhouette_score
    d = correlation_distance(profiles)
    tree = linkage(squareform(d, checks=False), method='average')
    if k is not None:
        return fcluster(tree, k, criterion='maxclust'), int(k), {}
    scores = {}
    for kk in range(k_range[0], min(k_range[1], len(d) - 1) + 1):
        lab = fcluster(tree, kk, criterion='maxclust')
        if len(np.unique(lab)) > 1:
            scores[kk] = float(silhouette_score(d, lab, metric='precomputed'))
    best = max(scores, key=lambda kk: (scores[kk], -kk))
    return fcluster(tree, best, criterion='maxclust'), best, scores


def bump(u, centre, width):
    """Gaussian of full width at half maximum ``width`` and height 1."""
    return np.exp(-4 * np.log(2) * (np.asarray(u, float) - centre) ** 2 / width ** 2)


def ramp(u, centre, width):
    """Logistic from 0 to 1 with 10–90% width ``width`` and half-point ``centre``."""
    return 1 / (1 + np.exp(-(np.asarray(u, float) - centre) * np.log(81) / width))


def matched_step(values, label, specimen):
    """Each unit's value replaced by the mean of ``values`` over its specimen's units with the same label: the cluster
    view of the result equals that of ``values`` by construction."""
    frame = pd.DataFrame({'v': np.asarray(values, float), 'l': np.asarray(label), 's': np.asarray(specimen)})
    return frame.groupby(['s', 'l']).v.transform('mean').to_numpy()


def simulate_counts(library, rate, shape, effect, alpha, rng):
    """NB counts (gamma–Poisson) for genes x units: mean library x rate x exp(shape + effect), rescaled per gene so
    the expected total equals library-weighted rate. ``rate`` (genes), ``shape`` (units, natural log), ``effect``
    (genes x units, natural log), ``alpha`` (genes). Returns units x genes."""
    lib = np.asarray(library, float)
    lin = np.exp(np.asarray(shape, float)[None, :] + np.asarray(effect, float))
    mu = lin / (lin * lib[None, :]).sum(axis=1, keepdims=True) * (lib.sum() * np.asarray(rate, float))[:, None]
    mu = mu * lib[None, :]
    a = np.maximum(np.asarray(alpha, float), 1e-8)[:, None]
    lam = rng.gamma(1 / a, mu * a)
    return rng.poisson(lam).T


def segment_call(values, *, effect=EFFECT, segments=SEGMENTS):
    """The cluster view's call from specimens x segments values (log2): the segments >= ``effect`` in every specimen,
    as 'none', 'S2' or 'S1+S2'."""
    v = np.asarray(values, float)
    hit = [s for j, s in enumerate(segments) if np.all(v[:, j] >= effect)]
    return '+'.join(hit) if hit else 'none'


def centre_by_median(delta, usable):
    """Subtract, per bin, the median over ``usable`` rows (76's ``centre_rows``)."""
    D = np.asarray(delta, float)
    return D - np.nanmedian(D[np.asarray(usable, bool)], axis=0)


def thin_to(counts, library, target, rng):
    """Binomial thinning of summed counts to library ``target`` (target <= library)."""
    if not library > 0:
        raise ValueError('Need a positive library to thin.')
    p = min(1., float(target) / float(library))
    return rng.binomial(np.asarray(counts, np.int64), p)


def pair_r2(a, b):
    """Squared Pearson correlation of two vectors over finite entries (NaN if < 3)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return np.nan
    return float(np.corrcoef(a[ok], b[ok])[0, 1] ** 2)
