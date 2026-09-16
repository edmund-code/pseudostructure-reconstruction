"""Curve-module discovery: positional programs and response programs from fitted GAM curves.

Each measured feature is a curve over the shared pseudospace grid, so a module is a group of features
that share a *shape*. Two catalogs are built from the same machinery:

* **positional** modules cluster the reference-condition curves (where is a gene expressed along PT),
* **response** modules cluster the fitted difference curves ``Delta f(s) = f_comparison - f_reference``
  (how a gene changes), optionally mean-centred so a broad suppression and a peak relocation that
  share an overall offset are not forced together.

Clustering is on the correlation distance between grid-standardised curves, so level and amplitude
do not decide membership; amplitude and peak descriptors annotate a module instead. Peaks are never
the grouping key: a narrow spike and a broad plateau both peaking in mid-PT are different programs.
Dynamic time warping is deliberately not used -- aligning peaks would erase the positional
difference the analysis is trying to detect.
"""
from __future__ import annotations

import warnings
from typing import Iterable, Mapping

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.stats import hypergeom, spearmanr
from sklearn.metrics import adjusted_rand_score

from .stats_gam import bh_adjust

__all__ = [
    'difference_curves',
    'curve_descriptors',
    'discover_curve_modules',
    'module_stability',
    'enrich_modules',
]


def difference_curves(curve_reference, curve_comparison, *, center=False):
    """``f_comparison(s) - f_reference(s)``; ``center=True`` removes the mean so only the pattern of
    change remains (broad suppression vs selective loss vs relocation)."""
    diff = np.asarray(curve_comparison, dtype=float) - np.asarray(curve_reference, dtype=float)
    if center:
        diff = diff - np.nanmean(diff, axis=1, keepdims=True)
    return diff


def curve_descriptors(curves, grid):
    """Describe each curve by peak position, width, monotonicity, amplitude and direction.

    All descriptors are annotations of a discovered module, not selection criteria.
    """
    values = np.asarray(curves, dtype=float)
    positions = np.asarray(grid, dtype=float)
    rows = []
    for index in range(values.shape[0]):
        curve = values[index]
        finite = np.isfinite(curve)
        if finite.sum() < 3:
            rows.append({
                'peak_position': np.nan, 'peak_value': np.nan, 'half_max_width': np.nan,
                'monotonicity': np.nan, 'amplitude': np.nan, 'direction': 'unassigned',
            })
            continue
        x = positions[finite]
        y = curve[finite]
        peak_position = float(x[int(np.argmax(y))])
        amplitude = float(np.max(y) - np.min(y))
        half = np.min(y) + amplitude / 2.0
        above = x[y >= half]
        half_width = float(above.max() - above.min()) if above.size else 0.0
        # A constant curve has no rank correlation; report NaN rather than a spurious 0. The
        # tolerance is scale-relative because float noise can make a constant curve non-constant.
        monotonicity = (float(spearmanr(x, y).correlation)
                        if y.size > 2 and np.ptp(y) > 1e-9 * max(1.0, float(np.max(np.abs(y))))
                        else np.nan)
        if amplitude < 1e-12:
            direction = 'flat'
        elif peak_position > 0.66 * x.max():
            direction = 'late-rising' if monotonicity and monotonicity > 0.5 else 'late-peaked'
        elif peak_position < 0.33 * x.max():
            direction = 'early-declining' if monotonicity and monotonicity < -0.5 else 'early-peaked'
        else:
            direction = 'mid-plateau' if half_width > 0.5 else 'mid-peaked'
        rows.append({
            'peak_position': peak_position, 'peak_value': float(np.max(y)),
            'half_max_width': half_width, 'monotonicity': monotonicity,
            'amplitude': amplitude, 'direction': direction,
        })
    return pd.DataFrame(rows)


def _standardize_curves(curves, min_column_finite_fraction=0.5):
    """Standardise each curve over the grid this run actually covers.

    A run whose curves come from one specimen is undefined outside that specimen's own pseudospace
    support, and demanding a fully finite row would then reject every feature - which is how a
    leave-one-specimen-out probe can report zero modules. Columns covered by at least
    ``min_column_finite_fraction`` of the features define the run's support; features are then
    required to be finite on that support.
    """
    values = np.asarray(curves, dtype=float)
    finite = np.isfinite(values)
    column_ok = finite.mean(axis=0) >= float(min_column_finite_fraction)
    if column_ok.sum() < 3:                     # too little overlap to describe any shape
        column_ok = np.ones(values.shape[1], dtype=bool)
    subset = values[:, column_ok]
    subset_finite = finite[:, column_ok]
    # An all-NaN row is a legitimate outcome of a narrow run; report it as unusable, do not warn.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        mean = np.nanmean(subset, axis=1, keepdims=True)
        spread = np.nanstd(subset, axis=1, keepdims=True)
    standardized = np.where(spread > 0, (subset - mean) / np.where(spread > 0, spread, 1.0), np.nan)
    usable = subset_finite.all(axis=1) & (spread.ravel() > 0)
    return standardized, usable, column_ok


def discover_curve_modules(curves, grid, *, n_modules=None, max_distance=0.5, min_amplitude=0.05,
                           method='average', min_features=3, feature_names=None,
                           min_column_finite_fraction=0.5):
    """Cluster curves into modules by shape, using correlation distance on standardised curves.

    ``n_modules`` (if given) cuts the tree at that many clusters, otherwise ``max_distance`` is the
    cut height on the 1-Pearson distance scale. Curves that are all-NaN, constant, or below
    ``min_amplitude`` are labelled ``unassigned`` instead of being forced into a module, and modules
    smaller than ``min_features`` are also reported as ``unassigned`` (they are not programs yet).
    Returns ``(labels, module_table)``.
    """
    values = np.asarray(curves, dtype=float)
    names = (list(feature_names) if feature_names is not None
             else [f'feature_{i}' for i in range(values.shape[0])])
    positions = np.asarray(grid, dtype=float)
    standardized, usable, column_ok = _standardize_curves(
        values, min_column_finite_fraction=min_column_finite_fraction)
    # Descriptors are read on the same support the curves were standardised on, so peak positions are
    # never reported from a stretch of pseudospace this run cannot see.
    descriptors = curve_descriptors(values[:, column_ok], positions[column_ok])
    descriptors.insert(0, 'feature', names)
    usable &= descriptors['amplitude'].to_numpy() >= min_amplitude

    labels = np.array(['unassigned'] * len(names), dtype=object)
    if usable.sum() >= 2:
        matrix = standardized[usable]
        corr = np.corrcoef(matrix)
        corr = np.nan_to_num(corr, nan=0.0)
        distance = np.clip(1.0 - corr, 0.0, 2.0)
        np.fill_diagonal(distance, 0.0)
        condensed = distance[np.triu_indices(distance.shape[0], k=1)]
        tree = linkage(condensed, method=method)
        if n_modules is not None:
            raw = fcluster(tree, t=int(n_modules), criterion='maxclust')
        else:
            raw = fcluster(tree, t=float(max_distance), criterion='distance')
        usable_index = np.flatnonzero(usable)
        counts = pd.Series(raw).value_counts()
        for cluster, position in zip(raw, usable_index):
            size = int(counts[cluster])
            labels[position] = f'M{int(cluster)}' if size >= min_features else 'unassigned'

    labels = pd.Series(labels, index=names, name='module')

    rows = []
    for module, members in labels.groupby(labels, observed=True):
        if module == 'unassigned':
            continue
        subset = descriptors.set_index('feature').loc[members.index]
        rows.append({
            'module': module,
            'n_features': int(len(members)),
            'median_peak_position': float(subset['peak_position'].median()),
            'median_half_max_width': float(subset['half_max_width'].median()),
            'median_monotonicity': float(subset['monotonicity'].median()),
            'median_amplitude': float(subset['amplitude'].median()),
            'dominant_direction': subset['direction'].mode().iloc[0] if len(subset) else '',
            'features': '; '.join(map(str, members.index)),
        })
    module_table = pd.DataFrame(rows).sort_values('n_features', ascending=False).reset_index(drop=True) \
        if rows else pd.DataFrame(columns=['module', 'n_features'])
    return labels, module_table


def module_stability(labels_by_run: Mapping[str, pd.Series]) -> pd.DataFrame:
    """Agreement of module assignments across runs (leave-one-specimen-out, other smoothing).

    Returns one row per run with the adjusted Rand index against the first run (label invariant), the
    fraction of features labelled at all, and ``mean_reference_module_retention``: for each reference
    module, the fraction of its members that land together in this run's best-overlapping module, then
    averaged over reference modules. Retention is the label-invariant replacement for counting
    identical module names, which is meaningless because module names are arbitrary.
    """
    runs = {name: pd.Series(labels).astype(str) for name, labels in labels_by_run.items()}
    reference_name = next(iter(runs))
    reference = runs[reference_name]
    reference_modules = [module for module in reference.unique() if module != 'unassigned']

    rows = []
    for name, labels in runs.items():
        shared = reference.index.intersection(labels.index)
        left, right = reference.loc[shared], labels.loc[shared]
        retentions = []
        for module in reference_modules:
            members = set(left.index[left == module])
            if not members:
                continue
            best = right.loc[sorted(members)].value_counts()
            if best.empty:
                continue
            best_module = best.index[0]
            if best_module == 'unassigned':
                retentions.append(0.0)
                continue
            best_members = set(right.index[right == best_module])
            retentions.append(len(members & best_members) / len(members))
        rows.append({
            'run': name,
            'n_features': int(len(shared)),
            'n_modules': int(labels[labels != 'unassigned'].nunique()),
            'adjusted_rand_index_vs_reference': float(adjusted_rand_score(left, right)),
            'fraction_labelled': float((labels != 'unassigned').mean()),
            'mean_reference_module_retention': float(np.mean(retentions)) if retentions else np.nan,
        })
    return pd.DataFrame(rows)


def enrich_modules(module_members: Mapping[str, Iterable[str]], gene_sets: Mapping[str, Iterable[str]],
                   *, background: Iterable[str], correction='fdr_bh') -> pd.DataFrame:
    """Hypergeometric enrichment of every module against every gene set.

    ``background`` must be the genes that were *eligible for module discovery*, not the whole
    matrix: otherwise the test answers "is this set over-represented among expressed genes", which
    is already known. p-values are corrected across all tested (module, gene-set) pairs.
    """
    background_upper = {str(gene).upper() for gene in background}
    universe = len(background_upper)
    set_upper = {str(name): {str(gene).upper() for gene in members} & background_upper
                 for name, members in gene_sets.items()}
    rows = []
    for module, members in module_members.items():
        members = {str(gene).upper() for gene in members}
        usable = members & background_upper
        for name, genes in set_upper.items():
            hits = usable & genes
            if not genes or not usable:
                p_value = np.nan
            else:
                p_value = float(hypergeom.sf(len(hits) - 1, universe, len(genes), len(usable)))
            rows.append({
                'module': module,
                'gene_set': name,
                'n_module_genes': len(usable),
                'n_set_genes_in_background': len(genes),
                'n_overlap': len(hits),
                'overlap_genes': '; '.join(sorted(hits)),
                'expected_overlap': (len(usable) * len(genes) / universe) if universe else np.nan,
                'p_value': p_value,
            })
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame['p_value_adjusted'] = bh_adjust(frame['p_value'].to_numpy()) if correction == 'fdr_bh' \
        else frame['p_value']
    return frame.sort_values(['p_value_adjusted', 'n_overlap'], ascending=[True, False]).reset_index(drop=True)
