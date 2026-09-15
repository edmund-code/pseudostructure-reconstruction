"""Specimen-level summaries: balanced curves and pseudobulk profiles.

The pooled GAMs weight tubules, so a specimen contributes more wherever it happens to have more
structures and the specimen mixture can drift along pseudospace. Two summaries put the specimen
back in charge:

* :func:`specimen_balanced_curves` averages the per-specimen fitted curves with equal weight,
* :func:`pseudobulk_profiles` collapses raw counts per (specimen, pseudospace bin), which is the
  input a count-based or specimen-level analysis expects.

Bins from the same specimen are repeated measurements along one coordinate, not replicates; the
number of independent units is the number of specimens.
"""
from __future__ import annotations

import warnings
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

__all__ = ['specimen_balanced_curves', 'pseudobulk_profiles', 'coordinate_agreement']


def specimen_balanced_curves(curves_by_specimen: Mapping[str, np.ndarray]) -> np.ndarray:
    """Equal-weight mean of per-specimen curves, one column per grid point.

    Each specimen contributes 1/n regardless of how many structures it has, and a specimen that has
    no support at a grid point (NaN) simply does not vote there. Returns NaN where no specimen has
    support.
    """
    if not curves_by_specimen:
        raise ValueError('no specimen curves supplied')
    stack = np.stack([np.asarray(curves, dtype=float) for curves in curves_by_specimen.values()])
    # A grid point where no specimen has support is legitimately all-NaN; report NaN, do not warn.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.nanmean(stack, axis=0)


def pseudobulk_profiles(counts, specimen: Sequence[str], pseudospace: Sequence[float], *,
                        n_bins: int = 10, gene_names: Sequence[str] | None = None,
                        min_structures_per_bin: int = 5) -> pd.DataFrame:
    """Sum raw counts within each (specimen, pseudospace bin).

    ``counts`` is a (structures x genes) sparse or dense matrix aligned with ``specimen`` and
    ``pseudospace``. Returns a long-form table with one row per (specimen, bin) that holds at least
    ``min_structures_per_bin`` structures, together with the number of structures and the summed
    counts of every gene. Use the specimen identifier as the experimental unit downstream.
    """
    matrix = counts.tocsr() if hasattr(counts, 'tocsr') else np.asarray(counts, dtype=float)
    specimen = np.asarray(specimen).astype(str)
    pseudospace = np.asarray(pseudospace, dtype=float)
    if matrix.shape[0] != specimen.size or specimen.size != pseudospace.size:
        raise ValueError('counts, specimen and pseudospace must share their first axis')
    genes = list(gene_names) if gene_names is not None else [f'gene_{i}' for i in range(matrix.shape[1])]

    finite = np.isfinite(pseudospace)
    edges = np.linspace(0.0, 1.0, int(n_bins) + 1)
    codes = np.full(pseudospace.size, -1, dtype=int)
    codes[finite] = np.clip(np.digitize(pseudospace[finite], edges) - 1, 0, int(n_bins) - 1)

    rows = []
    for name in sorted(set(specimen[finite])):
        for bin_index in range(int(n_bins)):
            mask = (specimen == name) & (codes == bin_index)
            n_structures = int(mask.sum())
            if n_structures < int(min_structures_per_bin):
                continue
            block = matrix[mask]
            summed = np.asarray(block.sum(axis=0)).ravel()
            entry = {
                'specimen': name,
                'pseudospace_bin': bin_index,
                'pseudospace_bin_centre': float((edges[bin_index] + edges[bin_index + 1]) / 2.0),
                'n_structures': n_structures,
            }
            entry.update({f'count_{gene}': float(value) for gene, value in zip(genes, summed)})
            rows.append(entry)
    return pd.DataFrame(rows)


def coordinate_agreement(reference: Sequence[float], alternative: Sequence[float]) -> dict:
    """Summarise how much two versions of the same coordinate agree.

    Used for coordinate-robustness checks (feature exclusion, specimen omission): Pearson and
    Spearman correlation over finite pairs, the median absolute difference, and the fraction of
    pairs ordered identically, which is what a rank-based pseudospace actually relies on.
    """
    from scipy.stats import pearsonr, spearmanr

    left = np.asarray(reference, dtype=float)
    right = np.asarray(alternative, dtype=float)
    finite = np.isfinite(left) & np.isfinite(right)
    if finite.sum() < 3:
        return {'n_compared': int(finite.sum()), 'pearson': np.nan, 'spearman': np.nan,
                'median_abs_difference': np.nan, 'concordant_pairs': np.nan}
    left, right = left[finite], right[finite]
    # Order concordance over distinct pairs only: a structural zero from comparing a cell with
    # itself would otherwise count as agreement.
    off_diagonal = ~np.eye(left.size, dtype=bool)
    concordant = (np.sign(left[:, None] - left[None, :])
                  == np.sign(right[:, None] - right[None, :]))
    return {
        'n_compared': int(finite.sum()),
        'pearson': float(pearsonr(left, right)[0]),
        'spearman': float(spearmanr(left, right).correlation),
        'median_abs_difference': float(np.median(np.abs(left - right))),
        'concordant_pairs': float(concordant[off_diagonal].mean()),
    }
