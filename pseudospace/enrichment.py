"""Signed gene rankings and competitive (correlation-aware) set enrichment.

Three pathway questions need three different answers, and one unsigned ``shape_rms`` ranking cannot
give any of them:

1. *does this predefined program change along pseudospace* - answered by the fitted module/curve
   summaries,
2. *are its genes unusually affected relative to other genes* - needs a signed ranking and an
   explicit background, which is what this module provides,
3. *is a signalling pathway active* - transcriptional scores of pathway components do not establish
   activity; :func:`score_signed_signatures` scores signed downstream signatures instead (PROGENy's
   matrix is not bundled - supply it as a gene-to-weight mapping).

Competitive testing here works on a signed per-gene statistic. ``camera_like_enrichment`` reports
two p-values: a label-permutation p (random gene sets of the same size, which ignores within-set
correlation) and a correlation-aware z-test that inflates the variance of the set mean by the
average squared correlation of the set's genes measured in the expression matrix - the point CAMERA
makes, that genes in a pathway are not independent.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr

from .stats_gam import bh_adjust

__all__ = [
    'signed_gene_rankings',
    'camera_like_enrichment',
    'score_signed_signatures',
]


def _varies(values) -> bool:
    """True when a curve is not constant, using a scale-relative tolerance.

    A difference of two nearly identical fitted curves can be constant up to float noise; ranking
    such a curve would report a spurious correlation, so it is reported as undefined instead.
    """
    values = np.asarray(values, dtype=float)
    if values.size < 3:
        return False
    scale = max(1.0, float(np.max(np.abs(values))))
    return bool(np.ptp(values) > 1e-9 * scale)


def _as_matrix(expression):
    return expression.tocsr() if hasattr(expression, 'tocsr') else np.asarray(expression, dtype=float)


def signed_gene_rankings(curve_reference, curve_comparison, grid, *, gene_names=None,
                         thirds=3):
    """Signed per-gene statistics for the questions the unsigned metrics cannot separate.

    Returns one row per gene with

    * ``level_effect`` - mean of the fitted difference (a species/treatment offset),
    * ``amplitude_log2_ratio`` - gradient strength change,
    * ``redistribution`` - Spearman correlation between the difference curve and pseudospace
      (positive = the change shifts towards late pseudospace, a *signed* pattern statistic),
    * ``early_delta`` / ``mid_delta`` / ``late_delta`` - mean difference within each pseudospace
      third, so "lost the early program" and "gained late expression" stop being the same number.
    """
    reference = np.asarray(curve_reference, dtype=float)
    comparison = np.asarray(curve_comparison, dtype=float)
    delta = comparison - reference
    positions = np.asarray(grid, dtype=float)
    names = list(gene_names) if gene_names is not None else [f'feature_{i}' for i in range(delta.shape[0])]

    edges = np.quantile(positions, np.linspace(0.0, 1.0, int(thirds) + 1))
    masks = [(positions >= edges[i]) & (positions <= edges[i + 1]) for i in range(int(thirds))]
    frame = pd.DataFrame({
        'gene': names,
        'level_effect': np.nanmean(delta, axis=1),
        'amplitude_reference': np.nanmax(reference, axis=1) - np.nanmin(reference, axis=1),
        'amplitude_comparison': np.nanmax(comparison, axis=1) - np.nanmin(comparison, axis=1),
    })
    frame['amplitude_log2_ratio'] = np.where(
        (frame['amplitude_reference'] > 0) & (frame['amplitude_comparison'] > 0),
        np.log2(frame['amplitude_comparison']) - np.log2(frame['amplitude_reference']),
        np.nan,
    )
    frame['redistribution'] = [
        float(spearmanr(positions[np.isfinite(row)], row[np.isfinite(row)]).correlation)
        if np.isfinite(row).sum() > 2 and _varies(row[np.isfinite(row)]) else np.nan
        for row in delta
    ]
    for name, mask in zip(('early_delta', 'mid_delta', 'late_delta'), masks):
        with np.errstate(invalid='ignore'):
            frame[name] = np.nanmean(delta[:, mask], axis=1)
    return frame


def camera_like_enrichment(gene_statistics, gene_sets, *, expression=None, gene_names=None,
                           background=None, design=None, n_permutations=1000, seed=0,
                           min_set_size=5, correction='fdr_bh',
                           correlation_gene_cap=200, correlation_structure_cap=2000,
                           validate_with_permutations=False):
    """Competitive enrichment of signed gene statistics, with a correlation-aware variance.

    The test compares a set's mean statistic with the **background** mean, not with zero: the
    question is whether the members of this set behave differently from an average testable gene,
    which is what a competitive null means. For a set of ``k`` genes drawn from ``N``:

    * ``set_mean_statistic_centred`` = mean(member statistics) - mean(background statistics)
    * ``p_value_competitive_normal`` = two-sided normal tail of that mean under simple random
      sampling without replacement, i.e. sd = sd(background)/sqrt(k) * sqrt((N - k)/(N - 1)),
      multiplied by the variance inflation ``correlation_inflation`` estimated from the set's genes
    * ``p_value_permutation_with_replacement`` = an optional sanity check from random gene sets of
      the same size, which ignores both the finite-population correction and the correlation

    The variance inflation is the CAMERA idea - genes in a pathway are not independent - estimated
    from the **residual** correlations when a ``design`` is supplied (the fitted trend is removed, so
    shared pseudospace structure does not inflate it), and from mean-centred correlations otherwise.
    ``correlation_gene_cap`` members and ``correlation_structure_cap`` structures bound the cost.

    These are exploratory annotations. Nothing here is a confirmatory test: the design has two mice
    and one human donor, so the units of replication are specimens, not genes or sets.
    """
    statistics = dict(gene_statistics) if not isinstance(gene_statistics, pd.Series) else gene_statistics.to_dict()
    statistics = {str(gene).upper(): float(value) for gene, value in statistics.items()
                  if value is not None and np.isfinite(value)}
    allowed = ({str(gene).upper() for gene in background} if background is not None
               else set(statistics))

    ordered_genes = list(statistics)
    values = np.array([statistics[gene] for gene in ordered_genes], dtype=float)
    gene_index = {gene: index for index, gene in enumerate(ordered_genes)}
    background_mask = np.array([gene in allowed for gene in ordered_genes], dtype=bool)
    if background_mask.sum() < max(2, int(min_set_size)):
        raise ValueError('fewer testable background genes than the minimum set size')
    background_values = values[background_mask]
    background_mean = float(background_values.mean())
    background_sd = float(background_values.std(ddof=1))
    n_available = int(background_mask.sum())

    matrix = None
    column_by_gene = None
    if expression is not None and gene_names is not None:
        matrix = _as_matrix(expression)
        column_by_gene = {str(gene).upper(): index for index, gene in enumerate(gene_names)}

    rng = np.random.default_rng(seed)
    rows = []
    for name, members in gene_sets.items():
        present = np.array(sorted({gene_index[gene] for gene in
                                   {str(member).upper() for member in members}
                                   if gene in gene_index}), dtype=int)
        present = np.array([index for index in present if background_mask[index]], dtype=int)
        size = int(present.size)
        if size < int(min_set_size):
            rows.append({'gene_set': name, 'n_genes_tested': size,
                         'set_mean_statistic': np.nan, 'set_mean_statistic_centred': np.nan,
                         'variance_inflation': np.nan, 'correlation_basis': '',
                         'p_value_competitive_normal': np.nan,
                         'p_value_permutation_with_replacement': np.nan})
            continue
        observed = float(values[present].mean())
        centred = observed - background_mean

        inflation = 1.0
        basis = 'none'
        if matrix is not None:
            columns = [column_by_gene[ordered_genes[index]] for index in present
                       if ordered_genes[index] in column_by_gene]
            if len(columns) >= max(3, int(min_set_size)):
                capped = columns
                if len(capped) > int(correlation_gene_cap):
                    capped = list(rng.choice(capped, size=int(correlation_gene_cap), replace=False))
                block = (np.asarray(matrix[:, capped].todense()) if hasattr(matrix, 'todense')
                         else np.asarray(matrix[:, capped], dtype=float))
                if block.shape[0] > int(correlation_structure_cap):
                    rows_kept = rng.choice(block.shape[0],
                                           size=int(correlation_structure_cap), replace=False)
                    block = block[rows_kept]
                if design is not None:
                    design_block = np.asarray(design, dtype=float)
                    if design_block.shape[0] > block.shape[0]:
                        design_block = design_block[rows_kept] if 'rows_kept' in dir() else design_block
                    keep_rows = min(design_block.shape[0], block.shape[0])
                    design_block, block = design_block[:keep_rows], block[:keep_rows]
                    fitted = design_block @ (np.linalg.pinv(design_block) @ block)
                    block = block - fitted
                    basis = 'residual'
                else:
                    basis = 'mean-centred'
                block = block - block.mean(axis=0, keepdims=True)
                with np.errstate(invalid='ignore', divide='ignore'):
                    corr = np.corrcoef(block, rowvar=False)
                corr = np.nan_to_num(corr, nan=0.0)
                off_diagonal = ~np.eye(corr.shape[0], dtype=bool)
                if off_diagonal.any():
                    mean_squared = float(np.mean(corr[off_diagonal] ** 2))
                    inflation = 1.0 + (corr.shape[0] - 1) * mean_squared

        # competitive normal tail with the finite-population correction and the correlation inflation
        if background_sd > 0 and size < n_available:
            sd = background_sd / np.sqrt(size) * np.sqrt((n_available - size) / (n_available - 1))
            sd *= np.sqrt(max(inflation, 1.0))
            z_value = centred / sd if sd > 0 else np.nan
            p_competitive = 2 * (1 - norm.cdf(abs(z_value))) if np.isfinite(z_value) else np.nan
        else:
            p_competitive = np.nan

        p_permutation = np.nan
        if validate_with_permutations:
            draw = rng.choice(n_available, size=(int(n_permutations), size), replace=True)
            null_means = background_values[draw].mean(axis=1)
            tail = int(np.sum(np.abs(null_means - background_mean) >= abs(centred) - 1e-12))
            p_permutation = (tail + 1) / (int(n_permutations) + 1)

        rows.append({
            'gene_set': name, 'n_genes_tested': size, 'set_mean_statistic': observed,
            'set_mean_statistic_centred': centred, 'variance_inflation': inflation,
            'correlation_basis': basis, 'p_value_competitive_normal': float(p_competitive)
            if np.isfinite(p_competitive) else np.nan,
            'p_value_permutation_with_replacement': p_permutation,
        })
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame['p_value_competitive_normal_adjusted'] = (
        bh_adjust(frame['p_value_competitive_normal'].to_numpy()) if correction == 'fdr_bh'
        else frame['p_value_competitive_normal'])
    return frame.sort_values('p_value_competitive_normal_adjusted').reset_index(drop=True)


def score_signed_signatures(expression, signatures, *, gene_names=None, min_genes=5):
    """Score signed downstream signatures (e.g. a PROGENy matrix supplied by the caller).

    A signature maps a name to a gene -> weight mapping (positive = the signature is up when the gene
    is up, negative = the gene represses it). Genes are standardised across structures, so a score is
    comparable between signatures, and the score is the weighted mean over the genes actually present.
    Returns ``(scores, coverage)``: ``scores`` is (structures x signatures), ``coverage`` records how
    many members were present and whether the signature met ``min_genes``.
    """
    matrix = _as_matrix(expression)
    names = list(gene_names) if gene_names is not None else [f'gene_{i}' for i in range(matrix.shape[1])]
    index_by_gene = {str(gene).upper(): index for index, gene in enumerate(names)}
    scores = pd.DataFrame(index=pd.RangeIndex(matrix.shape[0]))
    coverage = []
    for signature, weights in signatures.items():
        entries = list(weights.items()) if isinstance(weights, dict) else [(gene, 1.0) for gene in weights]
        columns, signs = [], []
        for gene, weight in entries:
            index = index_by_gene.get(str(gene).upper())
            if index is not None:
                columns.append(index)
                signs.append(float(weight))
        coverage.append({
            'signature': signature,
            'n_members_requested': len(entries),
            'n_members_present': len(columns),
            'usable': len(columns) >= int(min_genes),
        })
        if len(columns) < int(min_genes):
            continue
        block = np.asarray(matrix[:, columns].todense()) if hasattr(matrix, 'todense') \
            else np.asarray(matrix[:, columns], dtype=float)
        mean = block.mean(axis=0, keepdims=True)
        spread = block.std(axis=0, keepdims=True)
        spread[spread == 0] = 1.0
        standardized = (block - mean) / spread
        weights = np.asarray(signs, dtype=float)
        scores[signature] = standardized @ (weights / np.abs(weights).sum())
    return scores, pd.DataFrame(coverage)
