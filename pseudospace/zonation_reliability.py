"""Same-species reliability ceilings and symmetric flat confirmation for zonation (notebook 43).

A cross-species correlation of gene gradients means little without the correlation two datasets
of the *same* species reach for the same contrast. Platform, dissection and label definitions all
lower it even when the biology is identical. Each correlation is first disattenuated for sampling
noise, using within-dataset replicate variance. The conservation index then divides the
cross-species value by the within-species cross-dataset ceilings:

    index = r*(human, mouse) / sqrt(r*(human, human') * r*(mouse, mouse'))

Here r* is the noise-corrected correlation, and human' and mouse' are independent datasets of the
same species. An index near 1 means the species agree as well as the data allow; near 0 means
different genes are zonated.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def paired_gradient(per_replicate):
    """Mean gradient and its sampling variance from replicates x genes contrasts (paired design)."""
    frame = pd.DataFrame(per_replicate)
    n = frame.shape[0]
    if n < 2:
        raise ValueError('Need at least two replicates.')
    return frame.mean(), frame.var(ddof=1) / n


def unpaired_gradient(late, early):
    """Mean gradient and sampling variance when late and early come from different replicates."""
    late, early = pd.DataFrame(late), pd.DataFrame(early)
    if late.shape[0] < 2 or early.shape[0] < 2:
        raise ValueError('Need at least two replicates per segment.')
    return late.mean() - early.mean(), late.var(ddof=1) / late.shape[0] + early.var(ddof=1) / early.shape[0]


def winsorize(values, q=.01):
    values = pd.Series(values, dtype=float)
    low, high = values.quantile([q, 1 - q])
    return values.clip(low, high)


def reliability(mean, noise):
    """Share of across-gene variance in the mean gradient that is not sampling noise (floored at 0)."""
    total = float(np.var(mean, ddof=1))
    return max(1 - float(np.mean(noise)) / total, 0.) if total > 0 else np.nan


def disattenuated(a, b, genes, *, q=.01):
    """Raw and noise-corrected Pearson correlation of two (mean, noise) gradients over genes."""
    (ma, na), (mb, nb) = a, b
    x, y = winsorize(ma.loc[genes], q), winsorize(mb.loc[genes], q)
    r = float(np.corrcoef(x, y)[0, 1])
    rel_a, rel_b = reliability(ma.loc[genes], na.loc[genes]), reliability(mb.loc[genes], nb.loc[genes])
    corrected = r / np.sqrt(rel_a * rel_b) if rel_a > 0 and rel_b > 0 else np.nan
    return {'r': r, 'reliability_a': rel_a, 'reliability_b': rel_b, 'r_corrected': corrected,
            'spearman': float(stats.spearmanr(ma.loc[genes], mb.loc[genes]).statistic), 'genes': len(genes)}


def conservation_index(human, mouse, human_ref, mouse_ref, genes, *, n_boot=1000, seed=0, q=.01):
    """Cross-species noise-corrected correlation divided by the two within-species ceilings.

    Returns the point estimate, its components and a percentile interval from resampling genes.
    """
    genes = pd.Index(genes)

    def compute(g):
        cross = disattenuated(human, mouse, g, q=q)['r_corrected']
        ceiling_h = disattenuated(human, human_ref, g, q=q)['r_corrected']
        ceiling_m = disattenuated(mouse, mouse_ref, g, q=q)['r_corrected']
        index = cross / np.sqrt(ceiling_h * ceiling_m) if ceiling_h > 0 and ceiling_m > 0 else np.nan
        return index, cross, ceiling_h, ceiling_m

    point = compute(genes)
    rng = np.random.default_rng(seed)
    draws = np.array([compute(genes[rng.integers(0, len(genes), len(genes))]) for _ in range(n_boot)])
    low, high = np.nanpercentile(draws[:, 0], [2.5, 97.5])
    return {'index': point[0], 'ci_low': low, 'ci_high': high, 'cross_corrected': point[1],
            'ceiling_human': point[2], 'ceiling_mouse': point[3], 'genes': len(genes)}


def confirm_flat(classes, calls_mouse, calls_human, reference_mouse, reference_human, *, alpha=.05, effect=.5):
    """Require the 'flat' side of a species-only class to be flat-compatible in its external atlas.

    For a mouse-only gene, the human reference must not be significant in the mouse direction
    (one-sided p > alpha), and its |mean| must be below ``effect``; human-only is symmetric. Genes
    that fail, or are missing from the reference, become 'indeterminate'.
    """
    out = classes.copy()
    for klass, zonated_calls, reference in (('mouse-only', calls_mouse, reference_human),
                                            ('human-only', calls_human, reference_mouse)):
        genes = out.index[out.eq(klass)]
        ref = reference.reindex(genes)
        direction = zonated_calls.reindex(genes)
        p_up, p_down = stats.t.sf(ref.t, ref.df), stats.t.cdf(ref.t, ref.df)
        same_direction = np.where(direction.eq('up'), p_up, p_down) <= alpha
        ok = ~same_direction & (ref['mean'].abs() < effect) & ref['mean'].notna()
        out.loc[genes[~ok.to_numpy()]] = 'indeterminate'
    return out
