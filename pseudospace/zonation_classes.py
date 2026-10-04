"""Noise-aware per-gene zonation calls and cross-species zonation classes (notebook 40).

Each species is called separately, from per-specimen (or per-donor) segment contrasts. A
variance-moderated t scales every call to that species' own noise. The gene's between-replicate
variance is shrunk toward genes of similar expression, as in limma's squeezeVar, applied within
expression bins. A gene is:
- **zonated**: significant and at least ``effect`` in size;
- **flat**: its confidence interval lies inside ±``margin``;
- **indeterminate**: anything else.

A gene is therefore called mouse-only only when the human contrast is confidently small, not
merely non-significant. This keeps a noisier species from inflating the other species' "only"
class.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .split_plot import squeeze_var

CLASSES = ('conserved', 'reversal', 'mouse-only', 'human-only', 'indeterminate', 'neither')


def moderated_contrast(replicates, expression, *, n_bins=10):
    """Moderated one-sample t of the replicate mean, shrinking variances within expression bins.

    ``replicates``: DataFrame (replicates x genes) of contrasts; ``expression``: Series of a
    per-gene expression level used only to form bins. Returns a DataFrame with mean, se, t, df, p.
    """
    values = replicates.to_numpy(float)
    n = values.shape[0]
    if n < 2:
        raise ValueError('Need at least two replicates.')
    genes = replicates.columns
    mean = values.mean(axis=0)
    s2 = np.maximum(values.var(axis=0, ddof=1), 1e-12)
    level = pd.Series(expression, dtype=float).reindex(genes).to_numpy()
    bins = pd.qcut(pd.Series(level).rank(method='first'), min(n_bins, len(genes)), labels=False).to_numpy()
    posterior, total_df = np.empty_like(s2), np.empty_like(s2)
    for b in np.unique(bins):
        mask = bins == b
        post, d0, _ = squeeze_var(s2[mask], n - 1)
        posterior[mask] = post
        total_df[mask] = min(d0 + n - 1, 1e6)
    se = np.sqrt(posterior / n)
    t = mean / se
    p = 2 * stats.t.sf(np.abs(t), total_df)
    return pd.DataFrame({'mean': mean, 'se': se, 't': t, 'df': total_df, 'p': p}, index=genes)


def zonation_calls(contrast, *, alpha=.05, effect=.5, margin=.3, z=1.96):
    """'up', 'down', 'flat' or 'indeterminate' per gene from a moderated contrast table (BH over genes)."""
    q = multipletests(contrast.p, method='fdr_bh')[1]
    up = (q <= alpha) & (contrast['mean'] >= effect)
    down = (q <= alpha) & (contrast['mean'] <= -effect)
    flat = (contrast['mean'].abs() + z * contrast.se) < margin
    calls = np.select([up, down, flat], ['up', 'down', 'flat'], 'indeterminate')
    return pd.Series(calls, index=contrast.index), pd.Series(q, index=contrast.index)


def confirm(calls, reference, *, alpha=.05):
    """Keep a zonated call only if the same-species reference agrees (one-sided p <= alpha).

    ``reference``: moderated contrast table (mean, t, df) for the reference; missing genes lose the
    call (become 'indeterminate'). Flat and indeterminate calls are unchanged.
    """
    ref = reference.reindex(calls.index)
    one_sided_up = stats.t.sf(ref.t, ref.df)
    one_sided_down = stats.t.cdf(ref.t, ref.df)
    keep = ((calls == 'up') & (one_sided_up <= alpha)) | ((calls == 'down') & (one_sided_down <= alpha))
    zonated = calls.isin(['up', 'down'])
    out = calls.copy()
    out[zonated & ~keep.fillna(False)] = 'indeterminate'
    return out


def cross_species_class(mouse_call, human_call):
    """Combine two species' calls into a zonation class."""
    zon_m, zon_h = mouse_call in ('up', 'down'), human_call in ('up', 'down')
    if zon_m and zon_h:
        return 'conserved' if mouse_call == human_call else 'reversal'
    if zon_m:
        return 'mouse-only' if human_call == 'flat' else 'indeterminate'
    if zon_h:
        return 'human-only' if mouse_call == 'flat' else 'indeterminate'
    return 'neither' if (mouse_call == 'flat' or human_call == 'flat') else 'indeterminate'


def union_class(per_contrast):
    """One class per gene across contrasts, with precedence reversal > conserved > species-only.

    ``per_contrast``: DataFrame (genes x contrasts) of classes (NaN where not testable). A gene
    that is mouse-only in one contrast and human-only in another counts as 'reversal' only if a
    contrast is a reversal; otherwise its first species-only class, by contrast order, is kept.
    """
    order = ['reversal', 'conserved', 'mouse-only', 'human-only', 'indeterminate', 'neither']

    def pick(row):
        values = [v for v in row if isinstance(v, str)]
        for klass in order:
            if klass in values:
                return klass
        return np.nan

    return per_contrast.apply(pick, axis=1)
