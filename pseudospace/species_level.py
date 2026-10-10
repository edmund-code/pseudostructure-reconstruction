"""The constant (level) human-vs-mouse difference as a positive control (notebook 68).

Protocol: the "Notebook 68 protocol" section of the results notes. Level differences are confounded with
per-gene probe efficiency (the species use different probe panels) and, with one human donor, have no species
replication, so these helpers measure agreement and effect sizes and never test a species difference on our data.

* ``grid_centred_level``: L_g = (1 / ln 2) * mean over grid x of [Delta_g(x) - m(x)], m(x) the median of Delta at x
  over reference genes. Delta is a natural-log rate ratio on a common grid (``pathway_pipelines.GeneModel.delta``);
  the result is in log2 units relative to the typical gene. ``pair_levels`` applies it to specimen-curve pairs.
* ``segment_pseudobulk`` / ``equal_segment_log_cpm``: raw-count sums per unit (specimen or donor) x segment, and the
  equal-weight mean over segments of voom-style log-CPM, log2((y + 0.5) / (N + 1) * 1e6). The same construction
  serves our specimens (P_g), the atlas donors and family sums (``family_sums``).
* ``limma_trend_species``: limma lmFit + eBayes(trend = TRUE) of donor-level values on a two-group design (R, through
  rpy2), with the coefficient centred at its median over reference genes and tested against that median.
* ``band``, ``clears_band``, ``wilson_interval`` and ``bootstrap_spearman``: the yardsticks of the report.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import norm, rankdata
from scipy.stats import t as t_dist

from .stats_gam import bh_adjust

LN2 = np.log(2.0)
PRIOR_COUNT = .5


def _reference_mask(n, reference):
    mask = np.ones(n, bool) if reference is None else np.asarray(reference, bool)
    if mask.shape != (n,):
        raise ValueError('The reference mask needs one entry per gene.')
    return mask


def grid_centred_level(delta, *, reference=None, weights=None):
    """Per gene, the weighted grid mean of Delta_g(x) - m(x), divided by ln 2 (log2 units).

    ``delta``: genes x grid natural-log ratios. ``reference``: genes whose median at each grid point is the typical
    gene (default: all); reference rows with a non-finite value are ignored there. ``weights``: one per grid point,
    normalised to sum 1 (default uniform). A gene with any non-finite Delta gets NaN. Returns a Series when ``delta``
    is a DataFrame, else an array.
    """
    index = delta.index if isinstance(delta, pd.DataFrame) else None
    d = np.asarray(delta, float)
    if d.ndim != 2 or d.shape[1] < 1:
        raise ValueError('Delta must be genes x grid.')
    ref = _reference_mask(d.shape[0], reference)
    w = np.full(d.shape[1], 1.0 / d.shape[1]) if weights is None else np.asarray(weights, float)
    if w.shape != (d.shape[1],) or not np.isfinite(w).all() or (w < 0).any() or w.sum() <= 0:
        raise ValueError('Need one nonnegative finite weight per grid point.')
    w = w / w.sum()
    block = np.where(np.isfinite(d[ref]), d[ref], np.nan)
    if not np.isfinite(block).any(axis=0).all():
        raise ValueError('Every grid point needs at least one finite reference gene.')
    centre = np.nanmedian(block, axis=0)
    finite = np.isfinite(d).all(axis=1)
    level = np.full(d.shape[0], np.nan)
    level[finite] = ((d[finite] - centre) @ w) / LN2
    return pd.Series(level, index=index) if index is not None else level


def pair_levels(curves, pairs, *, valid=None, reference=None, weights=None, index=None):
    """``grid_centred_level`` of curve_a - curve_b for each pair (a, b), centred by its own median.

    ``curves``: {specimen: genes x grid natural-log rates}. ``valid``: {specimen: bool per gene} (converged and not
    separated in that specimen's fit; default all True). The pair's median uses reference genes valid in both
    specimens; a gene invalid in either specimen gets NaN. Returns a DataFrame genes x 'a - b'.
    """
    names = {n for pair in pairs for n in pair}
    missing = names - set(curves)
    if missing:
        raise ValueError(f'No curves for {sorted(missing)}.')
    shapes = {np.shape(curves[n]) for n in names}
    if len(shapes) != 1:
        raise ValueError('All specimen curves must have the same genes x grid shape.')
    n_genes = next(iter(shapes))[0]
    ok = {n: np.ones(n_genes, bool) if valid is None else np.asarray(valid[n], bool) for n in names}
    ref = _reference_mask(n_genes, reference)
    out = {}
    for a, b in pairs:
        both = ok[a] & ok[b]
        level = grid_centred_level(np.asarray(curves[a], float) - np.asarray(curves[b], float),
                                   reference=ref & both, weights=weights)
        out[f'{a} - {b}'] = np.where(both, level, np.nan)
    return pd.DataFrame(out, index=index)


def pair_summary(levels):
    """Per gene: minimum, maximum, number of pairs present and whether all pairs are present and share the sign."""
    values = np.asarray(levels, float)
    present = np.isfinite(values)
    n = present.sum(axis=1)
    complete = n == values.shape[1]
    with np.errstate(invalid='ignore'):
        low = np.where(n > 0, np.nanmin(np.where(present, values, np.inf), axis=1), np.nan)
        high = np.where(n > 0, np.nanmax(np.where(present, values, -np.inf), axis=1), np.nan)
    same_sign = complete & ((low > 0) | (high < 0))
    index = levels.index if isinstance(levels, pd.DataFrame) else None
    return pd.DataFrame({'pair_min': low, 'pair_max': high, 'pairs': n, 'pairs_same_sign': same_sign}, index=index)


def segment_pseudobulk(counts, library, unit, segment, *, segments=('S1', 'S2', 'S3')):
    """Raw-count and library sums per unit x segment.

    ``counts``: structures x genes (dense or sparse) raw counts; ``library``: per structure; ``unit``/``segment``:
    labels per structure. Every unit must have every segment. Returns (units, sums: units x segments x genes,
    totals: units x segments), units sorted.
    """
    unit, segment = np.asarray(unit).astype(str), np.asarray(segment).astype(str)
    library = np.asarray(library, float)
    n = counts.shape[0]
    if unit.shape != (n,) or segment.shape != (n,) or library.shape != (n,):
        raise ValueError('Counts, library, unit and segment must align.')
    units = sorted(set(unit))
    sums = np.zeros((len(units), len(segments), counts.shape[1]))
    totals = np.zeros((len(units), len(segments)))
    matrix = sparse.csr_matrix(counts) if sparse.issparse(counts) else np.asarray(counts, float)
    for i, u in enumerate(units):
        for j, s in enumerate(segments):
            rows = (unit == u) & (segment == s)
            if not rows.any():
                raise ValueError(f'Unit {u} has no structure in segment {s}.')
            sums[i, j] = np.asarray(matrix[rows].sum(axis=0)).ravel()
            totals[i, j] = library[rows].sum()
    return units, sums, totals


def log_cpm(sums, totals, prior_count=PRIOR_COUNT):
    """voom-style log2 counts per million: log2((y + prior) / (N + 1) * 1e6); ``totals`` broadcast over genes."""
    y, n = np.asarray(sums, float), np.asarray(totals, float)
    if (y < 0).any() or not np.isfinite(y).all() or (n <= 0).any() or not np.isfinite(n).all():
        raise ValueError('Need nonnegative finite counts and positive finite totals.')
    return np.log2((y + prior_count) / (n[..., None] + 1.0) * 1e6)


def equal_segment_log_cpm(sums, totals, *, weights=None, prior_count=PRIOR_COUNT):
    """Per unit, the weighted mean over segments of ``log_cpm`` (default equal weights). Returns units x genes."""
    values = log_cpm(sums, totals, prior_count)
    k = values.shape[1]
    w = np.full(k, 1.0 / k) if weights is None else np.asarray(weights, float)
    if w.shape != (k,) or (w < 0).any() or w.sum() <= 0:
        raise ValueError('Need one nonnegative weight per segment.')
    return np.einsum('usg,s->ug', values, w / w.sum())


def species_difference(values, case):
    """Mean over case units minus mean over the other units (units x genes values). Uncentred."""
    values, case = np.asarray(values, float), np.asarray(case, bool)
    if values.ndim != 2 or case.shape != (values.shape[0],) or case.all() or not case.any():
        raise ValueError('Need units x genes values and a case flag per unit with both groups present.')
    return values[case].mean(axis=0) - values[~case].mean(axis=0)


def family_sums(counts, genes, families):
    """Structures x families raw-count sums over each family's member genes.

    ``counts``: structures x genes (dense or sparse); ``genes``: the column names; ``families``: {family: members}.
    A member absent from ``genes`` raises. Returns (sums as a CSR matrix, family names in ``families`` order).
    """
    genes = pd.Index(np.asarray(genes).astype(str))
    if genes.has_duplicates or len(genes) != counts.shape[1]:
        raise ValueError('Need one unique gene name per count column.')
    rows, cols = [], []
    for k, members in enumerate(families.values()):
        idx = genes.get_indexer(list(map(str, members)))
        if (idx < 0).any() or len(idx) == 0:
            raise ValueError(f'Family {list(families)[k]!r} has members outside the gene list or none.')
        rows.extend(idx)
        cols.extend([k] * len(idx))
    indicator = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(genes), len(families)))
    return sparse.csr_matrix(sparse.csr_matrix(counts, dtype=float) @ indicator), list(families)


_LIMMA_TREND = """
function(E, case, trend, robust) {
  design <- cbind(Intercept = 1, case = case)
  fit <- limma::eBayes(limma::lmFit(E, design), trend = trend, robust = robust)
  list(coef = unname(fit$coefficients[, 2]), su = unname(fit$stdev.unscaled[, 2]), sigma = unname(fit$sigma),
       s2post = unname(fit$s2.post), dftotal = unname(fit$df.total), dfresid = unname(fit$df.residual),
       amean = unname(fit$Amean), s2prior = unname(fit$s2.prior), dfprior = unname(fit$df.prior))
}
"""


def limma_trend_species(values, case, *, genes=None, reference=None, trend=True, robust=True):
    """limma's moderated two-group coefficient on units x genes values (R, through rpy2), centred and tested.

    ``values``: units (donors) x genes, e.g. ``equal_segment_log_cpm``; ``case``: bool per unit (group 1). Fits
    lmFit(~ case) and eBayes(trend, robust) (limma's defaults otherwise). The coefficient (log2 when the values are)
    equals the difference of group means. Centred = coefficient - its median over ``reference`` genes; t_centred =
    centred / (stdev.unscaled * sqrt(s2.post)) on df.total (two-sided p against the median, BH q over all genes).
    """
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr

    values, case = np.asarray(values, float), np.asarray(case, bool)
    if values.ndim != 2 or case.shape != (values.shape[0],) or case.sum() < 2 or (~case).sum() < 2:
        raise ValueError('Need units x genes values and at least two units per group.')
    if not np.isfinite(values).all():
        raise ValueError('Values must be finite.')
    n_units, n_genes = values.shape
    ref = _reference_mask(n_genes, reference)
    importr('limma')
    matrix = ro.r['matrix'](ro.FloatVector(values.T.ravel(order='F')), nrow=n_genes, ncol=n_units)
    fit = ro.r(_LIMMA_TREND)(matrix, ro.FloatVector(case.astype(float)), bool(trend), bool(robust))
    get = lambda key: np.array(np.broadcast_to(np.asarray(fit.rx2(key), float), (n_genes,)))   # noqa: E731
    out = pd.DataFrame({key: get(r_key) for key, r_key in (
        ('coefficient', 'coef'), ('stdev_unscaled', 'su'), ('sigma', 'sigma'), ('s2_post', 's2post'),
        ('df_total', 'dftotal'), ('df_residual', 'dfresid'), ('amean', 'amean'), ('s2_prior', 's2prior'),
        ('df_prior', 'dfprior'))}, index=genes)
    out['centred'] = out.coefficient - np.median(out.coefficient.to_numpy()[ref])
    out['se'] = out.stdev_unscaled * np.sqrt(out.s2_post)
    out['t_centred'] = out.centred / out.se
    df = out.df_total.to_numpy()
    tail = np.where(np.isinf(df), norm.sf(np.abs(out.t_centred)), t_dist.sf(np.abs(out.t_centred), np.where(np.isinf(df), 1, df)))
    out['p'] = 2 * tail
    out['q'] = bh_adjust(out.p.to_numpy())
    return out


def band(values, quantiles=(.025, .975)):
    """The (low, high) quantiles of the finite values."""
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        raise ValueError('No finite values for the band.')
    return tuple(float(x) for x in np.quantile(v, quantiles))


def clears_band(values, expected_sign, bounds):
    """True where a value lies beyond the band on its expected side (+1: above high; -1: below low)."""
    v, s = np.asarray(values, float), np.asarray(expected_sign, float)
    low, high = bounds
    return np.where(s > 0, v > high, np.where(s < 0, v < low, False))


def wilson_interval(k, n, level=.95):
    """Wilson score interval of a proportion k / n; (nan, nan) when n = 0."""
    if n == 0:
        return np.nan, np.nan
    z = norm.ppf(.5 + level / 2)
    p = k / n
    centre = (p + z ** 2 / (2 * n)) / (1 + z ** 2 / n)
    half = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / (1 + z ** 2 / n)
    return float(centre - half), float(centre + half)


def _spearman_rows(x, y):
    """Row-wise Spearman correlation of two equally shaped 2-D arrays (average ranks)."""
    rx, ry = rankdata(x, axis=1), rankdata(y, axis=1)
    rx -= rx.mean(axis=1, keepdims=True)
    ry -= ry.mean(axis=1, keepdims=True)
    return (rx * ry).sum(axis=1) / np.sqrt((rx ** 2).sum(axis=1) * (ry ** 2).sum(axis=1))


def bootstrap_spearman(x, y, *, n_boot=2000, seed=68, level=.95, chunk=200):
    """Spearman's rho of the finite pairs and its percentile interval over ``n_boot`` resamples of the pairs."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if x.size < 3:
        return {'rho': np.nan, 'low': np.nan, 'high': np.nan, 'n': int(x.size)}
    rho = float(_spearman_rows(x[None], y[None])[0])
    rng = np.random.default_rng(seed)
    draws = []
    for start in range(0, n_boot, chunk):
        idx = rng.integers(0, x.size, size=(min(chunk, n_boot - start), x.size))
        draws.append(_spearman_rows(x[idx], y[idx]))
    draws = np.concatenate(draws)
    low, high = np.nanquantile(draws, [(1 - level) / 2, (1 + level) / 2])
    return {'rho': rho, 'low': float(low), 'high': float(high), 'n': int(x.size)}
