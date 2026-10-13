"""Sharp or gradual: the width of a profile's changes along a coordinate (notebook 85).

Protocol: ``docs/results/pt-transition-sharpness.md``.

- ``equal_count_bins`` and ``binned_log_rates``: P1–P2, profiles on equal-count bins (no spline).
- ``fit_edges``: E1–E3, one or two logistic edges per profile, each with a midpoint and a 10–90% width.
- ``edge_class``, ``edge_location``, ``agreed_edges``: E4–E5 and R1–R2, the class and place of edges two specimens
  share.
- ``spline_step_response`` and ``crossing_width``: Q1, what a regression spline makes of a true step.
- ``edge_composition``: G2, the member genes behind a pathway edge.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import least_squares

from .continuous_descriptors import _crossing

LN81 = np.log(81.)       # 10–90% width of a logistic = ln(81) x its scale
W_BOUNDS = (.05, 6.)     # E1: widths in segment units
FLAT = .25               # E3: log2; an edge smaller than this does not count (84's flat threshold)
SHARP, GRADUAL = .5, 1.  # E4: 81's thresholds
BIC_MARGIN = 6.          # E2
MIN_PLATEAU = 3          # E2 (amended before any fit): bins between two M2 edges


def equal_count_bins(position, n_bins):
    """Bin index (0..n_bins-1) of each unit after sorting by ``position``: equal-count bins (P1)."""
    position = np.asarray(position, float)
    if n_bins < 2 or n_bins > len(position):
        raise ValueError('Need 2 <= n_bins <= number of units.')
    order = np.argsort(position, kind='stable')
    bins = np.empty(len(position), int)
    bins[order] = np.arange(len(position)) * n_bins // len(position)
    return bins


def binned_log_rates(counts, library, bins, n_bins, *, pseudo=.5):
    """Genes x bins of log2((sum of counts + pseudo) / sum of library) (P2). ``counts`` is units x genes."""
    member = sparse.csr_matrix((np.ones(len(bins)), (np.asarray(bins), np.arange(len(bins)))), shape=(n_bins, len(bins)))
    sums = member @ counts
    sums = sums.toarray() if sparse.issparse(sums) else np.asarray(sums)
    lib = member @ np.asarray(library, float)
    return (np.log2(sums + pseudo) - np.log2(lib)[:, None]).T


def _sig(x, m, w):
    return 1. / (1. + np.exp(-np.clip((x - m) * LN81 / w, -60, 60)))


def _design(x, params):
    """Intercept plus one logistic column per (m, w) pair."""
    return np.column_stack([np.ones(len(x))] + [_sig(x, params[i], params[i + 1]) for i in range(0, len(params), 2)])


def _rss_grid(x, Y, combos):
    """RSS of every profile (rows of Y) under every (m, w, ...) combination, amplitudes solved by least squares."""
    q = np.stack([np.linalg.qr(_design(x, c))[0] for c in combos])          # combos x bins x k
    yy = (Y ** 2).sum(axis=1)
    out = np.empty((len(combos), len(Y)))
    for start in range(0, len(Y), 500):
        proj = np.einsum('cnk,pn->cpk', q, Y[start:start + 500])
        out[:, start:start + 500] = yy[start:start + 500] - (proj ** 2).sum(axis=2)
    return np.maximum(out, 0.)


def _refine(x, y, start, lo, hi, w_bounds):
    """Variable projection: bounded least squares over (m, w, ...) with the amplitudes solved inside."""
    k = len(start) // 2
    lower = np.tile([lo, np.log(w_bounds[0])], k)
    upper = np.tile([hi, np.log(w_bounds[1])], k)
    theta0 = np.clip(np.array([v if i % 2 == 0 else np.log(v) for i, v in enumerate(start)]), lower + 1e-9, upper - 1e-9)

    def unpack(theta):
        return [v if i % 2 == 0 else np.exp(v) for i, v in enumerate(theta)]

    def resid(theta):
        X = _design(x, unpack(theta))
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        return y - X @ beta

    fit = least_squares(resid, theta0, bounds=(lower, upper), x_scale='jac')
    params = unpack(fit.x)
    X = _design(x, params)
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    return params, beta, float(((y - X @ beta) ** 2).sum())


def fit_edges(x, profiles, *, w_bounds=W_BOUNDS, flat=FLAT, bic_margin=BIC_MARGIN, min_plateau=MIN_PLATEAU, n_m1=41,
              n_w1=31, n_m2=15, n_w2=9):
    """One- and two-edge logistic fits of each profile on ``x`` (E1–E3).

    ``profiles`` is profiles x len(x) (no NaN). Midpoints are kept between x[1] and x[-2] (sorted x). M2 also needs
    its edges separated: >= ``min_plateau`` bins between the midpoints and a distance >= their mean width. Returns
    (models, edges): one row per profile (model 'M1' / 'M2' / 'flat', BICs, counted edges) and one row per counted
    edge (profile, position, width, amplitude, direction), edges of a profile sorted by position."""
    x = np.asarray(x, float)
    Y = np.atleast_2d(np.asarray(profiles, float))
    if Y.shape[1] != len(x) or not np.isfinite(Y).all():
        raise ValueError('Need finite profiles x len(x).')
    xs = np.sort(x)
    lo, hi = xs[1], xs[-2]
    n = len(x)
    m1 = np.linspace(lo, hi, n_m1)
    w1 = np.geomspace(*w_bounds, n_w1)
    combos1 = [(m, w) for m in m1 for w in w1]
    m2 = np.linspace(lo, hi, n_m2)
    w2 = np.geomspace(*w_bounds, n_w2)
    combos2 = [(a, wa, b, wb) for i, a in enumerate(m2) for b in m2[i + 1:] for wa in w2 for wb in w2
               if ((x > a) & (x < b)).sum() >= min_plateau and b - a >= (wa + wb) / 2]
    best1 = np.argmin(_rss_grid(x, Y, combos1), axis=0)
    best2 = np.argmin(_rss_grid(x, Y, combos2), axis=0)
    model_rows, edge_rows = [], []

    def separated(par):
        """Two edges enclose >= min_plateau bins and lie at least their mean width apart (else a spike, not a step)."""
        (a, wa), (b, wb) = sorted([(par[0], par[1]), (par[2], par[3])])
        return ((x > a) & (x < b)).sum() >= min_plateau and b - a >= (wa + wb) / 2

    for p in range(len(Y)):
        y = Y[p]
        par1, beta1, rss1 = _refine(x, y, combos1[best1[p]], lo, hi, w_bounds)
        par2, beta2, rss2 = _refine(x, y, combos2[best2[p]], lo, hi, w_bounds)
        if not separated(par2):                                    # refinement collapsed the edges: keep the grid fit
            par2 = list(combos2[best2[p]])
            X2 = _design(x, par2)
            beta2 = np.linalg.lstsq(X2, y, rcond=None)[0]
            rss2 = float(((y - X2 @ beta2) ** 2).sum())
        bic1 = n * np.log(max(rss1, 1e-300) / n) + 4 * np.log(n)
        bic2 = n * np.log(max(rss2, 1e-300) / n) + 7 * np.log(n)
        two = [(par2[0], par2[1], beta2[1]), (par2[2], par2[3], beta2[2])]
        use_two = (bic2 <= bic1 - bic_margin and all(abs(b) >= flat for _, _, b in two) and separated(par2))
        edges = two if use_two else [(par1[0], par1[1], beta1[1])]
        edges = sorted([e for e in edges if abs(e[2]) >= flat], key=lambda e: e[0])
        model = ('M2' if use_two else 'M1') if edges else 'flat'
        model_rows.append({'profile': p, 'model': model, 'bic M1': bic1, 'bic M2': bic2, 'rss M1': rss1, 'rss M2': rss2,
                           'edges': len(edges)})
        for k, (m, w, b) in enumerate(edges):
            edge_rows.append({'profile': p, 'edge': k, 'position': float(m), 'width': float(w), 'amplitude': float(b),
                              'direction': int(np.sign(b))})
    edges = pd.DataFrame(edge_rows, columns=['profile', 'edge', 'position', 'width', 'amplitude', 'direction'])
    return pd.DataFrame(model_rows), edges


def edge_class(width, *, sharp=SHARP, gradual=GRADUAL):
    """'sharp' (w <= sharp), 'gradual' (w >= gradual) or 'graded' (E4); '' for NaN."""
    if not np.isfinite(width):
        return ''
    return 'sharp' if width <= sharp else 'gradual' if width >= gradual else 'graded'


def edge_location(position, *, window=.15, segments=('S1', 'S2', 'S3')):
    """'at S1/S2' or 'at S2/S3' within ``window`` of a boundary, else 'inside S<k>' (E5; segment units)."""
    if not np.isfinite(position):
        return ''
    for b in range(1, len(segments)):
        if abs(position - b) <= window:
            return f'at {segments[b - 1]}/{segments[b]}'
    return f'inside {segments[int(np.clip(np.floor(position), 0, len(segments) - 1))]}'


def match_edges(a, b, *, max_distance=.25):
    """One-to-one pairs (i, j) of rows of edge frames ``a`` and ``b`` with the same direction and midpoints within
    ``max_distance``, smallest distance first (R1)."""
    cand = sorted((abs(ea.position - eb.position), i, j) for i, ea in a.iterrows() for j, eb in b.iterrows()
                  if ea.direction == eb.direction and abs(ea.position - eb.position) <= max_distance)
    used_a, used_b, pairs = set(), set(), []
    for _, i, j in cand:
        if i not in used_a and j not in used_b:
            pairs.append((i, j))
            used_a.add(i)
            used_b.add(j)
    return pairs


def agreed_edges(a, b, *, max_distance=.25, window=.15):
    """Edges two specimens share (R1–R2): mean position, width and amplitude, per-specimen widths, class, location.

    The class is that of the mean width, or 'mice disagree on width' when one specimen is sharp and the other
    gradual."""
    rows = []
    for i, j in match_edges(a, b, max_distance=max_distance):
        ea, eb = a.loc[i], b.loc[j]
        width = (ea.width + eb.width) / 2
        classes = {edge_class(ea.width), edge_class(eb.width)}
        position = (ea.position + eb.position) / 2
        rows.append({'position': position, 'width': width, 'amplitude': (ea.amplitude + eb.amplitude) / 2,
                     'direction': int(ea.direction), 'width A': ea.width, 'width B': eb.width,
                     'position A': ea.position, 'position B': eb.position,
                     'class': 'mice disagree on width' if classes == {'sharp', 'gradual'} else edge_class(width),
                     'location': edge_location(position, window=window)})
    cols = ['position', 'width', 'amplitude', 'direction', 'width A', 'width B', 'position A', 'position B', 'class',
            'location']
    return pd.DataFrame(rows, columns=cols).sort_values('position').reset_index(drop=True)


def crossing_width(x, y, *, lo=.1, hi=.9):
    """Distance between the first points where a curve rising from 0 to 1 reaches ``lo`` and ``hi``."""
    return _crossing(np.asarray(x, float), np.asarray(y, float), hi) - _crossing(np.asarray(x, float), np.asarray(y, float), lo)


def spline_step_response(position, step_at, knots, grid):
    """Least-squares fit of notebook 12's cubic B-spline basis (``stats_gam.make_gam_design``) to a 0/1 step at
    ``step_at`` on the units' ``position``, evaluated on ``grid`` (all on the raw coordinate) (Q1)."""
    from .stats_gam import make_gam_design

    position = np.asarray(position, float)
    X = make_gam_design(position, knots)
    beta = np.linalg.lstsq(X, (position >= step_at).astype(float), rcond=None)[0]
    return make_gam_design(np.asarray(grid, float), knots) @ beta


def edge_composition(edge, gene_edges, core, *, near=.5, share=.5, spread=.25, min_genes=3):
    """The core genes behind one pathway edge (G2).

    ``edge`` has position and direction; ``gene_edges`` is the agreed gene edges with a ``gene`` column. Genes with
    an edge of the same direction within ``near`` of the pathway edge count. Returns the counts per class, the
    interquartile range of their positions and a label for a graded or gradual edge."""
    g = gene_edges[gene_edges.gene.isin(list(core)) & (gene_edges.direction == edge['direction'])
                   & ((gene_edges.position - edge['position']).abs() <= near)]
    g = g[g['class'].isin(['sharp', 'graded', 'gradual'])].drop_duplicates('gene')
    counts = g['class'].value_counts().reindex(['sharp', 'graded', 'gradual'], fill_value=0)
    out = {'genes': len(g), **{f'genes {k}': int(v) for k, v in counts.items()},
           'position IQR': float(np.subtract(*np.percentile(g.position, [75, 25]))) if len(g) else np.nan}
    if len(g) < min_genes:
        out['composition'] = 'not evaluable'
    elif counts['sharp'] / len(g) >= share and out['position IQR'] >= spread:
        out['composition'] = 'staggered gene switches'
    elif counts['gradual'] / len(g) >= share:
        out['composition'] = 'gene gradient'
    else:
        out['composition'] = 'mixed'
    return out
