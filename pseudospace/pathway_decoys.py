"""Decoy calibration, the joint matched-correlation test and descriptive grouping of PT pathways.

With two specimens per species, the only exchangeable specimen relabelings are the two balanced
mouse+human groupings of notebook 31. Under independent specimen deviations each one is a draw of
the section-level noise that the species contrast also carries (notebook 37 derives this). Pooling
those decoys over pathways estimates a false-discovery proportion for a whole screen (target-decoy,
Elias & Gygi 2007). It is not a per-pathway p-value, and it cannot see noise shared by both human
sections (one donor).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


def joint_matched_test(tests, vif, *, keys=('partition', 'statistic')):
    """Covariate-matched z deflated by a per-pathway variance inflation; one-sided p and BH q.

    ``tests`` holds matched rank-AUC rows (auc, null_auc_mean, null_auc_sd, pathway_id). ``vif``
    maps pathway_id to a variance inflation (CAMERA-style, >= 0). The matched null keeps the
    measurability matching; the inflation adds inter-gene correlation, which the matched null
    ignores. A missing inflation gives p = 1 (a reserved BH slot), never an untested pass.
    """
    need = {'pathway_id', 'auc', 'null_auc_mean', 'null_auc_sd'}
    if not need.issubset(tests):
        raise ValueError(f'Missing columns: {sorted(need - set(tests))}')
    out = tests.copy()
    inflation = out.pathway_id.map(pd.Series(vif, dtype=float))
    if (inflation <= 0).any():
        raise ValueError('Variance inflation must be positive.')
    out['z_matched'] = (out.auc - out.null_auc_mean) / out.null_auc_sd.clip(lower=1e-12)
    out['vif'] = inflation
    out['z_joint'] = out.z_matched / np.sqrt(inflation)
    out['p_joint'] = norm.sf(out.z_joint).astype(float)
    out.loc[inflation.isna(), 'p_joint'] = 1.
    present = [k for k in keys if k in out]
    if present:
        out['q_joint'] = out.groupby(present).p_joint.transform(lambda p: multipletests(p, method='fdr_bh')[1])
    else:
        out['q_joint'] = multipletests(out.p_joint, method='fdr_bh')[1]
    out.loc[inflation.isna(), 'q_joint'] = np.nan
    return out


def decoy_fdp(target, decoys):
    """Target-decoy false-discovery estimate at every threshold that the target scores define.

    ``target``: Series of pathway scores (larger = stronger). ``decoys``: DataFrame on the same
    index, one column per decoy. FDP(t) = mean decoy count >= t / target count >= t, capped at 1.
    ``q`` is the smallest FDP at or below each threshold (monotone, like a q-value).
    """
    target = pd.Series(target, dtype=float)
    decoys = pd.DataFrame(decoys).reindex(target.index)
    if target.isna().any() or decoys.isna().any().any() or decoys.shape[1] < 1:
        raise ValueError('Target and decoy scores must be finite on the same pathways.')
    thresholds = np.unique(target.to_numpy())[::-1]
    sorted_target = np.sort(target.to_numpy())
    sorted_decoys = [np.sort(decoys[c].to_numpy()) for c in decoys]
    n = len(target)
    target_calls = n - np.searchsorted(sorted_target, thresholds, side='left')
    decoy_calls = np.mean([n - np.searchsorted(d, thresholds, side='left') for d in sorted_decoys], axis=0)
    fdp = np.minimum(decoy_calls / target_calls, 1.)
    q = np.minimum.accumulate(fdp[::-1])[::-1]
    return pd.DataFrame({'threshold': thresholds, 'target_calls': target_calls,
                         'decoy_calls': decoy_calls, 'fdp': fdp, 'q': q})


def decoy_calls(target, decoys, alpha=.1):
    """Pathways called at decoy-estimated FDP <= alpha (empty when none qualify)."""
    curve = decoy_fdp(target, decoys)
    passing = curve[curve.q <= alpha]
    if passing.empty:
        return pd.Index([], dtype=object)
    target = pd.Series(target, dtype=float)
    return target.index[target >= passing.threshold.min()]


def bootstrap_decoy_calls(target, decoys, alpha=.1, *, n_boot=500, seed=0):
    """Percentile interval of decoy-calibrated call counts, resampling pathways with replacement."""
    target = pd.Series(target, dtype=float)
    decoys = pd.DataFrame(decoys).reindex(target.index)
    rng = np.random.default_rng(seed)
    counts = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(target), len(target))
        t = pd.Series(target.to_numpy()[idx])
        d = pd.DataFrame(decoys.to_numpy()[idx], columns=decoys.columns)
        counts.append(len(decoy_calls(t, d, alpha)))
    return np.percentile(counts, [2.5, 50, 97.5])


def decoy_inflation(null_z, set_size):
    """Least-squares fit of E[z^2] = a * (1 + (m - 1) * rho) to decoy z of every pathway.

    ``null_z``: DataFrame (pathways x decoys); ``set_size``: Series of member counts. The form is
    the variance of a sum of m equicorrelated scores (CAMERA). Returns (a, rho, vif Series >= 1).
    """
    z = pd.DataFrame(null_z)
    m = pd.Series(set_size, dtype=float).reindex(z.index)
    if m.isna().any() or z.isna().any().any():
        raise ValueError('Every pathway needs a size and finite decoy z.')
    x = np.column_stack([np.ones(len(m)), m.to_numpy() - 1])
    design = np.vstack([x] * z.shape[1])
    response = (z.to_numpy() ** 2).T.ravel()
    intercept, slope = np.linalg.lstsq(design, response, rcond=None)[0]
    a, rho = intercept, slope / intercept if intercept > 0 else np.nan
    return a, rho, pd.Series(np.maximum(x @ [intercept, slope], 1.), index=z.index)


def overlap_programs(gene_sets, ids, *, contributing=None, cut=.5, min_contributing=3):
    """Average-linkage groups of pathways by the overlap coefficient of their contributing genes.

    Contributing genes are members in ``contributing`` (for example the top decile of the gene
    statistic); a pathway with fewer than ``min_contributing`` falls back to all members. Distance
    is 1 - |A & B| / min(|A|, |B|). Programs are numbered by size. Descriptive grouping after the
    list is frozen; no test uses it.
    """
    ids = list(ids)
    if not ids:
        return pd.DataFrame(columns=['pathway_id', 'program'])
    keep = set(contributing) if contributing is not None else None
    genes = {}
    for p in ids:
        members = set(gene_sets[p])
        chosen = members & keep if keep is not None else members
        genes[p] = chosen if len(chosen) >= min_contributing else members
    if len(ids) == 1:
        return pd.DataFrame({'pathway_id': ids, 'program': [1]})
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    d = np.zeros((len(ids), len(ids)))
    for i, a in enumerate(ids):
        for j in range(i):
            b = ids[j]
            d[i, j] = d[j, i] = 1 - len(genes[a] & genes[b]) / min(len(genes[a]), len(genes[b]))
    raw = fcluster(linkage(squareform(d), 'average'), cut, criterion='distance')
    order = pd.Series(raw).value_counts().index  # largest first; ties keep first appearance
    relabel = {old: new for new, old in enumerate(order, start=1)}
    return pd.DataFrame({'pathway_id': ids, 'program': [relabel[r] for r in raw]})


def averaging_loss(delta):
    """Share of a curve's absolute difference that cancels in its average: 1 - |mean| / mean |.|.

    0 = one-signed constant difference; 1 = fully cancelling. Rows are genes (or pathways).
    """
    delta = np.atleast_2d(np.asarray(delta, float))
    magnitude = np.abs(delta).mean(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        loss = 1 - np.abs(delta.mean(axis=1)) / magnitude
    loss[magnitude == 0] = np.nan
    return loss


def _longest_run(mask):
    best = current = 0
    for value in mask:
        current = current + 1 if value else 0
        best = max(best, current)
    return best


def shape_class(curve, *, relative=.25, run=3, local_width=1 / 3):
    """Display rule for a positional difference curve: 'reversing', 'localized' or 'graded'.

    Reversing: runs of at least ``run`` grid points on both signs with |curve| >= relative * max.
    Localized: otherwise, |curve| >= half its maximum over less than ``local_width`` of the grid.
    Graded: the rest. A display rule on a fitted curve, not a test.
    """
    curve = np.asarray(curve, float)
    peak = np.abs(curve).max()
    if not np.isfinite(peak) or peak == 0:
        return 'flat'
    strong = np.abs(curve) >= relative * peak
    if _longest_run(strong & (curve > 0)) >= run and _longest_run(strong & (curve < 0)) >= run:
        return 'reversing'
    return 'localized' if np.mean(np.abs(curve) >= .5 * peak) < local_width else 'graded'


def welch_t(a, b):
    """Column-wise Welch t of group a minus group b (rows are replicates). Returns (t, df)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.shape[0] < 2 or b.shape[0] < 2 or a.shape[1:] != b.shape[1:]:
        raise ValueError('Need at least two replicates per group and matching columns.')
    va, vb = a.var(axis=0, ddof=1) / a.shape[0], b.var(axis=0, ddof=1) / b.shape[0]
    se = np.sqrt(va + vb)
    with np.errstate(invalid='ignore', divide='ignore'):
        t = (a.mean(axis=0) - b.mean(axis=0)) / se
        df = (va + vb) ** 2 / (va ** 2 / (a.shape[0] - 1) + vb ** 2 / (b.shape[0] - 1))
    t[se == 0] = 0.
    return t, df
