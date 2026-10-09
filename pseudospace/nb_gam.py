"""Negative-binomial GLM counterparts of the gene-level Gaussian GAMs, on raw counts.

Log link, offset log(library), and the same designs and prior weights as the Gaussian screens
(``pathway_calibration.contrast_designs``, ``specimen_null.segment_designs``,
``pathway_remodeling.fit_nested_trajectories``): w_i = n / (2 J_s n_j), each group half the weight and each
specimen an equal share. The weights enter the log-likelihood and the deviance,
D = 2 sum_i w_i [y log(y/mu) - (y + 1/a) log((1 + a y)/(1 + a mu))].
They balance specimens and are not precision weights, so the default uncertainty is HC3.

The per-gene dispersion a is estimated once, by maximum likelihood under a full design (no Cox-Reid adjustment),
and passed explicitly to every nested fit. The deviance difference LR replaces the partial F as the ranking
statistic. LR is chi-square only with unit weights (equal specimen sizes) and a known a; otherwise use it to rank.
Each design is fitted on an orthonormal basis of its column space, so redundant columns are allowed. The fit is a
vectorised IRLS with step-halving and R glm's deviance convergence criterion.
A gene without any count is reported as not converged: it carries no information. ``separated`` flags genes
whose design loses rank on the structures with a positive count, such as a gene with no counts in one species.
Their ML fit lies on the boundary (some fitted rates go to zero), as in R glm. The limiting deviances still give a
valid LR, so separated genes are neither failed nor filled. For a species without counts, the spatial LR is near 0
(no positional information) and the level LR is large. Coefficients, curves and SEs on the empty side are extreme,
or NaN where the sandwich is undefined, so mask z-scores by ``separated``.
"""
from __future__ import annotations

import numpy as np
from joblib import Parallel, delayed, parallel_config
from scipy import sparse
from scipy.special import digamma, gammaln, polygamma, xlogy
from threadpoolctl import threadpool_limits

from .pathway_calibration import COMPARISONS, contrast_designs
from .stats_gam import gam_internal_knots, make_gam_design

GENE_MODEL = 'nb.1'
A_MIN, A_MAX = 1e-8, 1e4                    # dispersion bounds; A_MIN is Poisson to working precision
MU_MIN, MU_MAX = 1e-12, 1e8                 # fitted means are kept in this range
# Ponytail: per-gene ML dispersion without shrinkage. DESeq2-style trend shrinkage left the SCF13 joint-test calls
# unchanged (same 31; gene ranks rho 0.999), so it is not implemented; add it here if rare genes start to matter.


def _orthobasis(x):
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    return u[:, s > np.finfo(float).eps * max(x.shape) * s[0]]


def _nested(inner, outer):
    return np.linalg.norm(inner - outer @ (outer.T @ inner)) <= 1e-8 * max(1., np.sqrt(inner.shape[1]))


def _outer(q):
    """Row outer products q_i q_i' as (n, r * r): Q' diag(v) Q for many v is then one matrix product."""
    n, r = q.shape
    return (q[:, :, None] * q[:, None, :]).reshape(n, r * r)


def _mean(eta):
    with np.errstate(over='ignore'):                # exp overflow is clipped to MU_MAX
        return np.clip(np.exp(eta), MU_MIN, MU_MAX)


def _deviance(y, mu, alpha, w):
    """Weighted NB deviance per gene (column); alpha: (G,)."""
    term = xlogy(y, y / mu) - (y + 1. / alpha) * (np.log1p(alpha * y) - np.log1p(alpha * mu))
    return 2. * (w[:, None] * term).sum(axis=0)


def _loglik(y, mu, alpha, w):
    theta = 1. / alpha
    ll = (gammaln(y + theta) - gammaln(theta) - gammaln(y + 1) + theta * np.log(theta / (theta + mu))
          + xlogy(y, mu / (theta + mu)))
    return (w[:, None] * ll).sum(axis=0)


def _solve(qq, q, wts, rhs):
    """Batched weighted least squares: for each gene g solve (Q' W_g Q) b = Q' rhs_g. Returns (G, r)."""
    r = q.shape[1]
    a = (qq.T @ wts).T.reshape(-1, r, r)
    b = (q.T @ rhs).T
    try:
        return np.linalg.solve(a, b[:, :, None])[:, :, 0]
    except np.linalg.LinAlgError:           # an exactly singular gene: the others are solved as in the batch
        out = np.empty_like(b)
        for g in range(len(b)):
            try:
                out[g] = np.linalg.solve(a[g], b[g])
            except np.linalg.LinAlgError:
                out[g] = np.linalg.lstsq(a[g], b[g], rcond=None)[0]
        return out


def _irls(y, q, qq, offset, w, alpha, mu_start, *, maxit=100, tol=1e-8, max_halving=30):
    """NB GLM fits of all columns of y on basis q at fixed alpha. Returns dict(beta, mu, dev, converged).

    mu_start should lie in the model space (a nested model's fit, or a constant rate), so the starting beta is exact.
    Convergence: |dev - dev_old| / (|dev| + 0.1) < tol (R glm's criterion), with step-halving on deviance increase.
    A gene without any count never counts as converged.
    """
    mu = np.clip(mu_start, MU_MIN, MU_MAX)
    work = w[:, None] * mu / (1. + alpha * mu)
    beta = _solve(qq, q, work, work * (np.log(mu) - offset[:, None]))
    mu = _mean(q @ beta.T + offset[:, None])
    dev = _deviance(y, mu, alpha, w)
    converged = np.zeros(y.shape[1], bool)
    active = np.arange(y.shape[1])
    for _ in range(maxit):
        ya, mua, aa, ba, da = y[:, active], mu[:, active], alpha[active], beta[active], dev[active]
        wts = w[:, None] * mua / (1. + aa * mua)
        z = np.log(mua) - offset[:, None] + (ya - mua) / mua
        new = _solve(qq, q, wts, wts * z)
        mu_new = _mean(q @ new.T + offset[:, None])
        d_new = _deviance(ya, mu_new, aa, w)
        bad = ~np.isfinite(d_new) | (d_new > da + 1e-10 * (np.abs(da) + 1.))
        for _ in range(max_halving):
            if not bad.any():
                break
            new[bad] = .5 * (new[bad] + ba[bad])
            mu_new[:, bad] = _mean(q @ new[bad].T + offset[:, None])
            d_new[bad] = _deviance(ya[:, bad], mu_new[:, bad], aa[bad], w)
            bad = ~np.isfinite(d_new) | (d_new > da + 1e-10 * (np.abs(da) + 1.))
        done = np.abs(d_new - da) / (np.abs(d_new) + .1) < tol
        beta[active], mu[:, active], dev[active] = new, mu_new, d_new
        converged[active[done]] = True
        active = active[~done]
        if not len(active):
            break
    return {'beta': beta, 'mu': mu, 'dev': dev, 'converged': converged & np.isfinite(dev) & (y > 0).any(axis=0)}


def _separated(y, qq, r):
    """The design loses rank on the rows with a positive count (smallest eigenvalue of that Gram block)."""
    gram = ((y > 0).T.astype(float) @ qq).reshape(-1, r, r)
    return np.linalg.eigvalsh(gram)[:, 0] <= 1e-10


def _ml_alpha(y, mu, w, alpha0, *, maxit=60, tol=1e-7, gtol=1e-6):
    """Per-gene ML of the NB dispersion given mu (weighted log-likelihood), safeguarded Newton on phi = log(1/alpha).

    A gene stops when its step is below ``tol``, its gradient in phi is below ``gtol`` (log-likelihood units), or it
    sits on a bound with the gradient pointing outward. Converged genes are frozen.
    """
    phi = -np.log(np.clip(alpha0, A_MIN, A_MAX))
    lo, hi = -np.log(A_MAX), -np.log(A_MIN)
    converged = np.zeros(len(phi), bool)
    active = np.arange(len(phi))
    for _ in range(maxit):
        ya, ma, ph = y[:, active], mu[:, active], phi[active]
        th = np.exp(ph)[None, :]
        tm = th + ma
        score = (w[:, None] * (digamma(ya + th) - digamma(th) + np.log(th) + 1. - np.log(tm)
                               - (ya + th) / tm)).sum(0)
        hess = (w[:, None] * (polygamma(1, ya + th) - polygamma(1, th) + 1. / th - 2. / tm
                              + (ya + th) / tm ** 2)).sum(0)
        t = th[0]
        g1 = t * score
        g2 = t * t * hess + t * score
        step = np.where(g2 < 0, -g1 / np.where(g2 < 0, g2, -1.), np.sign(g1))
        new = np.clip(ph + np.clip(step, -2., 2.), lo, hi)
        ll_old = _loglik(ya, ma, np.exp(-ph), w)
        ll_new = _loglik(ya, ma, np.exp(-new), w)
        for _ in range(20):
            worse = ll_new < ll_old - 1e-10 * (np.abs(ll_old) + 1.)
            if not worse.any():
                break
            new[worse] = .5 * (new[worse] + ph[worse])
            ll_new[worse] = _loglik(ya[:, worse], ma[:, worse], np.exp(-new[worse]), w)
        done = (np.abs(new - ph) < tol) | (np.abs(g1) < gtol) | ((ph <= lo) & (g1 < 0)) | ((ph >= hi) & (g1 > 0))
        phi[active] = new
        converged[active[done]] = True
        active = active[~done]
        if not len(active):
            break
    return np.exp(-phi), converged


def _constant_rate(y, offset, w):
    """Weighted constant-rate fit: the start of a design without a nested, already-fitted design."""
    rate = (w[:, None] * y).sum(axis=0) / (w * np.exp(offset)).sum()
    return np.exp(offset)[:, None] * rate[None, :]


def _fit_dispersion(y, q, qq, offset, w, *, rounds=40, tol=1e-6):
    """Joint ML of (beta, alpha) under basis q: near-Poisson start, moment alpha, then alternate IRLS(beta | alpha)
    and ML(alpha | mu) per gene until |change in log alpha| < tol (genes that settle are frozen)."""
    fit = _irls(y, q, qq, offset, w, np.full(y.shape[1], 1e-4), _constant_rate(y, offset, w))
    mu = fit['mu']
    num = (w[:, None] * ((y - mu) ** 2 - mu)).sum(0)
    alpha = np.clip(num / (w[:, None] * mu ** 2).sum(0), 1e-4, A_MAX)
    ml_ok = np.zeros(y.shape[1], bool)
    stable = np.zeros(y.shape[1], bool)
    active = np.arange(y.shape[1])
    for _ in range(rounds):
        f = _irls(y[:, active], q, qq, offset, w, alpha[active], mu[:, active])
        mu[:, active] = f['mu']
        new, ok = _ml_alpha(y[:, active], f['mu'], w, alpha[active])
        settled = np.abs(np.log(new) - np.log(alpha[active])) < tol
        alpha[active], ml_ok[active] = new, ok
        stable[active[settled]] = True
        active = active[~settled]
        if not len(active):
            break
    return alpha, _irls(y, q, qq, offset, w, alpha, mu), ml_ok & stable


def _contrast_se(y, mu, q, qq, w, alpha, lmat, *, hc3=True):
    """Sandwich and model-based SEs of lmat @ b per gene (b on basis q). Returns two (G, m) arrays.

    Bread A = Q' diag(w mu/(1+a mu)) Q, the expected information. The HC3 meat sums
    [w (y-mu)/(1+a mu)/(1-h)]^2 q q' with h = w mu/(1+a mu) q'A^-1 q; hc3=False gives HC0 (tests only).
    A gene with any h >= .99 gets a NaN sandwich SE (the Gaussian code refuses such designs).
    Model-based: A^-1 (sum w^2 mu/(1+a mu) q q') A^-1, again because the weights are not precision weights.
    """
    r = q.shape[1]
    info = w[:, None] * mu / (1. + alpha * mu)
    inv = np.linalg.pinv((qq.T @ info).T.reshape(-1, r, r), hermitian=True)
    score = w[:, None] * (y - mu) / (1. + alpha * mu)
    unstable = np.zeros(y.shape[1], bool)
    if hc3:
        leverage = info * (qq @ inv.reshape(-1, r * r).T)
        unstable = (leverage >= .99).any(axis=0)
        score = score / (1. - np.minimum(leverage, .99))
    proj = lmat @ inv                       # (G, m, r)
    out = []
    for middle in (score ** 2, w[:, None] * info):
        meat = (qq.T @ middle).T.reshape(-1, r, r)
        out.append(np.sqrt(np.maximum(np.einsum('gmr,grs,gms->gm', proj, meat, proj), 1e-12)))
    out[0][unstable] = np.nan
    return out


def _fit_block(y, alpha, start, plan, offset, w, contrasts, pearson):
    """Every design of ``plan`` in rank order, warm-started from its parent's fit; then contrasts and residuals."""
    out, fits, bases = {}, {}, {}
    for name, q, back, parent in plan:
        qq = _outer(q)
        if parent is not None:
            mu0 = fits[parent]['mu']
        elif start is not None:
            mu0 = _mean(q @ start.T + offset[:, None])
        else:
            mu0 = _constant_rate(y, offset, w)
        fit = fits[name] = _irls(y, q, qq, offset, w, alpha, mu0)
        bases[name] = q
        out[f'{name}|dev'], out[f'{name}|beta'] = fit['dev'], fit['beta'] @ back.T
        out[f'{name}|ok'], out[f'{name}|separated'] = fit['converged'], _separated(y, qq, q.shape[1])
    for key, name, lmat in contrasts:
        out[f'{key}|value'] = fits[name]['beta'] @ lmat.T
        q = bases[name]
        se = _contrast_se(y, fits[name]['mu'], q, _outer(q), w, alpha, lmat)
        out[f'{key}|se_hc3'], out[f'{key}|se_model'] = se
    if pearson is not None:
        mu = fits[pearson]['mu']
        out['pearson'] = (y - mu) / np.sqrt(mu + alpha * mu ** 2)
    return out


def _dispersion_block(y, q, offset, w):
    qq = _outer(q)
    alpha, fit, ok = _fit_dispersion(y, q, qq, offset, w)
    return {'alpha': alpha, 'converged': ok & fit['converged'], 'separated': _separated(y, qq, q.shape[1]),
            'deviance': fit['dev']}


def _single_thread(worker, *args):
    with threadpool_limits(1):
        return worker(*args)


def _blocks(worker, y, per_gene, shared, block_size, n_jobs):
    """Run ``worker`` on dense gene blocks (loky, one BLAS thread each); outputs concatenated in gene order."""
    if int(block_size) < 1:
        raise ValueError('block_size must be positive.')
    spans = [(a, min(a + block_size, y.shape[1])) for a in range(0, y.shape[1], block_size)]
    tasks = (delayed(_single_thread)(worker, y[:, a:b].toarray() if sparse.issparse(y) else np.array(y[:, a:b]),
                                     *[None if v is None else v[a:b] for v in per_gene], *shared)
             for a, b in spans)
    with parallel_config(backend='loky', inner_max_num_threads=1):
        parts = Parallel(n_jobs=n_jobs)(tasks)
    return {key: np.concatenate([part[key] for part in parts], axis=1 if key == 'pearson' else 0)
            for key in parts[0]}


def _inputs(counts, weights, offset):
    y = sparse.csc_matrix(counts, dtype=float) if sparse.issparse(counts) else np.asarray(counts, float)
    values = y.data if sparse.issparse(y) else y
    if y.ndim != 2 or not y.shape[1] or not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('Counts must be a finite, nonnegative structures x genes matrix.')
    w, off = np.asarray(weights, float), np.asarray(offset, float)
    if (w.shape != (y.shape[0],) or off.shape != w.shape or not np.isfinite(w).all() or (w <= 0).any() or
            not np.isfinite(off).all()):
        raise ValueError('Need positive finite weights and finite offsets aligned to rows.')
    return y, w, off


def _log_library(library, n):
    lib = np.asarray(library, float)
    if lib.shape != (n,) or not np.isfinite(lib).all() or (lib <= 0).any():
        raise ValueError('Library sizes must be positive, finite and aligned to rows.')
    return np.log(lib)


def _dispersion_of(dispersion, n_genes):
    """(alpha, converged) from nb_dispersion's dict, or a plain per-gene array (taken as converged)."""
    if isinstance(dispersion, dict):
        alpha, ok = dispersion['alpha'], dispersion['converged']
    else:
        alpha, ok = dispersion, True
    alpha = np.asarray(alpha, float)
    if alpha.shape != (n_genes,) or not np.isfinite(alpha).all() or (alpha <= 0).any():
        raise ValueError('Need one positive finite dispersion per gene.')
    return alpha, np.broadcast_to(np.asarray(ok, bool), alpha.shape).copy()


def _plan(designs, n):
    """(name, basis, pinv(X) @ basis, parent) in rank order; parent: the largest earlier design nested in it."""
    matrices, bases = {}, {}
    for name, x in designs.items():
        matrices[name] = x = np.asarray(x, float)
        if x.ndim != 2 or x.shape[0] != n or not x.shape[1] or not np.isfinite(x).all():
            raise ValueError(f'{name}: design must be a finite matrix aligned to rows.')
        bases[name] = _orthobasis(x)
    order = sorted(bases, key=lambda k: (bases[k].shape[1], k))
    plan = []
    for i, name in enumerate(order):
        inner = [m for m in order[:i] if _nested(bases[m], bases[name])]
        parent = max(inner, key=lambda m: bases[m].shape[1], default=None)
        plan.append((name, bases[name], np.linalg.pinv(matrices[name]) @ bases[name], parent))
    return plan


def _fill(values, ok, fill):
    """Failed genes get the converged genes' median (fill='median') or NaN (fill=None)."""
    out = np.array(values, float)
    out[~ok] = np.median(out[ok]) if fill == 'median' and ok.any() else np.nan
    return out


def _check_fill(fill):
    if fill not in ('median', None):
        raise ValueError("fill must be 'median' or None.")


def nb_dispersion(counts, design, weights, offset, *, block_size=200, n_jobs=1):
    """Per-gene ML dispersion under one (full) design, to be held fixed for every nested model.

    Alternates IRLS for the coefficients with safeguarded Newton for a (from a near-Poisson fit and a moment
    start), bounded to [A_MIN, A_MAX]. Returns alpha, converged (alternation settled, ML converged, final fit
    converged), separated (under ``design``), at_lower (a at A_MIN: no detectable overdispersion), deviance under
    ``design``, and gene_model.
    """
    y, w, off = _inputs(counts, weights, offset)
    (_, q, _, _), = _plan({'design': design}, y.shape[0])
    out = _blocks(_dispersion_block, y, [], (q, off, w), block_size, n_jobs)
    out['at_lower'] = out['alpha'] <= A_MIN * (1. + 1e-6)
    out['gene_model'] = GENE_MODEL
    return out


def nested_nb_lr(counts, designs, weights, offset, dispersion, comparisons=None, *, block_size=200, n_jobs=1,
                 quasi=False, fill='median'):
    """NB likelihood-ratio statistics for nested designs at fixed dispersion; drop-in for ``nested_partial_f``.

    Each comparison is (reduced, full) design names; comparisons whose designs are absent are skipped. Returns
    {name: LR} with LR = max(D_reduced - D_full, 0), with ``quasi`` also name + '_ql' = (LR/df1) / (D_full/df2),
    'converged' (every used fit converged, every LR finite, and the dispersion converged when ``dispersion`` is
    nb_dispersion's dict), 'separated' (in any used design; never filled) and 'df' = {name: (added, residual)} as
    in ``nested_partial_f``. Statistics of genes that did not converge are filled by ``fill``.
    """
    pairs = {k: v for k, v in (comparisons or COMPARISONS).items() if v[0] in designs and v[1] in designs}
    _check_fill(fill)
    if not pairs:
        raise ValueError('Need comparisons whose designs are supplied.')
    y, w, off = _inputs(counts, weights, offset)
    alpha, ok = _dispersion_of(dispersion, y.shape[1])
    n = y.shape[0]
    plan = _plan({name: designs[name] for name in sorted({x for pair in pairs.values() for x in pair})}, n)
    bases = {name: q for name, q, _, _ in plan}
    df = {}
    for name, (reduced, full) in pairs.items():
        if not _nested(bases[reduced], bases[full]):
            raise ValueError(f'{name}: {reduced} is not nested in {full}.')
        df[name] = (bases[full].shape[1] - bases[reduced].shape[1], n - bases[full].shape[1])
        if df[name][0] <= 0 or df[name][1] <= 2:
            raise ValueError(f'{name}: no added columns or no residual degrees of freedom.')
    fits = _blocks(_fit_block, y, [alpha, None], (plan, off, w, (), None), block_size, n_jobs)
    separated = np.zeros(len(ok), bool)
    for name in bases:
        ok &= fits[f'{name}|ok']
        separated |= fits[f'{name}|separated']
    scores = {}
    for name, (reduced, full) in pairs.items():
        scores[name] = np.maximum(fits[f'{reduced}|dev'] - fits[f'{full}|dev'], 0.)
        ok &= np.isfinite(scores[name])
        if quasi:
            noise = np.maximum(fits[f'{full}|dev'] / df[name][1], 1e-12)
            scores[name + '_ql'] = scores[name] / df[name][0] / noise
    scores = {name: _fill(values, ok, fill) for name, values in scores.items()}
    scores['converged'], scores['separated'], scores['df'] = ok, separated, df
    return scores


def nb_fit(counts, design, weights, offset, dispersion, *, contrasts=None, pearson=False, start=None,
           block_size=200, n_jobs=1):
    """NB GLM of every gene on one design at fixed dispersion.

    ``contrasts``: {name: C}, C (m, p) on the design's columns (estimable ones when columns are redundant).
    ``start``: optional (G, p) coefficients to start from, such as a previous fit's beta. Returns beta (G, p) in
    the design's columns (minimum norm when columns are redundant), deviance, converged, separated, contrasts
    {name: {value, se_hc3, se_model}} each (G, m), and with ``pearson`` the residuals (y - mu)/sqrt(mu + a mu^2),
    structures x genes.
    """
    y, w, off = _inputs(counts, weights, offset)
    alpha, ok = _dispersion_of(dispersion, y.shape[1])
    plan = _plan({'design': design}, y.shape[0])
    _, q, back, _ = plan[0]
    x = np.asarray(design, float)
    if start is not None:
        start = np.asarray(start, float)
        if start.shape != (y.shape[1], x.shape[1]) or not np.isfinite(start).all():
            raise ValueError('start must be finite coefficients, genes x design columns.')
        start = start @ (q.T @ x).T                     # the same linear predictor on the orthonormal basis
    lmats = {}
    for name, c in (contrasts or {}).items():
        c = np.atleast_2d(np.asarray(c, float))
        if c.ndim != 2 or c.shape[1] != x.shape[1]:
            raise ValueError(f'{name}: a contrast needs one column per design column.')
        lmats[name] = c @ back
    fits = _blocks(_fit_block, y, [alpha, start],
                   (plan, off, w, [(k, 'design', v) for k, v in lmats.items()], 'design' if pearson else None),
                   block_size, n_jobs)
    out = {'beta': fits['design|beta'], 'deviance': fits['design|dev'], 'converged': ok & fits['design|ok'],
           'separated': fits['design|separated'],
           'contrasts': {k: {s: fits[f'{k}|{s}'] for s in ('value', 'se_hc3', 'se_model')} for k in lmats}}
    if pearson:
        out['pearson'] = fits['pearson']
    return out


def nb_species_trajectories(counts, library, position, human, specimen, grid, dispersion, *, basis_df=6,
                            return_residuals=False, block_size=200, n_jobs=1, fill='median'):
    """NB counterpart of ``pathway_remodeling.fit_nested_trajectories``, with the same designs and keys.

    T_level/T_spatial/T_total are LRs (filled by ``fill`` for genes that did not converge). beta_human_level is the
    human-minus-mouse log rate ratio in M_level with its HC3 SE. mouse/human/delta are natural-log rates on
    ``grid`` (add log 1e4 for log CP10k) at the equal-specimen species mean, with HC3 se/z and model-based
    se_model. residuals are full-model Pearson residuals. Also returns converged, separated and gene_model.
    Estimate ``dispersion`` with nb_dispersion under contrast_designs(position, human, specimen,
    gam_internal_knots(position, basis_df))[0]['full'], its weights and log(library).
    """
    _check_fill(fill)
    s, c = np.asarray(position, float), np.asarray(human, float)
    samples, grid = np.asarray(specimen, str), np.asarray(grid, float)
    n = len(s)
    if counts.shape[0] != n or s.shape != c.shape or s.shape != samples.shape:
        raise ValueError('Count rows, coordinate, species and specimens must align.')
    if (not np.isfinite(s).all() or not np.isfinite(grid).all() or s.min() < 0 or s.max() > 1 or grid.size < 2 or
            grid.min() < s.min() or grid.max() > s.max() or np.any(np.diff(grid) <= 0)):
        raise ValueError('Need an increasing grid within the observed [0, 1] coordinate support.')
    if set(np.unique(c)) != {0., 1.} or basis_df < 3:
        raise ValueError('Need both species (mouse=0, human=1) and basis_df >= 3.')
    knots = gam_internal_knots(s, basis_df=basis_df)
    designs, weights = contrast_designs(s, c, samples, knots)
    for x in designs.values():
        if np.linalg.matrix_rank(x) != x.shape[1] or n <= x.shape[1] + 2:
            raise ValueError('Insufficient coordinate support for the requested spline/specimen design.')
    y, w, off = _inputs(counts, weights, _log_library(library, n))
    alpha, ok = _dispersion_of(dispersion, y.shape[1])
    base = make_gam_design(grid, knots)
    p_base = base.shape[1]
    mouse_design = np.zeros((len(grid), designs['full'].shape[1]))
    mouse_design[:, :p_base] = base
    human_design = mouse_design.copy()
    human_design[:, p_base] = 1
    human_design[:, p_base + 1:2 * p_base] = base[:, 1:]
    plan = _plan(designs, n)
    back = {name: b for name, _, b, _ in plan}
    contrasts = [('level', 'level', np.eye(designs['level'].shape[1])[[p_base]] @ back['level']),
                 ('delta', 'full', (human_design - mouse_design) @ back['full'])]
    fits = _blocks(_fit_block, y, [alpha, None],
                   (plan, off, w, contrasts, 'full' if return_residuals else None), block_size, n_jobs)
    result = {}
    separated = np.zeros(len(ok), bool)
    for name in ('base', 'level', 'full'):
        ok &= fits[f'{name}|ok']
        separated |= fits[f'{name}|separated']
    for key in ('T_level', 'T_spatial', 'T_total'):
        reduced, full = COMPARISONS[key]
        result[key] = np.maximum(fits[f'{reduced}|dev'] - fits[f'{full}|dev'], 0.)
        ok &= np.isfinite(result[key])
    result = {key: _fill(values, ok, fill) for key, values in result.items()}
    level, level_se = fits['level|value'][:, 0], fits['level|se_hc3'][:, 0]
    delta, se = fits['delta|value'], fits['delta|se_hc3']
    result.update({'beta_human_level': level, 'se_human_level': level_se, 'Z_level': level / level_se,
                   'mouse': fits['full|beta'] @ mouse_design.T, 'human': fits['full|beta'] @ human_design.T,
                   'delta': delta, 'se': se, 'z': delta / se, 'se_model': fits['delta|se_model'],
                   'knots': knots, 'grid': grid,
                   'design_df': np.array([designs[k].shape[1] for k in ('base', 'level', 'full')]),
                   'converged': ok, 'separated': separated, 'gene_model': GENE_MODEL})
    if return_residuals:
        result['residuals'] = fits['pearson']
    return result


def nb_group_difference_curves(counts, library, position, group, specimen, knots, grid, dispersion, *,
                               nuisance_shape=None, block_size=200, n_jobs=1, fill='median'):
    """NB counterpart of ``continuous_descriptors.group_difference_curves``.

    Group-1 minus group-0 log rate ratio Δ(s) on ``grid`` from contrast_designs' full model (genes x grid), its HC3
    SE, the LR T_spatial of level versus full (filled by ``fill``), converged and separated.
    """
    _check_fill(fill)
    designs, weights = contrast_designs(position, group, specimen, knots, nuisance_shape=nuisance_shape)
    designs = {name: designs[name] for name in ('level', 'full')}
    n = len(weights)
    if counts.shape[0] != n:
        raise ValueError('Count rows and positions must align.')
    y, w, off = _inputs(counts, weights, _log_library(library, n))
    alpha, ok = _dispersion_of(dispersion, y.shape[1])
    base = make_gam_design(np.asarray(grid, float), knots)
    p_b = base.shape[1]
    rows = np.zeros((len(base), designs['full'].shape[1]))
    rows[:, p_b], rows[:, p_b + 1:2 * p_b] = 1, base[:, 1:]
    plan = _plan(designs, n)
    back = {name: b for name, _, b, _ in plan}
    fits = _blocks(_fit_block, y, [alpha, None], (plan, off, w, [('delta', 'full', rows @ back['full'])], None),
                   block_size, n_jobs)
    lr = np.maximum(fits['level|dev'] - fits['full|dev'], 0.)
    ok &= fits['level|ok'] & fits['full|ok'] & np.isfinite(lr)
    return {'delta': fits['delta|value'], 'se': fits['delta|se_hc3'], 'T_spatial': _fill(lr, ok, fill),
            'converged': ok, 'separated': fits['level|separated'] | fits['full|separated']}
