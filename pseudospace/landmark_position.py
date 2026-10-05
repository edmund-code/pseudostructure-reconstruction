"""Landmark count-position model (LCP) for repeated PT cross-sections (notebooks 51-52).

Each structure sits at an unknown grid position z on one shared, ordered object. Landmark genes,
chosen from external segment references only, carry the positional signal. Counts follow a
negative binomial with the structure's library as exposure, so cut size and read depth enter as an
offset and cannot become position. One rate path per gene is shared by every structure of a species
and estimated with equal specimen weight at each grid point: this is the repeated-object prior. An
external segment-level shape initializes and weakly shrinks each rate path; it fixes orientation and
gauge, and segment labels are never used.

Nothing here reads tissue coordinates, depth or labels. Posteriors are conditional on the model,
including the tempered independent-gene likelihood; they are not calibrated anatomical intervals.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.special import gammaln, logsumexp

from pseudospace.stats_gam import gam_internal_knots, make_gam_design


# ---------------------------------------------------------------------------------------- basics
def grid_points(size=50):
    if size < 3:
        raise ValueError('Need at least three grid points.')
    return np.linspace(0.0, 1.0, int(size))


def piecewise_shape(points, values, grid):
    """Linear interpolation of external segment levels; constant beyond the first and last point.

    ``values`` is genes x points (any finite units). Returns genes x grid, centred over the grid.
    """
    points = np.asarray(points, float)
    values = np.atleast_2d(np.asarray(values, float))
    if values.shape[1] != len(points) or np.any(np.diff(points) <= 0) or not np.isfinite(values).all():
        raise ValueError('Need increasing points and finite values with one column per point.')
    shape = np.vstack([np.interp(grid, points, row) for row in values])
    return shape - shape.mean(axis=1, keepdims=True)


def spline_basis(grid, df=6):
    """Centred cubic B-spline basis on [0, 1] with ``df`` columns (uniform knots, no intercept)."""
    grid = np.asarray(grid, float)
    if int(df) < 4 or grid.min() < 0 or grid.max() > 1:
        raise ValueError('Need df >= 4 and a grid inside [0, 1].')
    knots = np.linspace(0, 1, int(df) - 1)[1:-1]
    basis = make_gam_design(grid, knots)[:, 1:]
    return basis - basis.mean(axis=0)


def landmark_ratio_score(late_counts, early_counts):
    """log((late + 0.5) / (early + 0.5)): exposure-free in expectation; larger means later."""
    late, early = np.asarray(late_counts, float), np.asarray(early_counts, float)
    if late.shape != early.shape or (late < 0).any() or (early < 0).any():
        raise ValueError('Need nonnegative late and early counts of equal shape.')
    return np.log((late + .5) / (early + .5))


def _dense(counts):
    y = counts.toarray() if sparse.issparse(counts) else np.asarray(counts)
    y = np.asarray(y, float)
    if y.ndim != 2 or not np.isfinite(y).all() or (y < 0).any():
        raise ValueError('Counts must be a finite nonnegative 2-D array.')
    return y


def nb_grid_loglik(Y, log_library, eta, theta, *, chunk=400):
    """Negative-binomial log-likelihood of each row at every grid point, up to row constants.

    sum_g y_ig * eta_gk - (theta_g + y_ig) * log(theta_g + l_i * exp(eta_gk)); terms that do not
    depend on the grid point are dropped (they cancel in posteriors).
    """
    n = Y.shape[0]
    out = Y @ eta  # n x K
    log_theta = np.log(theta)[None, :, None]
    for start in range(0, n, chunk):
        rows = slice(start, min(n, start + chunk))
        log_mu = log_library[rows, None, None] + eta[None, :, :]          # r x G x K
        out[rows] -= np.einsum('rg,rgk->rk', theta[None, :] + Y[rows], np.logaddexp(log_theta, log_mu))
    return out


def nb_row_constant(Y, log_library, theta):
    """The grid-independent part of the NB log-pmf, so predictive likelihoods are complete."""
    t = theta[None, :]
    return (gammaln(Y + t) - gammaln(t) - gammaln(Y + 1) + t * np.log(t) + Y * log_library[:, None]).sum(axis=1)


def posterior_summaries(posterior, grid):
    grid = np.asarray(grid, float)
    cdf = np.cumsum(posterior, axis=1)
    quantile = lambda p: grid[np.minimum((cdf < p).sum(axis=1), len(grid) - 1)]
    with np.errstate(divide='ignore', invalid='ignore'):
        entropy = -np.nansum(np.where(posterior > 0, posterior * np.log(posterior), 0.0), axis=1)
    return {'mean': posterior @ grid, 'q10': quantile(.10), 'q90': quantile(.90),
            'map': grid[np.argmax(posterior, axis=1)], 'entropy': entropy}


# ---------------------------------------------------------------------------------------- model
def _balanced_aggregates(Y, library, posterior, specimen, min_mass):
    """Per-specimen posterior-weighted rates averaged with equal specimen weight at each grid point.

    Returns pseudo-counts S (G x K) and exposure E (K): S = mean_s(rate_s) * E, E = pooled exposure.
    Grid points that no specimen occupies with at least ``min_mass`` structures get zero exposure.
    """
    names = np.unique(specimen)
    G, K = Y.shape[1], posterior.shape[1]
    rates, used = np.zeros((G, K)), np.zeros(K)
    exposure = posterior.T @ library
    for name in names:
        rows = specimen == name
        q = posterior[rows]
        mass, e_s = q.sum(axis=0), q.T @ library[rows]
        ok = (mass >= min_mass) & (e_s > 0)
        r = (q.T @ Y[rows]).T  # G x K counts
        rates[:, ok] += r[:, ok] / e_s[ok]
        used += ok
    has = used > 0
    rates[:, has] /= used[has]
    exposure = np.where(has, exposure, 0.0)
    return rates * exposure[None, :], exposure


def _penalized_poisson_paths(S, E, X, theta0, lam, *, start=None, iterations=30):
    """Newton fit of eta_g = X @ beta_g to grid pseudo-counts, penalty lam/2 |beta_g[1:] - beta0_g[1:]|^2."""
    G, p = theta0.shape
    beta = (theta0 if start is None else start).copy()
    P = np.diag(np.r_[0.0, np.ones(p - 1)])

    def objective(b):
        eta = b @ X.T
        return (S * eta - E[None, :] * np.exp(np.clip(eta, -50, 30))).sum(axis=1) \
            - .5 * lam * (((b - theta0)[:, 1:]) ** 2).sum(axis=1)

    current = objective(beta)
    for _ in range(iterations):
        eta = np.clip(beta @ X.T, -50, 30)
        mu = E[None, :] * np.exp(eta)
        grad = (S - mu) @ X - lam * (beta - theta0) @ P
        hess = np.einsum('kp,gk,kq->gpq', X, mu, X) + lam * P[None] + 1e-8 * np.eye(p)[None]
        step = np.linalg.solve(hess, grad[..., None])[..., 0]
        scale = np.ones(G)
        for _ in range(20):  # per-gene step halving keeps every gene's objective nondecreasing
            trial = beta + scale[:, None] * step
            value = objective(trial)
            worse = value < current - 1e-9
            if not worse.any():
                break
            scale[worse] *= .5
        trial = beta + scale[:, None] * step
        value = objective(trial)
        better = value >= current - 1e-9
        beta[better], current[better] = trial[better], value[better]
        if np.max(np.abs(scale[:, None] * step)) < 1e-7:
            break
    return beta


def _moment_dispersion(Y, mean, lo=1e-4, hi=10.0):
    """Per-gene 1/theta from E[(y - mu)^2] = mu + mu^2 / theta, clipped to [lo, hi]."""
    inv = ((Y - mean) ** 2 - mean).sum(axis=0) / np.maximum((mean ** 2).sum(axis=0), 1e-12)
    return 1.0 / np.clip(inv, lo, hi)


def fit_lcp(counts, library, specimen, shape, *, grid_size=50, df=6, lam=10.0, temperature=1.0,
            max_iter=100, tol=1e-4, min_mass=.5, init_noise=0.0, seed=0):
    """Fit the landmark count-position model by tempered, specimen-balanced EM.

    counts: structures x landmark genes (raw counts); library: exposure per structure; specimen:
    replicate labels; shape: genes x grid centred external log-shape (``piecewise_shape``). The
    shape sets each gene's initial path and the shrinkage target (lam), hence orientation and gauge.
    Returns the rate paths, dispersions, posteriors and their summaries, plus the objective trace.
    """
    Y = _dense(counts)
    library = np.asarray(library, float).ravel()
    specimen = np.asarray(specimen).astype(str)
    grid = grid_points(grid_size)
    shape = np.asarray(shape, float)
    if len(library) != len(Y) or len(specimen) != len(Y) or (library <= 0).any():
        raise ValueError('Need one positive library and one specimen label per structure.')
    if shape.shape != (Y.shape[1], len(grid)) or not np.isfinite(shape).all():
        raise ValueError('shape must be genes x grid and finite.')
    if not (temperature >= 1.0):
        raise ValueError('temperature must be >= 1.')
    B = spline_basis(grid, df)
    X = np.column_stack([np.ones(len(grid)), B])
    b0 = np.linalg.lstsq(B, shape.T, rcond=None)[0].T  # G x df
    level = np.log(np.maximum(Y.sum(axis=0), .5) / library.sum())
    theta0 = np.column_stack([level, b0])
    beta = theta0.copy()
    if init_noise > 0:
        beta[:, 1:] += np.random.default_rng(seed).normal(0, init_noise, size=b0.shape)
    theta = np.full(Y.shape[1], 10.0)
    log_l = np.log(library)
    trace, previous = [], None
    converged = False
    for iteration in range(max_iter):
        eta = beta @ X.T
        ll = nb_grid_loglik(Y, log_l, eta, theta)
        tempered = ll / temperature
        norm = logsumexp(tempered, axis=1, keepdims=True)
        posterior = np.exp(tempered - norm)
        trace.append(float(norm.sum() + nb_row_constant(Y, log_l, theta).sum() / temperature
                           - len(Y) * np.log(len(grid))))
        mean = posterior @ grid
        if previous is not None and np.max(np.abs(mean - previous)) < tol:
            converged = True
            break
        previous = mean
        S, E = _balanced_aggregates(Y, library, posterior, specimen, min_mass)
        beta = _penalized_poisson_paths(S, E, X, theta0, lam, start=beta)
        rate_mean = library[:, None] * (posterior @ np.exp(np.clip(beta @ X.T, -50, 30)).T)
        theta = _moment_dispersion(Y, rate_mean)
    eta = beta @ X.T
    trace_arr = np.asarray(trace)
    return {'grid': grid, 'eta': eta, 'theta': theta, 'beta': beta, 'temperature': float(temperature),
            'posterior': posterior, **posterior_summaries(posterior, grid), 'objective': trace_arr,
            'iterations': iteration + 1, 'converged': converged,
            'objective_decreases': int(np.sum(np.diff(trace_arr) < -1e-6 * np.abs(trace_arr[:-1]).clip(1)))}


def project_lcp(model, counts, library, *, temperature=None):
    """Posterior positions of new structures under frozen rate paths and dispersions."""
    Y = _dense(counts)
    library = np.asarray(library, float).ravel()
    if Y.shape[1] != model['eta'].shape[0] or len(library) != len(Y) or (library <= 0).any():
        raise ValueError('Counts must use the model genes; libraries must be positive.')
    t = model['temperature'] if temperature is None else float(temperature)
    ll = nb_grid_loglik(Y, np.log(library), model['eta'], model['theta']) / t
    posterior = np.exp(ll - logsumexp(ll, axis=1, keepdims=True))
    return {'posterior': posterior, **posterior_summaries(posterior, model['grid'])}


def predictive_loglik(model, counts_b, library_b):
    """Half-B log predictive density under the half-A posterior mixture of the fitted model."""
    Y = _dense(counts_b)
    library_b = np.asarray(library_b, float).ravel()
    log_l = np.log(library_b)
    ll = nb_grid_loglik(Y, log_l, model['eta'], model['theta']) + nb_row_constant(Y, log_l, model['theta'])[:, None]
    with np.errstate(divide='ignore'):
        return float(logsumexp(np.log(model['posterior']) + ll, axis=1).sum())


# ---------------------------------------------------------------------------------------- metrics
def residual_spline_gain(train_y, train_z, train_segment, train_base, test_y, test_z, test_segment, test_base,
                         *, basis_df=4, ridge=1e-6):
    """P4c: does a within-segment spline of the coordinate improve on a baseline prediction?

    ``*_base`` are baseline predictions (for example the segment + coverage oracle). Within each
    segment a cubic spline of the coordinate is fitted to training residuals; the coordinate is
    rescaled to [0, 1] by the training range of that segment and test values are clipped to it.
    Returns per-gene MSEs of the baseline and of baseline + spline on the test rows, the training
    variance, and per-segment test row masks.
    """
    ytr, yte = (np.asarray(v.toarray() if sparse.issparse(v) else v, float) for v in (train_y, test_y))
    ztr, zte = np.asarray(train_z, float), np.asarray(test_z, float)
    seg_tr, seg_te = np.asarray(train_segment).astype(str), np.asarray(test_segment).astype(str)
    resid = ytr - np.asarray(train_base, float)
    correction = np.zeros_like(yte)
    for s in np.unique(seg_te):
        tr, te = seg_tr == s, seg_te == s
        if tr.sum() < 10 or not np.isfinite(ztr[tr]).all():
            continue
        lo, hi = ztr[tr].min(), ztr[tr].max()
        if hi <= lo:
            continue
        x_tr = (ztr[tr] - lo) / (hi - lo)
        x_te = np.clip((zte[te] - lo) / (hi - lo), 0, 1)
        knots = gam_internal_knots(x_tr, basis_df=basis_df)
        D_tr, D_te = make_gam_design(x_tr, knots), make_gam_design(x_te, knots)
        coef = np.linalg.solve(D_tr.T @ D_tr + ridge * np.eye(D_tr.shape[1]), D_tr.T @ resid[tr])
        correction[te] = D_te @ coef
    base = np.asarray(test_base, float)
    return {'base_mse': ((yte - base) ** 2), 'coordinate_mse': ((yte - base - correction) ** 2),
            'train_var': ytr.var(axis=0), 'segments': seg_te}


def scaled_gain(base_sq, new_sq, train_var, rows=None):
    """1 - sum_g mean(new)/var_g / sum_g mean(base)/var_g over the selected test rows."""
    rows = slice(None) if rows is None else rows
    var = np.maximum(np.asarray(train_var, float), 1e-12)
    return float(1 - np.sum(new_sq[rows].mean(axis=0) / var) / np.sum(base_sq[rows].mean(axis=0) / var))


def _poisson_glm(Y, D, offset, *, iterations=25, ridge=1e-3):
    """Vectorized Poisson Newton fit with a shared design (first column = intercept) and log offset.

    A small ridge on the non-intercept coefficients keeps sparse genes finite; returns G x p.
    """
    n, G = Y.shape
    p = D.shape[1]
    rate0 = np.log(np.maximum(Y.sum(axis=0), .5) / np.exp(offset).sum())
    beta = np.zeros((G, p))
    beta[:, 0] = rate0
    penalty = ridge * np.diag(np.r_[0.0, np.ones(p - 1)])
    for _ in range(iterations):
        eta = np.clip(offset[:, None] + D @ beta.T, -50, 30)
        mu = np.exp(eta)
        grad = D.T @ (Y - mu) - (beta @ penalty).T  # p x G
        hess = np.einsum('np,ng,nq->gpq', D, mu, D) + penalty[None] + 1e-9 * np.eye(p)[None]
        step = np.linalg.solve(hess, grad.T[..., None])[..., 0]
        beta += step
        if np.max(np.abs(step)) < 1e-6:
            break
    return beta


def poisson_deviance(Y, mu):
    Y = np.asarray(Y, float)
    mu = np.maximum(np.asarray(mu, float), 1e-12)
    with np.errstate(divide='ignore', invalid='ignore'):
        term = np.where(Y > 0, Y * np.log(Y / mu), 0.0)
    return 2 * (term - (Y - mu))


def deviance_gain(train_counts, train_library, train_z, train_segment, test_counts, test_library, test_z,
                  test_segment, *, basis_df=4):
    """P4c-dev: Poisson deviance of segment rates versus a within-segment spline of the coordinate.

    Both models use the library as offset and are fitted on the training specimen only. Returns the
    test deviance matrices of both models and the test segment labels.
    """
    ytr, yte = (np.asarray(v.toarray() if sparse.issparse(v) else v, float) for v in (train_counts, test_counts))
    ltr, lte = np.asarray(train_library, float), np.asarray(test_library, float)
    ztr, zte = np.asarray(train_z, float), np.asarray(test_z, float)
    seg_tr, seg_te = np.asarray(train_segment).astype(str), np.asarray(test_segment).astype(str)
    mu_a, mu_b = np.zeros_like(yte), np.zeros_like(yte)
    for s in np.unique(seg_te):
        tr, te = seg_tr == s, seg_te == s
        rate = np.maximum(ytr[tr].sum(axis=0), .5) / ltr[tr].sum()
        mu_a[te] = lte[te, None] * rate[None, :]
        lo, hi = ztr[tr].min(), ztr[tr].max()
        x_tr = (ztr[tr] - lo) / (hi - lo)
        x_te = np.clip((zte[te] - lo) / (hi - lo), 0, 1)
        knots = gam_internal_knots(x_tr, basis_df=basis_df)
        D_tr, D_te = make_gam_design(x_tr, knots), make_gam_design(x_te, knots)
        beta = _poisson_glm(ytr[tr], D_tr, np.log(ltr[tr]))
        mu_b[te] = np.exp(np.clip(np.log(lte[te])[:, None] + D_te @ beta.T, -50, 30))
        silent = ytr[tr].sum(axis=0) == 0  # no training counts: no information, keep the segment rate
        mu_b[np.ix_(te, silent)] = mu_a[np.ix_(te, silent)]
    return {'segment_deviance': poisson_deviance(yte, mu_a), 'coordinate_deviance': poisson_deviance(yte, mu_b),
            'segments': seg_te}
