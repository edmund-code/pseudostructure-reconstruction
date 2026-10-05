"""Landmark reconstruction, phase 3 (notebooks 52-53): a gauge-anchored landmark model, zone-derived
landmarks, bootstrap uncertainty for the landmark ratio, and the within-segment residual transfer test.

Notebook 51's ``pseudospace.landmark_position`` is left unchanged for its reproducibility; this
module builds on its likelihood helpers.

Model. Each structure sits at an unknown position z on one shared ordered object. Landmark counts
follow a negative binomial with the structure's library as exposure. Rate paths (spline or linear in
z) are shared by all structures of a block (species); the M-step can weight specimens equally at
every grid point. The likelihood is invariant to monotone re-parameterizations of z, so after every
M-step the grid is re-mapped to make (specimen-balanced or pooled) occupancy uniform: this anchors
the gauge and removes the flat direction that kept notebook 51's EM from converging. A coarse grid
is followed by a fine one. An optional zone prior uses coarse discrete labels as soft order.

Nothing here reads tissue coordinates. Posteriors are conditional on the model; their calibration is
checked separately (design-effect temperature, standardized differences between disjoint halves).
"""
from __future__ import annotations

from itertools import permutations

import numpy as np
from scipy import sparse
from scipy.special import logsumexp

from pseudospace.landmark_position import (_moment_dispersion, _penalized_poisson_paths, landmark_ratio_score,
                                           nb_grid_loglik, nb_row_constant, piecewise_shape, poisson_deviance)
from pseudospace.stats_gam import gam_internal_knots, make_gam_design

REFERENCE_GRID = np.linspace(0.0, 1.0, 201)


# ---------------------------------------------------------------------------------------- basis
class PathBasis:
    """Centred basis functions of z on [0, 1]: 'spline' (cubic B-spline, df columns) or 'linear'."""

    def __init__(self, kind='spline', df=6):
        if kind not in ('spline', 'linear'):
            raise ValueError("kind must be 'spline' or 'linear'")
        if kind == 'spline' and int(df) < 4:
            raise ValueError('spline df must be >= 4')
        self.kind, self.df = kind, (int(df) if kind == 'spline' else 1)
        self.knots = np.linspace(0, 1, self.df - 1)[1:-1] if kind == 'spline' else None
        self.center = self._raw(REFERENCE_GRID).mean(axis=0)

    def _raw(self, x):
        x = np.clip(np.asarray(x, float), 0, 1)
        if self.kind == 'linear':
            return x[:, None]
        return make_gam_design(x, self.knots)[:, 1:]

    def at(self, x):
        return self._raw(x) - self.center

    def design(self, x):
        return np.column_stack([np.ones(len(np.atleast_1d(x))), self.at(x)])


def zone_log_prior(labels, edges, grid, eps=0.10):
    """log p(z | zone label): (1 - eps) uniform inside the zone, eps uniform elsewhere; -1 = no label."""
    labels = np.asarray(labels, int)
    edges = np.asarray(edges, float)
    grid = np.asarray(grid, float)
    if not (0 < eps < 1) or edges[0] != 0 or edges[-1] != 1 or np.any(np.diff(edges) <= 0):
        raise ValueError('Need 0 < eps < 1 and increasing zone edges from 0 to 1.')
    zone_of_grid = np.clip(np.searchsorted(edges, grid, side='right') - 1, 0, len(edges) - 2)
    out = np.zeros((len(labels), len(grid)))
    for j in range(len(edges) - 1):
        inside = zone_of_grid == j
        n_in, n_out = inside.sum(), (~inside).sum()
        if n_in == 0:
            raise ValueError('A zone has no grid points; use a finer grid.')
        row = np.where(inside, (1 - eps) / n_in, eps / max(n_out, 1))
        out[labels == j] = np.log(row / row.sum())
    out[labels < 0] = -np.log(len(grid))
    return out


# ---------------------------------------------------------------------------------------- fitting
def _dense(counts):
    y = counts.toarray() if sparse.issparse(counts) else np.asarray(counts)
    y = np.asarray(y, float)
    if y.ndim != 2 or not np.isfinite(y).all() or (y < 0).any():
        raise ValueError('Counts must be a finite nonnegative 2-D array.')
    return y


def _aggregates(Y, library, posterior, specimen, balanced, min_mass=.5):
    if not balanced:
        return (posterior.T @ Y).T, posterior.T @ library
    G, K = Y.shape[1], posterior.shape[1]
    rates, used = np.zeros((G, K)), np.zeros(K)
    exposure = posterior.T @ library
    for name in np.unique(specimen):
        rows = specimen == name
        q = posterior[rows]
        mass, e_s = q.sum(axis=0), q.T @ library[rows]
        ok = (mass >= min_mass) & (e_s > 0)
        r = (q.T @ Y[rows]).T
        rates[:, ok] += r[:, ok] / e_s[ok]
        used += ok
    has = used > 0
    rates[:, has] /= used[has]
    exposure = np.where(has, exposure, 0.0)
    return rates * exposure[None, :], exposure


def _occupancy(posteriors, specimens, balanced):
    """Occupancy over the grid: mean of per-specimen normalized mass (balanced) or pooled mass."""
    parts = []
    for q, spec in zip(posteriors, specimens):
        if balanced:
            parts += [q[spec == s].sum(axis=0) / (spec == s).sum() for s in np.unique(spec)]
        else:
            parts.append(q.sum(axis=0))
    occ = np.mean(parts, axis=0) if balanced else np.sum(parts, axis=0)
    occ = occ / occ.sum() + 1e-9
    return occ / occ.sum()


def _anchor_map(occ, grid):
    """Old grid positions -> new gauge with uniform occupancy (endpoints fixed)."""
    cdf = np.cumsum(occ) - occ / 2
    cdf = (cdf - cdf[0]) / (cdf[-1] - cdf[0])
    return np.maximum.accumulate(cdf)


def _warp_beta(beta, basis, grid, mapping):
    """Re-express paths in the new gauge: eta_new(z') = eta_old(u) with mapping(u) = z'."""
    u = np.interp(grid, mapping, grid)  # old location of each new grid point
    eta_new = beta @ basis.design(u).T
    X = basis.design(grid)
    return np.linalg.lstsq(X, eta_new.T, rcond=None)[0].T


def _shared_paths(S0, E0, S1, E1, X, start0, start1, b0, lam, iterations=30):
    """Joint penalized Poisson fit: shared shape b, block-specific levels a0, a1 (genes x ...)."""
    G, p = start0.shape
    B = X[:, 1:]
    theta = np.column_stack([start0[:, 0], start1[:, 0], 0.5 * (start0[:, 1:] + start1[:, 1:])])

    def objective(t):
        e0 = np.clip(t[:, [0]] + t[:, 2:] @ B.T, -50, 30)
        e1 = np.clip(t[:, [1]] + t[:, 2:] @ B.T, -50, 30)
        return ((S0 * e0 - E0[None] * np.exp(e0)).sum(1) + (S1 * e1 - E1[None] * np.exp(e1)).sum(1)
                - .5 * lam * ((t[:, 2:] - b0) ** 2).sum(1))

    current = objective(theta)
    D0 = np.column_stack([np.ones(len(X)), np.zeros(len(X)), B])
    D1 = np.column_stack([np.zeros(len(X)), np.ones(len(X)), B])
    P = np.diag(np.r_[0., 0., np.ones(p - 1)])
    for _ in range(iterations):
        e0 = np.clip(theta @ D0.T, -50, 30)
        e1 = np.clip(theta @ D1.T, -50, 30)
        mu0, mu1 = E0[None] * np.exp(e0), E1[None] * np.exp(e1)
        pen = np.zeros_like(theta)
        pen[:, 2:] = theta[:, 2:] - b0
        grad = (S0 - mu0) @ D0 + (S1 - mu1) @ D1 - lam * pen
        hess = (np.einsum('kp,gk,kq->gpq', D0, mu0, D0) + np.einsum('kp,gk,kq->gpq', D1, mu1, D1)
                + lam * P[None] + 1e-8 * np.eye(p + 1)[None])
        step = np.linalg.solve(hess, grad[..., None])[..., 0]
        scale = np.ones(G)
        for _ in range(20):
            value = objective(theta + scale[:, None] * step)
            worse = value < current - 1e-9
            if not worse.any():
                break
            scale[worse] *= .5
        trial = theta + scale[:, None] * step
        value = objective(trial)
        better = value >= current - 1e-9
        theta[better], current[better] = trial[better], value[better]
        if np.max(np.abs(scale[:, None] * step)) < 1e-7:
            break
    return np.column_stack([theta[:, 0], theta[:, 2:]]), np.column_stack([theta[:, 1], theta[:, 2:]])


def posterior_summary(posterior, grid):
    grid = np.asarray(grid, float)
    mean = posterior @ grid
    sd = np.sqrt(np.maximum(posterior @ grid ** 2 - mean ** 2, 0.0))
    cdf = np.cumsum(posterior, axis=1)
    quantile = lambda p: grid[np.minimum((cdf < p).sum(axis=1), len(grid) - 1)]
    with np.errstate(divide='ignore', invalid='ignore'):
        entropy = -np.sum(np.where(posterior > 0, posterior * np.log(posterior), 0.0), axis=1)
    return {'mean': mean, 'sd': sd, 'q10': quantile(.1), 'q90': quantile(.9),
            'map': grid[np.argmax(posterior, axis=1)], 'entropy': entropy}


def fit_landmark_model(blocks, *, basis='spline', df=6, lam=10.0, balanced=True, anchor=True, shared=None,
                       coarse=50, fine=200, coarse_max=300, fine_max=500, tol_mean=1e-4, tol_rel=1e-9,
                       patience=5, min_mass=.5, init_noise=0.0, seed=0):
    """Fit one or two blocks (species) by tempered EM with gauge anchoring, coarse then fine grid.

    Each block is a dict with counts (structures x genes), library, specimen, shape_points,
    shape_values (genes x points, log scale), temperature, and optionally labels + zone_edges + eps
    for a zone prior. ``shared`` lists (gene in block 0, gene in block 1) pairs whose shape
    coefficients are shared (species-specific levels); with sharing, one joint occupancy anchors
    both blocks. Returns one result dict per block plus fit diagnostics.
    """
    pb = PathBasis(basis, df)
    data = []
    for blk in blocks:
        Y = _dense(blk['counts'])
        lib = np.asarray(blk['library'], float).ravel()
        spec = np.asarray(blk['specimen']).astype(str)
        if len(lib) != len(Y) or len(spec) != len(Y) or (lib <= 0).any():
            raise ValueError('Need one positive library and one specimen per structure.')
        T = float(blk.get('temperature', 1.0))
        if T < 1:
            raise ValueError('temperature must be >= 1')
        shape = piecewise_shape(blk['shape_points'], blk['shape_values'], REFERENCE_GRID)
        if shape.shape[0] != Y.shape[1]:
            raise ValueError('shape_values must have one row per gene.')
        b0 = np.linalg.lstsq(pb.at(REFERENCE_GRID), shape.T, rcond=None)[0].T
        level = np.log(np.maximum(Y.sum(axis=0), .5) / lib.sum())
        theta0 = np.column_stack([level, b0])
        beta = theta0.copy()
        if init_noise > 0:
            beta[:, 1:] += np.random.default_rng(seed).normal(0, init_noise, size=b0.shape)
        data.append({'Y': Y, 'lib': lib, 'log_l': np.log(lib), 'spec': spec, 'T': T, 'theta0': theta0,
                     'beta': beta, 'disp': np.full(Y.shape[1], 10.0),
                     'labels': blk.get('labels'), 'edges': blk.get('zone_edges'), 'eps': blk.get('eps', .10),
                     'const_rows': None})
    if shared is not None:
        if len(blocks) != 2:
            raise ValueError('Sharing needs exactly two blocks.')
        shared = np.asarray(shared, int).reshape(-1, 2)
    trace, stage_iters = [], {}
    converged = False
    for stage, K, max_iter in (('coarse', coarse, coarse_max), ('fine', fine, fine_max)):
        grid = np.linspace(0, 1, K)
        X = pb.design(grid)
        for d in data:
            d['log_prior'] = (zone_log_prior(d['labels'], d['edges'], grid, d['eps']) if d['labels'] is not None
                              else np.full((len(d['Y']), K), -np.log(K)))
        previous, stable, last_obj = None, 0, None
        for it in range(max_iter):
            objective = 0.0
            for d in data:
                eta = d['beta'] @ X.T
                ll = nb_grid_loglik(d['Y'], d['log_l'], eta, d['disp']) / d['T'] + d['log_prior']
                norm = logsumexp(ll, axis=1, keepdims=True)
                d['post'] = np.exp(ll - norm)
                objective += float(norm.sum() + nb_row_constant(d['Y'], d['log_l'], d['disp']).sum() / d['T'])
            trace.append(objective)
            means = np.concatenate([d['post'] @ grid for d in data])
            if previous is not None:
                rel = abs(objective - last_obj) / max(abs(last_obj), 1.0)
                small = np.max(np.abs(means - previous)) < tol_mean
                if stage == 'coarse' and small:
                    break
                stable = stable + 1 if (small and rel < tol_rel) else 0
                if stage == 'fine' and stable >= patience:
                    converged = True
                    break
            previous, last_obj = means, objective
            # M-step
            aggs = [_aggregates(d['Y'], d['lib'], d['post'], d['spec'], balanced, min_mass) for d in data]
            for d, (S, E) in zip(data, aggs):
                d['beta'] = _penalized_poisson_paths(S, E, X, d['theta0'], lam, start=d['beta'])
            if shared is not None:
                (S0, E0), (S1, E1) = aggs
                g0, g1 = shared[:, 0], shared[:, 1]
                b0 = 0.5 * (data[0]['theta0'][g0, 1:] + data[1]['theta0'][g1, 1:])
                n0, n1 = _shared_paths(S0[g0], E0, S1[g1], E1, X, data[0]['beta'][g0], data[1]['beta'][g1], b0, lam)
                data[0]['beta'][g0], data[1]['beta'][g1] = n0, n1
            for d in data:
                rate = d['lib'][:, None] * (d['post'] @ np.exp(np.clip(d['beta'] @ X.T, -50, 30)).T)
                d['disp'] = _moment_dispersion(d['Y'], rate)
            if anchor:
                if shared is not None:
                    maps = [_anchor_map(_occupancy([d['post'] for d in data], [d['spec'] for d in data], balanced), grid)] * 2
                else:
                    maps = [_anchor_map(_occupancy([d['post']], [d['spec']], balanced), grid) for d in data]
                for d, mapping in zip(data, maps):
                    d['beta'] = _warp_beta(d['beta'], pb, grid, mapping)
        stage_iters[stage] = it + 1
    out = []
    for d in data:
        out.append({'grid': grid, 'eta': d['beta'] @ X.T, 'beta': d['beta'], 'theta': d['disp'], 'temperature': d['T'],
                    'posterior': d['post'], **posterior_summary(d['post'], grid)})
    diag = {'objective': np.asarray(trace), 'iterations_coarse': stage_iters['coarse'],
            'iterations_fine': stage_iters['fine'], 'converged': converged}
    return out, diag


def project_landmark_model(result, counts, library, *, temperature=None, log_prior=None):
    """Posterior of new structures under frozen paths and dispersions (fine grid)."""
    Y = _dense(counts)
    library = np.asarray(library, float).ravel()
    T = result['temperature'] if temperature is None else float(temperature)
    grid = result['grid']
    ll = nb_grid_loglik(Y, np.log(library), result['eta'], result['theta']) / T
    ll = ll + (log_prior if log_prior is not None else -np.log(len(grid)))
    post = np.exp(ll - logsumexp(ll, axis=1, keepdims=True))
    return {'posterior': post, **posterior_summary(post, grid)}


def combine_half_posteriors(post_a, post_b, log_prior=None):
    """Normalized product of two independent half-likelihood posteriors (each includes the prior once)."""
    with np.errstate(divide='ignore'):
        log = np.log(np.maximum(post_a, 1e-300)) + np.log(np.maximum(post_b, 1e-300))
    if log_prior is not None:
        log = log - log_prior
    return np.exp(log - logsumexp(log, axis=1, keepdims=True))


def design_effect_temperature(mean_1, sd_1, mean_2, sd_2):
    """T = max(1, (Q0.8(|d|) / 1.2816)^2) with d the standardized difference of two gene-half fits."""
    d = (np.asarray(mean_1) - np.asarray(mean_2)) / np.sqrt(np.asarray(sd_1) ** 2 + np.asarray(sd_2) ** 2 + 1e-12)
    return float(max(1.0, (np.quantile(np.abs(d), .8) / 1.2815516) ** 2)), d


def standardized_difference(mean_1, sd_1, mean_2, sd_2):
    return (np.asarray(mean_1) - np.asarray(mean_2)) / np.sqrt(np.asarray(sd_1) ** 2 + np.asarray(sd_2) ** 2 + 1e-12)


# ---------------------------------------------------------------------------------------- LRS-U
def lrs_bootstrap_sd(late_counts, early_counts, *, n_boot=200, seed=0):
    """SD of the landmark ratio under Poisson count resampling plus landmark-gene resampling."""
    L, E = _dense(late_counts), _dense(early_counts)
    rng = np.random.default_rng(seed)
    draws = np.empty((n_boot, len(L)))
    for b in range(n_boot):
        gl = rng.integers(0, L.shape[1], L.shape[1])
        ge = rng.integers(0, E.shape[1], E.shape[1])
        late = rng.poisson(L[:, gl]).sum(axis=1)
        early = rng.poisson(E[:, ge]).sum(axis=1)
        draws[b] = landmark_ratio_score(late, early)
    return draws.std(axis=0, ddof=1)


# ---------------------------------------------------------------------------------------- zones
def zone_order(profiles, start):
    """Shortest Hamiltonian path through zone profiles (zones x features) beginning at ``start``."""
    P = np.asarray(profiles, float)
    n = len(P)
    if n > 8:
        raise ValueError('Brute force is limited to 8 zones.')
    dist = np.sqrt(((P[:, None, :] - P[None, :, :]) ** 2).sum(-1))
    best, best_len = None, np.inf
    for rest in permutations([i for i in range(n) if i != start]):
        path = (start, *rest)
        length = sum(dist[a, b] for a, b in zip(path[:-1], path[1:]))
        if length < best_len - 1e-12:
            best, best_len = path, length
    return list(best), float(best_len)


def zone_rates(counts, library, zone, specimen, n_zones):
    """Specimen-balanced zone rates (genes x zones): mean over specimens of sum(y) / sum(library)."""
    Y = counts.tocsr() if sparse.issparse(counts) else sparse.csr_matrix(np.asarray(counts))
    out = np.zeros((Y.shape[1], n_zones))
    for j in range(n_zones):
        per = []
        for s in np.unique(specimen):
            rows = (zone == j) & (specimen == s)
            if rows.sum() == 0:
                continue
            per.append(np.asarray(Y[rows].sum(axis=0)).ravel() / library[rows].sum())
        out[:, j] = np.mean(per, axis=0)
    return out


def zone_f_statistic(counts, library, zone, specimen, n_zones):
    """Quasi-Poisson F of zone-specific rates against one constant rate (library offset).

    Fitted means use specimen-balanced zone rates; the dispersion is the Pearson statistic of the
    zone model per gene. Returns F (genes) and the zone rates (genes x zones).
    """
    Y = counts.tocsr() if sparse.issparse(counts) else sparse.csr_matrix(np.asarray(counts))
    rates = zone_rates(Y, library, zone, specimen, n_zones)
    const = np.asarray(Y.sum(axis=0)).ravel() / library.sum()
    dense = Y.toarray().astype(float)
    mu_z = library[:, None] * rates[:, zone].T
    mu_c = library[:, None] * const[None, :]
    with np.errstate(divide='ignore', invalid='ignore'):
        dev = lambda mu: 2 * np.sum(np.where(dense > 0, dense * np.log(dense / np.maximum(mu, 1e-12)), 0.0)
                                    - (dense - mu), axis=0)
        phi = np.sum((dense - mu_z) ** 2 / np.maximum(mu_z, 1e-12), axis=0) / max(len(dense) - n_zones, 1)
        F = (dev(mu_c) - dev(mu_z)) / (n_zones - 1) / np.maximum(phi, 1e-12)
    return np.nan_to_num(F, nan=0.0), rates


# ---------------------------------------------------------------------------------------- tests
def _segment_spline(z_tr, z_te, y_tr, df=4, ridge=1e-6):
    lo, hi = np.min(z_tr), np.max(z_tr)
    if not hi > lo:
        return np.zeros((len(z_te), y_tr.shape[1]))
    x_tr = (z_tr - lo) / (hi - lo)
    x_te = np.clip((z_te - lo) / (hi - lo), 0, 1)
    knots = gam_internal_knots(x_tr, basis_df=df)
    D_tr, D_te = make_gam_design(x_tr, knots), make_gam_design(x_te, knots)
    coef = np.linalg.solve(D_tr.T @ D_tr + ridge * np.eye(D_tr.shape[1]), D_tr.T @ y_tr)
    return D_te @ coef


def wsr_correlations(train_resid, train_z, train_segment, test_resid, test_z, test_segment, perms, *, df=4):
    """Within-segment residual transfer: per segment, per-gene Pearson r between the training spline's
    prediction at test positions and the observed test residual, plus the same under permutations.

    ``perms`` maps segment -> (n_perm x n_test_rows_in_segment) index arrays applied to the test
    coordinate. Returns segment -> (r_observed (genes), r_null (n_perm x genes)).
    """
    Rtr, Rte = np.asarray(train_resid, float), np.asarray(test_resid, float)
    ztr, zte = np.asarray(train_z, float), np.asarray(test_z, float)
    seg_tr, seg_te = np.asarray(train_segment).astype(str), np.asarray(test_segment).astype(str)
    out = {}
    for s, perm in perms.items():
        tr, te = seg_tr == s, seg_te == s
        y = Rte[te]
        yc = y - y.mean(axis=0)
        ynorm = np.sqrt((yc ** 2).sum(axis=0))

        def corr(pred):
            pc = pred - pred.mean(axis=0)
            with np.errstate(invalid='ignore', divide='ignore'):
                return (pc * yc).sum(axis=0) / (np.sqrt((pc ** 2).sum(axis=0)) * ynorm)

        pred = _segment_spline(ztr[tr], zte[te], Rtr[tr], df)
        observed = corr(pred)
        # The prediction at a permuted coordinate is the prediction matrix with permuted rows.
        null = np.vstack([corr(pred[p]) for p in perm])
        out[s] = (observed, null)
    return out


def wsr_statistic(r):
    r = np.clip(np.asarray(r, float), -.999999, .999999)
    return float(np.nanmean(np.arctanh(r), axis=-1)) if r.ndim == 1 else np.nanmean(np.arctanh(r), axis=-1)


def _poisson_glm_ridge(Y, D, offset, *, ridge=1.0, iterations=50):
    n, G = Y.shape
    p = D.shape[1]
    beta = np.zeros((G, p))
    beta[:, 0] = np.log(np.maximum(Y.sum(axis=0), .5) / np.exp(offset).sum())
    penalty = ridge * np.diag(np.r_[0.0, np.ones(p - 1)])
    for _ in range(iterations):
        mu = np.exp(np.clip(offset[:, None] + D @ beta.T, -50, 30))
        grad = D.T @ (Y - mu) - (beta @ penalty).T
        hess = np.einsum('np,ng,nq->gpq', D, mu, D) + penalty[None] + 1e-9 * np.eye(p)[None]
        step = np.linalg.solve(hess, grad.T[..., None])[..., 0]
        beta += step
        if np.max(np.abs(step)) < 1e-6:
            break
    return beta


def deviance_gain_repaired(train_counts, train_library, train_z, train_segment, test_counts, test_library,
                           test_z, test_segment, *, basis_df=4, min_train_counts=10, ridge=1.0):
    """P4c-dev, repaired: genes with < ``min_train_counts`` training counts in a segment are excluded
    from both models there; spline coefficients carry a ridge. Returns both deviance matrices (excluded
    entries are zero) and the test segment labels."""
    ytr, yte = (np.asarray(v.toarray() if sparse.issparse(v) else v, float) for v in (train_counts, test_counts))
    ltr, lte = np.asarray(train_library, float), np.asarray(test_library, float)
    ztr, zte = np.asarray(train_z, float), np.asarray(test_z, float)
    seg_tr, seg_te = np.asarray(train_segment).astype(str), np.asarray(test_segment).astype(str)
    dev_a, dev_b = np.zeros_like(yte), np.zeros_like(yte)
    for s in np.unique(seg_te):
        tr, te = seg_tr == s, seg_te == s
        keep = ytr[tr].sum(axis=0) >= min_train_counts
        if not keep.any() or tr.sum() < 10:
            continue
        rate = ytr[tr][:, keep].sum(axis=0) / ltr[tr].sum()
        mu_a = lte[te, None] * rate[None, :]
        lo, hi = ztr[tr].min(), ztr[tr].max()
        x_tr = (ztr[tr] - lo) / (hi - lo)
        x_te = np.clip((zte[te] - lo) / (hi - lo), 0, 1)
        knots = gam_internal_knots(x_tr, basis_df=basis_df)
        beta = _poisson_glm_ridge(ytr[tr][:, keep], make_gam_design(x_tr, knots), np.log(ltr[tr]), ridge=ridge)
        mu_b = np.exp(np.clip(np.log(lte[te])[:, None] + make_gam_design(x_te, knots) @ beta.T, -50, 30))
        cols = np.flatnonzero(keep)
        dev_a[np.ix_(te, cols)] = poisson_deviance(yte[np.ix_(te, cols)], mu_a)
        dev_b[np.ix_(te, cols)] = poisson_deviance(yte[np.ix_(te, cols)], mu_b)
    return {'segment_deviance': dev_a, 'coordinate_deviance': dev_b, 'segments': seg_te}
