"""Exploratory atlas helpers for repeated cross-sectional structures.

These small numerical routines implement a specimen-cross-fitted spline atlas;
they are research scaffolding, not a validated trajectory method.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import BSpline
from scipy.stats import rankdata, spearmanr


def _matrix(X, name="X"):
    x = np.asarray(X, dtype=float)
    if x.ndim != 2 or not x.shape[0] or not x.shape[1]:
        raise ValueError(f"{name} must be a non-empty 2D matrix")
    if not np.isfinite(x).all():
        raise ValueError(f"{name} must contain only finite values")
    return x


def _vectors(z, specimens, n, min_specimens=2):
    z = np.asarray(z, dtype=float).reshape(-1)
    s = np.asarray(specimens).reshape(-1)
    if len(z) != n or len(s) != n or not np.isfinite(z).all():
        raise ValueError("z and specimens must have one finite value per structure")
    if np.any((z < 0) | (z > 1)):
        raise ValueError("z must lie in [0, 1]")
    if len(np.unique(s)) < min_specimens:
        raise ValueError(f"at least {min_specimens} specimen(s) are required")
    return z, s


def _anatomy(anatomy, n):
    if anatomy is None:
        return None
    raw = np.asarray(anatomy).reshape(-1)
    if len(raw) != n or not np.isfinite(raw.astype(float)).all() or not np.isin(raw.astype(float), [0, 1, 2]).all():
        raise ValueError("anatomy must have one integer label in {0, 1, 2} per structure")
    return raw.astype(int)


def _specimen_weights(specimens, specimen_weights=None):
    labels = np.unique(specimens)
    if specimen_weights is None:
        return {label: 1.0 for label in labels}
    weights = {}
    for label in labels:
        if label not in specimen_weights:
            raise ValueError(f"specimen_weights is missing {label!r}")
        value = float(specimen_weights[label])
        if not np.isfinite(value) or value <= 0:
            raise ValueError("specimen weights must be positive and finite")
        weights[label] = value
    return weights


def _basis(z, grid):
    # Cubic B-splines with equally spaced interior knots and clamped endpoints.
    interior = np.linspace(0, 1, max(0, len(grid) // 8))[1:-1]
    knots = np.r_[np.zeros(4), interior, np.ones(4)]
    return np.asarray(BSpline.design_matrix(np.clip(z, 0, 1), knots, 3).toarray())


@dataclass
class SplineAtlas:
    grid: np.ndarray
    mean: np.ndarray
    scale: np.ndarray
    beta: np.ndarray  # genes x basis coefficients

    def predict(self, z):
        zz = np.asarray(z, dtype=float).reshape(-1)
        return _basis(zz, self.grid) @ self.beta.T * self.scale + self.mean


def fit_atlas(X, z, specimens, grid=None, ridge=1.0, specimen_weights=None):
    """Fit gene curves with equal specimen weights, or supplied bootstrap weights."""
    x = _matrix(X)
    zz, ss = _vectors(z, specimens, x.shape[0], min_specimens=1)
    if not np.isfinite(ridge) or ridge < 0:
        raise ValueError("ridge must be finite and non-negative")
    grid = np.linspace(0, 1, 101) if grid is None else np.asarray(grid, dtype=float)
    if grid.ndim != 1 or len(grid) < 8 or not np.isfinite(grid).all() or np.any(np.diff(grid) <= 0):
        raise ValueError("grid must be a finite, strictly increasing vector with at least 8 points")
    if grid[0] != 0 or grid[-1] != 1:
        raise ValueError("grid must span [0, 1]")
    # Scaling and curve fitting share the original-specimen weighting scheme.
    sw = _specimen_weights(ss, specimen_weights)
    weights = np.zeros(len(ss), dtype=float)
    for label in np.unique(ss):
        mask = ss == label
        weights[mask] = sw[label] / mask.sum()
    weights /= sum(sw.values())
    mean = np.sum(x * weights[:, None], axis=0)
    scale = np.sqrt(np.sum((x - mean) ** 2 * weights[:, None], axis=0))
    scale[scale < 1e-8] = 1.0
    y = (x - mean) / scale
    B = _basis(zz, grid)
    # Normalize the data term to the number of original specimens; penalize
    # wiggliness while leaving constant curves in the penalty null space.
    weights *= len(np.unique(ss))
    rootw = np.sqrt(weights)
    lhs = (B * rootw[:, None]).T @ (B * rootw[:, None])
    d2 = np.diff(np.eye(lhs.shape[0]), n=2, axis=0)
    penalty = float(ridge) * (d2.T @ d2) + np.eye(lhs.shape[0]) * 1e-8
    beta = np.linalg.solve(lhs + penalty, (B * rootw[:, None]).T @ (y * rootw[:, None])).T
    return SplineAtlas(grid=grid, mean=mean, scale=scale, beta=beta)


def map_atlas(X, atlas, anatomy=None, initial=None, radius=None):
    """Map profiles to the best fitting atlas grid coordinate, with optional bounds."""
    x = _matrix(X)
    if x.shape[1] != len(atlas.mean):
        raise ValueError("X gene count does not match atlas")
    a = _anatomy(anatomy, x.shape[0])
    init = None if initial is None else np.asarray(initial, dtype=float).reshape(-1)
    if init is not None and (len(init) != len(x) or not np.isfinite(init).all() or np.any((init < 0) | (init > 1))):
        raise ValueError("initial must have one coordinate in [0, 1] per structure")
    if radius is not None and (not np.isfinite(radius) or radius < 0 or init is None):
        raise ValueError("radius requires initial coordinates and must be non-negative")
    standardized = (x - atlas.mean) / atlas.scale
    curves = (_basis(atlas.grid, atlas.grid) @ atlas.beta.T)
    out = np.empty(len(x))
    for i in range(len(x)):
        allowed = np.ones(len(atlas.grid), dtype=bool)
        if a is not None:
            lo, hi = a[i] / 3, (a[i] + 1) / 3
            allowed &= (atlas.grid >= lo) & ((atlas.grid <= hi) if a[i] == 2 else (atlas.grid < hi))
            if not allowed.any():
                raise ValueError(f"atlas grid has no points in anatomy interval {int(a[i])}")
        if init is not None and radius is not None:
            allowed &= np.abs(atlas.grid - init[i]) <= radius
        if not allowed.any() and a is not None:
            anatomy_allowed = (atlas.grid >= a[i] / 3) & ((atlas.grid <= (a[i] + 1) / 3) if a[i] == 2 else (atlas.grid < (a[i] + 1) / 3))
            choices = np.flatnonzero(anatomy_allowed)
            selected = choices[np.argmin(np.abs(atlas.grid[choices] - init[i]))] if init is not None else choices[0]
            allowed[selected] = True
        elif not allowed.any():
            allowed[np.argmin(np.abs(atlas.grid - init[i]))] = True
        error = np.mean((curves[allowed] - standardized[i]) ** 2, axis=1)
        out[i] = atlas.grid[allowed][np.argmin(error)]
    return out


def select_panel(X, z, specimens, anatomy, k, prior=None):
    """Rank genes by replicate curve agreement plus within-anatomy positional signal."""
    x = _matrix(X)
    zz, ss = _vectors(z, specimens, len(x))
    aa = _anatomy(anatomy, len(x))
    if aa is None or not isinstance(k, (int, np.integer)) or not 1 <= k <= x.shape[1]:
        raise ValueError("anatomy is required and k must be between 1 and n_genes")
    # Score smooth within-segment structure separately within each specimen.
    # This captures peaks as well as monotone trends and avoids between-segment
    # differences masquerading as fine positional information.
    position_parts = [[] for _ in range(x.shape[1])]
    grids = np.linspace(0, 1, 41)
    curves = []
    for specimen in np.unique(ss):
        mask = ss == specimen
        specimen_curves = _basis(grids, grids)
        for seg in (0, 1, 2):
            use = mask & (aa == seg)
            if use.sum() < 8:
                continue
            B = _basis(zz[use], grids)
            d2 = np.diff(np.eye(B.shape[1]), n=2, axis=0)
            beta = np.linalg.solve(B.T @ B + d2.T @ d2 + np.eye(B.shape[1]) * 1e-6, B.T @ x[use])
            fitted = B @ beta
            total = np.sum((x[use] - x[use].mean(axis=0)) ** 2, axis=0)
            explained = np.maximum(0, 1 - np.sum((x[use] - fitted) ** 2, axis=0) / np.maximum(total, 1e-12))
            for g, val in enumerate(explained):
                position_parts[g].append(float(val))
        B = _basis(zz[mask], grids)
        d2 = np.diff(np.eye(B.shape[1]), n=2, axis=0)
        beta = np.linalg.solve(B.T @ B + d2.T @ d2 + np.eye(B.shape[1]) * 1e-6, B.T @ x[mask])
        curves.append((specimen_curves @ beta))
    position = np.array([np.mean(vals) if vals else 0.0 for vals in position_parts])
    agreement = np.zeros(x.shape[1])
    for g in range(x.shape[1]):
        vals = []
        for i in range(len(curves)):
            for j in range(i):
                zi, zj = zz[ss == np.unique(ss)[i]], zz[ss == np.unique(ss)[j]]
                overlap = (grids >= max(zi.min(), zj.min())) & (grids <= min(zi.max(), zj.max()))
                if overlap.sum() < 8:
                    continue
                r = spearmanr(curves[i][overlap, g], curves[j][overlap, g]).statistic
                vals.append(max(0.0, float(r)) if np.isfinite(r) else 0.0)
        agreement[g] = np.mean(vals) if vals else 0.0
    def rank01(v):
        return (rankdata(v, method="average") - 1) / max(len(v) - 1, 1)
    score = 0.75 * rank01(agreement) + 0.25 * rank01(position)
    if prior is not None:
        p = np.asarray(prior, dtype=float).reshape(-1)
        if len(p) != len(score) or not np.isfinite(p).all():
            raise ValueError("prior must be a finite score per gene")
        score += 0.1 * rank01(p)
    inds = np.argsort(-score, kind="mergesort")[:k]
    return inds, score


def refine_coordinate(X, z, specimens, anatomy, iterations=3, radius=0.08, specimen_weights=None):
    """Update each specimen using atlases fitted to all other specimens.

    Later iterations are not fully independent cross-fitting: prior updates of
    a held specimen can enter the atlas used for another specimen. Use a true
    outer specimen holdout for unbiased generalization estimates.
    """
    x = _matrix(X)
    zz, ss = _vectors(z, specimens, len(x))
    aa = _anatomy(anatomy, len(x))
    if aa is None or not np.isfinite(iterations) or iterations < 0 or not np.isfinite(radius) or radius < 0:
        raise ValueError("anatomy is required; iterations and radius must be non-negative")
    _specimen_weights(ss, specimen_weights)  # validate all original identities once
    history = [zz.copy()]
    for _ in range(int(iterations)):
        updated = zz.copy()
        for specimen in np.unique(ss):
            hold = ss == specimen
            train = ~hold
            if len(np.unique(ss[train])) < 1:
                continue
            atlas = fit_atlas(x[train], zz[train], ss[train], specimen_weights=specimen_weights)
            updated[hold] = map_atlas(x[hold], atlas, aa[hold], zz[hold], radius)
        zz = updated
        history.append(zz.copy())
    return zz, history


def evaluate_prediction(X_train, z_train, s_train, X_test, z_test, gene_names=None):
    """Return per-gene held-out MSE and improvement over training-mean prediction."""
    train, test = _matrix(X_train, "X_train"), _matrix(X_test, "X_test")
    if train.shape[1] != test.shape[1]:
        raise ValueError("training and test matrices must have equal gene count")
    atlas = fit_atlas(train, z_train, s_train)
    pred = atlas.predict(z_test)
    train_specimens = np.asarray(s_train).reshape(-1)
    if len(train_specimens) != len(train):
        raise ValueError("s_train length mismatch")
    train_weights = np.zeros(len(train))
    for label in np.unique(train_specimens):
        mask = train_specimens == label
        train_weights[mask] = 1 / (mask.sum() * len(np.unique(train_specimens)))
    baseline_mean = np.sum(train * train_weights[:, None], axis=0)
    # The prediction API has no test specimen labels, so report structure-weighted
    # metrics; training centering remains specimen-balanced.
    mse = np.mean((test - pred) ** 2, axis=0)
    base = np.mean((test - baseline_mean) ** 2, axis=0)
    names = np.arange(train.shape[1]) if gene_names is None else np.asarray(gene_names)
    if len(names) != train.shape[1]:
        raise ValueError("gene_names length mismatch")
    return {"gene_metrics": {"gene": names, "mse": mse, "baseline_mse": base, "improvement": base - mse},
            "mean_mse": float(mse.mean()), "mean_baseline_mse": float(base.mean()),
            "mean_improvement": float((base - mse).mean())}
