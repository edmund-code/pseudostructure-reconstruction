"""Small kernel mean atlas for cross-sectional profile projection.

The atlas is a descriptive molecular reference curve. Coordinates are supplied by
the caller during fitting; projection only sees molecular profiles and the frozen
atlas, so it does not consume anatomy labels or a precomputed trajectory.
"""
from __future__ import annotations

import numpy as np


def _matrix(values, name):
    x = np.asarray(values, dtype=float)
    if x.ndim != 2 or not x.shape[0] or not x.shape[1]:
        raise ValueError(f"{name} must be a non-empty 2D matrix")
    if not np.isfinite(x).all():
        raise ValueError(f"{name} must contain only finite values")
    return x


def _coordinates(z, n, name="z"):
    values = np.asarray(z, dtype=float).reshape(-1)
    if len(values) != n or not np.isfinite(values).all():
        raise ValueError(f"{name} must contain one finite value per profile")
    if np.any((values < 0) | (values > 1)):
        raise ValueError(f"{name} must lie in [0, 1]")
    return values


def _grid(grid):
    values = np.linspace(0, 1, 101) if grid is None else np.asarray(grid, dtype=float)
    if (values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all()
            or np.any(np.diff(values) <= 0) or values[0] != 0 or values[-1] != 1):
        raise ValueError("grid must be a finite, increasing vector spanning [0, 1]")
    return values


def _bandwidth(value):
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError("bandwidth must be positive and finite")
    return value


def _kernel_weights(query, training, bandwidth):
    """Gaussian row weights, evaluated stably even for very small bandwidths."""
    distance = np.abs(np.asarray(query)[:, None] - np.asarray(training)[None, :])
    minimum = distance.min(axis=1, keepdims=True)
    gap = distance - minimum
    # Subtracting the row maximum of log weights is equivalent to subtracting
    # the smallest squared distance. Compute the difference of squares in log
    # space so distant points do not overflow before that subtraction.
    log_contrast = np.zeros_like(distance)
    mask = gap > 0
    if mask.any():
        log_value = (
            np.log(gap[mask]) + np.log((distance + minimum)[mask])
            - 2.0 * np.log(bandwidth) - np.log(2.0)
        )
        log_contrast[mask] = np.exp(np.minimum(log_value, np.log(745.0)))
    log_weights = -log_contrast
    log_weights -= log_weights.max(axis=1, keepdims=True)
    return np.exp(log_weights)


def _distances(y, atlas, offset=None):
    curve = atlas["mean"]
    if offset is not None:
        curve = curve + np.asarray(offset, dtype=float)[None, :]
    # Block over queries so a large projection does not materialize the full
    # query-by-grid-by-feature difference tensor.
    output = np.empty((len(y), len(curve)), dtype=float)
    curve_norm = np.sum(curve * curve, axis=1)
    for start in range(0, len(y), 256):
        block = y[start:start + 256]
        output[start:start + 256] = np.maximum(
            np.sum(block * block, axis=1)[:, None] + curve_norm[None, :]
            - 2 * block @ curve.T,
            0.0,
        )
    return output


def fit_mean_atlas(Ytrain, ztrain, grid=None, bandwidth=0.08, min_effective=12):
    """Fit a Gaussian-kernel-smoothed mean profile at each canonical grid point.

    ``min_effective`` is checked at every grid point using the normalized-kernel
    effective sample size, ``1 / sum(w**2)``. The residual variance is one
    isotropic estimate from training profiles' nearest fitted grid mean.
    """
    y = _matrix(Ytrain, "Ytrain")
    z = _coordinates(ztrain, len(y), "ztrain")
    grid = _grid(grid)
    bandwidth = _bandwidth(bandwidth)
    min_effective = float(min_effective)
    if not np.isfinite(min_effective) or min_effective <= 0:
        raise ValueError("min_effective must be positive and finite")

    weights = _kernel_weights(grid, z, bandwidth)
    totals = weights.sum(axis=1)
    normalized = weights / totals[:, None]
    n_effective = 1.0 / np.sum(normalized * normalized, axis=1)
    unsupported = np.flatnonzero(n_effective < min_effective)
    if unsupported.size:
        raise ValueError(
            f"insufficient effective training support at {len(unsupported)} grid points "
            f"(minimum {n_effective.min():.2f}, required {min_effective:g})"
        )

    mean = normalized @ y
    provisional = {"grid": grid, "mean": mean}
    nearest = np.argmin(_distances(y, provisional), axis=1)
    residual = y - mean[nearest]
    sigma2 = max(float(np.mean(residual * residual)), 1e-8)
    return {
        "grid": grid.copy(),
        "mean": mean,
        "bandwidth": bandwidth,
        "min_effective": min_effective,
        "n_effective": n_effective,
        "sigma2": sigma2,
    }


def project_mean_atlas(atlas, Y):
    """Project profiles to the nearest atlas mean and return grid posterior summaries."""
    y = _matrix(Y, "Y")
    if y.shape[1] != atlas["mean"].shape[1]:
        raise ValueError("Y feature count does not match atlas")
    d2 = _distances(y, atlas)
    map_index = np.argmin(d2, axis=1)
    logits = -d2 / (2.0 * float(atlas["sigma2"]))
    logits -= logits.max(axis=1, keepdims=True)
    probabilities = np.exp(logits)
    probabilities /= probabilities.sum(axis=1, keepdims=True)
    entropy = -np.sum(probabilities * np.log(np.maximum(probabilities, 1e-300)), axis=1)
    residual_norm = np.sqrt(d2[np.arange(len(y)), map_index])
    return {
        "z_MAP": atlas["grid"][map_index],
        "z_mean": probabilities @ atlas["grid"],
        "posterior": probabilities,
        "entropy": entropy,
        "residual_norm": residual_norm,
    }


def fit_global_offset(atlas, calibration_Y, penalty=4.0, iterations=5):
    """Estimate one shrunk global feature offset using calibration profiles only."""
    y = _matrix(calibration_Y, "calibration_Y")
    if y.shape[1] != atlas["mean"].shape[1]:
        raise ValueError("calibration_Y feature count does not match atlas")
    penalty = float(penalty)
    if not np.isfinite(penalty) or penalty < 0:
        raise ValueError("penalty must be finite and non-negative")
    if not isinstance(iterations, (int, np.integer)) or iterations < 1:
        raise ValueError("iterations must be a positive integer")

    offset = np.zeros(y.shape[1], dtype=float)
    history = []
    for _ in range(iterations):
        d2 = _distances(y, atlas, offset=offset)
        positions = np.argmin(d2, axis=1)
        residual = y - atlas["mean"][positions]
        proposed = residual.mean(axis=0) / (1.0 + penalty)
        history.append({"offset_norm": float(np.linalg.norm(proposed)),
                        "position_changes": int(np.count_nonzero(positions != history[-1]["positions"]))
                        if history else 0,
                        "positions": positions})
        offset = proposed
    # Keep the history compact and serializable: positions are useful during the
    # alternating fit, but are not part of the public history contract.
    for item in history:
        item.pop("positions")
    return offset, history


def conditional_gene_predictions(Xtrain, ztrain, query_positions, bandwidth=0.08):
    """Kernel-predict gene values from training profiles in bounded query blocks."""
    x = _matrix(Xtrain, "Xtrain")
    z = _coordinates(ztrain, len(x), "ztrain")
    query = np.asarray(query_positions, dtype=float).reshape(-1)
    if not np.isfinite(query).all() or np.any((query < 0) | (query > 1)):
        raise ValueError("query_positions must be finite and lie in [0, 1]")
    bandwidth = _bandwidth(bandwidth)
    output = np.empty((len(query), x.shape[1]), dtype=float)
    for start in range(0, len(query), 256):
        positions = query[start:start + 256]
        weights = _kernel_weights(positions, z, bandwidth)
        weights /= weights.sum(axis=1, keepdims=True)
        output[start:start + len(positions)] = weights @ x
    return output


def segment_means(Xtrain, labels_train, labels_query):
    """Return segment-label mean predictions and training counts for S1/S2/S3."""
    x = _matrix(Xtrain, "Xtrain")
    train = np.asarray(labels_train).reshape(-1)
    query = np.asarray(labels_query).reshape(-1)
    required = ("S1", "S2", "S3")
    if len(train) != len(x):
        raise ValueError("labels_train must contain one label per training profile")
    if not set(train).issubset(required) or not set(query).issubset(required):
        raise ValueError("segment labels must be drawn from S1, S2, and S3")
    counts = {label: int(np.sum(train == label)) for label in required}
    if any(count == 0 for count in counts.values()):
        raise ValueError("training data must support each of S1, S2, and S3")
    means = {label: x[train == label].mean(axis=0) for label in required}
    prediction = np.vstack([means[label] for label in query]) if len(query) else np.empty((0, x.shape[1]))
    return {"prediction": prediction, "means": means, "counts": counts}
