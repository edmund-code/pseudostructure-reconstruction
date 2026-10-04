"""Mean-path registration for partial observations of an ordered structure.

The query may occupy any ordered subset of the reference. No equal occupancy,
matched endpoints, physical arc length, or calibrated position uncertainty is
implied. Every query mean is forced to match; use only as a healthy pilot until
unmatched molecular states can be handled and validated explicitly.
"""

import numpy as np


def _means(value, name):
    value = np.asarray(value, dtype=float)
    if value.ndim != 2 or not all(value.shape) or not np.isfinite(value).all():
        raise ValueError(f'{name} must be a nonempty finite two-dimensional matrix')
    return value


def _grid(value, n=None):
    value = np.asarray(value, dtype=float)
    if (value.ndim != 1 or len(value) < 2 or
            (n is not None and len(value) != n) or
            not np.isfinite(value).all() or np.any(np.diff(value) <= 0) or
            np.any((value < 0) | (value > 1))):
        raise ValueError('grid must be strictly increasing, finite, within [0, 1], and match means')
    return value


def local_curve_means(Y, z, grid, *, bandwidth=.08):
    """Kernel means and effective support on a native curve, without occupancy matching."""
    Y = _means(Y, 'Y')
    z = np.asarray(z, dtype=float)
    grid = _grid(grid)
    if (z.ndim != 1 or len(z) != len(Y) or not np.isfinite(z).all() or
            np.any((z < 0) | (z > 1))):
        raise ValueError('z must match rows and be finite within [0, 1]')
    if (isinstance(bandwidth, (bool, np.bool_)) or not np.isscalar(bandwidth) or
            not np.isfinite(bandwidth) or bandwidth <= 0):
        raise ValueError('bandwidth must be a positive finite scalar')
    weights = np.exp(-.5 * ((grid[:, None] - z[None, :]) / bandwidth) ** 2)
    total = weights.sum(axis=1)
    if np.any(total <= 0):
        raise ValueError('no local support for one or more query knots')
    weights /= total[:, None]
    return weights @ Y, 1. / np.sum(weights ** 2, axis=1)


def calibration_support(calibration_Y, query_Y, *, k=15, quantile=.95):
    """Return an empirical molecular support gate, not physical overlap probability.

    The calibration radius uses its kth *other* neighbor, excluding one self
    match. Query distances use k calibration neighbors. This permits frozen
    reference fallback when a calibration cloud covers only part of an object.
    """
    from sklearn.neighbors import NearestNeighbors

    calibration_Y = _means(calibration_Y, 'calibration_Y')
    query_Y = np.asarray(query_Y, dtype=float)
    if (query_Y.ndim != 2 or query_Y.shape[1] != calibration_Y.shape[1] or
            not np.isfinite(query_Y).all()):
        raise ValueError('query_Y must be finite and match calibration features')
    if (not isinstance(k, (int, np.integer)) or isinstance(k, (bool, np.bool_)) or
            k < 1 or len(calibration_Y) <= k):
        raise ValueError('k must be positive and smaller than the calibration row count')
    if (isinstance(quantile, (bool, np.bool_)) or not np.isscalar(quantile) or
            not np.isfinite(quantile) or not 0 < quantile <= 1):
        raise ValueError('quantile must be finite within (0, 1]')
    model = NearestNeighbors().fit(calibration_Y)
    other_distances = model.kneighbors(calibration_Y, n_neighbors=k + 1)[0][:, -1]
    radius = float(np.quantile(other_distances, quantile))
    query_distances = (model.kneighbors(query_Y, n_neighbors=k)[0][:, -1]
                       if len(query_Y) else np.empty(0))
    return dict(radius=radius, supported=query_distances <= radius,
                neighbor_distance=query_distances)


def _ordered_indices(cost):
    """Exact minimum cost with nondecreasing indices and free reference endpoints."""
    n, m = cost.shape
    previous = cost[0].copy()
    parents = np.zeros((n, m), dtype=np.int64)
    for i in range(1, n):
        following = np.empty(m)
        best = 0
        for j in range(m):
            if previous[j] < previous[best]:
                best = j
            parents[i, j] = best
            following[j] = cost[i, j] + previous[best]
        previous = following
    indices = np.empty(n, dtype=np.int64)
    indices[-1] = int(np.argmin(previous))
    for i in range(n - 1, 0, -1):
        indices[i - 1] = parents[i, indices[i]]
    return indices


def register_ordered_means(query_means, reference_means, reference_grid):
    """Choose query orientation and an ordered, possibly partial reference match.

    One reference knot is selected for every query knot. Repeated knots and
    arbitrary forward jumps are allowed. Orientation and path use modeling
    expression only; the independent nearest-mean result is an order ablation.
    Ties choose the forward orientation and earliest available reference knot.
    """
    query_means = _means(query_means, 'query_means')
    reference_means = _means(reference_means, 'reference_means')
    if query_means.shape[1] != reference_means.shape[1]:
        raise ValueError('query and reference means must have matching features')
    reference_grid = _grid(reference_grid, len(reference_means))
    cost = np.mean((query_means[:, None, :] - reference_means[None, :, :]) ** 2, axis=2)
    if not np.isfinite(cost).all():
        raise ValueError('matching cost must be finite')
    forward = _ordered_indices(cost)
    reverse = _ordered_indices(cost[::-1])[::-1]
    rows = np.arange(len(cost))
    forward_cost = float(cost[rows, forward].mean())
    reverse_cost = float(cost[rows, reverse].mean())
    orientation = -1 if reverse_cost < forward_cost else 1
    indices = reverse if orientation == -1 else forward
    independent = np.argmin(cost, axis=1)
    z = reference_grid[indices]
    return dict(reference_index=indices, z=z, orientation=orientation,
                cost=min(forward_cost, reverse_cost), forward_cost=forward_cost,
                reverse_cost=reverse_cost, orientation_gap=abs(forward_cost-reverse_cost),
                independent_z=reference_grid[independent],
                independent_cost=float(cost[rows, independent].mean()),
                reference_span=float(np.ptp(z)),
                distinct_reference_knots=int(len(np.unique(indices))),
                plateau_fraction=float(np.mean(np.diff(indices) == 0)) if len(z) > 1 else 1.,
                max_reference_jump=float(np.max(np.abs(np.diff(z)))) if len(z) > 1 else 0.)
