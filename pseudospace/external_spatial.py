"""Placing spots of an external Visium v1 dataset on the PT axis (notebook 105; docs/results/pt-external-spatial-iri.md).

* ``visium_v1_um``: grid (row, column) -> µm; neighbouring spots in a row are two columns and 100 µm apart, rows
  86.6 µm apart.
* ``edge_distance``: each tissue spot's distance (µm) to the nearest grid position without tissue.
* ``type_markers``: per reference type, the genes with the highest ratio of its CP10k to the maximum over the other
  types.
* ``nnls_weights``: non-negative least squares of each spot's CP10k on the types' CP10k profiles, scaled to sum 1.
* ``rank_position``: within-group rank of a score, scaled to (0, 1).
* ``origin_slope``: least-squares slope through the origin.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.spatial import cKDTree

COLUMN_UM, ROW_UM = 50.0, 100.0 * np.sqrt(3) / 2


def visium_v1_um(rows, cols):
    """(x, y) in µm for Visium v1 array rows and columns."""
    return np.asarray(cols, float) * COLUMN_UM, np.asarray(rows, float) * ROW_UM


def edge_distance(tissue_rows, tissue_cols, grid_rows, grid_cols):
    """Distance (µm) from each tissue spot to the nearest grid position that is not a tissue spot."""
    tissue = set(zip(np.asarray(tissue_rows).tolist(), np.asarray(tissue_cols).tolist()))
    empty = np.array([(r, c) for r, c in zip(np.asarray(grid_rows).tolist(), np.asarray(grid_cols).tolist())
                      if (r, c) not in tissue])
    if len(empty) == 0:
        raise ValueError('every grid position holds tissue; no edge to measure from')
    ex, ey = visium_v1_um(empty[:, 0], empty[:, 1])
    tx, ty = visium_v1_um(tissue_rows, tissue_cols)
    d, _ = cKDTree(np.c_[ex, ey]).query(np.c_[tx, ty])
    return d


def type_markers(profiles, *, n=50, min_cp10k=1.0):
    """{type: genes} from a genes x types CP10k DataFrame: top ``n`` by type / max(other types) among CP10k >= min."""
    out = {}
    for t in profiles.columns:
        other = profiles.drop(columns=t).max(axis=1)
        ratio = (profiles[t] + 1e-3) / (other + 1e-3)
        ok = profiles[t] >= min_cp10k
        out[t] = ratio[ok].sort_values(ascending=False).index[:n].tolist()
    return out


def nnls_weights(spots, profiles):
    """Spots x genes CP10k and genes x types CP10k (same genes) -> spots x types weights summing to 1 (0 if none)."""
    a = np.asarray(profiles, float)
    w = np.vstack([nnls(a, np.asarray(row, float))[0] for row in np.asarray(spots, float)])
    total = w.sum(axis=1, keepdims=True)
    return np.divide(w, total, out=np.zeros_like(w), where=total > 0)


def rank_position(score, group):
    """Within each group, (rank - 0.5) / n of ``score``."""
    s = pd.Series(np.asarray(score, float))
    g = pd.Series(np.asarray(group))
    r = s.groupby(g).rank(method='average')
    n = s.groupby(g).transform('size')
    return ((r - .5) / n).to_numpy()


def origin_slope(y, x):
    """Least-squares slope of y on x through the origin."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    return float((x * y).sum() / (x * x).sum())
