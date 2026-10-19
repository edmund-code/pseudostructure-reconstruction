"""Bin-level contamination audit of segmented units (notebook 104; docs/results/pt-contamination-audit.md).

Units are sums of the 2 µm bins inside a segmentation polygon. Transcripts that diffuse or are mis-assigned across a
polygon edge add the neighbours' signal to a unit. These helpers rebuild units from their bins at a chosen depth from
the edge, so a statistic can be recomputed on the interior bins only.

* ``bin_grid``: the (row, column) of each Visium HD bin from its barcode ``s_002um_RRRRR_CCCCC-1``.
* ``label_depth``: each bin's Euclidean distance, in bins, to the nearest edge bin. An edge bin has a 4-neighbour with
  another label, or off the tissue; edge bins have depth 0.
* ``unit_aggregate``: units x genes sums of the bins kept.
* ``size_matched_sample``: for each unit, as many of its bins, drawn uniformly, as a given mask keeps in it.
* ``contact_share``: per unit, the share of its edge bins touching each neighbour class.
* ``class_markers``: genes much higher in a reference profile than in the PT interior (and than other references).
* ``nonnative_fraction``: the share of a unit's counts attributable to a class from its marker share.
* ``log_cp10k``, ``end_contrast``: pooled log2 CP10k and the late-minus-early contrast of pooled units.

Ponytail: depth is measured on the bin grid (to within one 2 µm bin), not on the polygon geometry; the upgrade is a
point-to-polygon distance from the GeoJSON.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import ndimage, sparse

_OFF = -2   # raster value off the tissue


def bin_grid(barcodes):
    """(rows, cols) integer arrays from barcodes like ``s_002um_01022_01191-1``."""
    parts = pd.Series(np.asarray(barcodes).astype(str)).str.extract(r's_\d+um_(\d+)_(\d+)')
    if parts.isna().any().any():
        raise ValueError('barcodes must look like s_002um_RRRRR_CCCCC-1')
    return parts[0].astype(int).to_numpy(), parts[1].astype(int).to_numpy()


def _raster(rows, cols, labels):
    grid = np.full((int(rows.max()) + 3, int(cols.max()) + 3), _OFF, dtype=np.int64)
    grid[rows + 1, cols + 1] = labels       # one-bin margin so the image border is off the tissue
    return grid


def label_depth(rows, cols, labels):
    """Distance (bins) from each bin to the nearest bin with a 4-neighbour of a different label (or off tissue)."""
    rows, cols, labels = (np.asarray(v) for v in (rows, cols, labels))
    grid = _raster(rows, cols, labels)
    edge = np.zeros(grid.shape, bool)
    for axis in (0, 1):
        change = np.diff(grid, axis=axis) != 0
        if axis == 0:
            edge[:-1] |= change
            edge[1:] |= change
        else:
            edge[:, :-1] |= change
            edge[:, 1:] |= change
    depth = ndimage.distance_transform_edt(~edge)
    return depth[rows + 1, cols + 1]


def unit_aggregate(counts, unit_of_bin, n_units, keep=None):
    """Units x genes CSR sums of the bins with ``unit_of_bin >= 0`` (and ``keep``)."""
    unit = np.asarray(unit_of_bin)
    use = unit >= 0 if keep is None else (unit >= 0) & np.asarray(keep, bool)
    idx = np.flatnonzero(use)
    a = sparse.csr_matrix((np.ones(len(idx)), (unit[idx], idx)), shape=(n_units, len(unit)))
    return (a @ sparse.csr_matrix(counts)).tocsr()


def size_matched_sample(unit_of_bin, keep, rng):
    """A mask drawing, within each unit, as many of its bins as ``keep`` holds there, uniformly at random."""
    unit, keep = np.asarray(unit_of_bin), np.asarray(keep, bool)
    out = np.zeros(len(unit), bool)
    inside = np.flatnonzero(unit >= 0)
    order = inside[np.lexsort((rng.random(len(inside)), unit[inside]))]   # bins grouped by unit, shuffled within
    u = unit[order]
    start = np.r_[0, np.flatnonzero(np.diff(u)) + 1]
    rank = np.arange(len(order)) - np.repeat(start, np.diff(np.r_[start, len(order)]))
    quota = np.bincount(unit[keep & (unit >= 0)], minlength=int(unit.max()) + 1)
    out[order[rank < quota[u]]] = True
    return out


def contact_share(rows, cols, labels, label_class, units, classes):
    """For each label in ``units``: share of its edge bins with a 4-neighbour of each class.

    ``label_class`` maps a label (int >= 0) to its class name; ``-1`` (no polygon) is class 'outside'.
    Returns a units x classes DataFrame; a bin touching two classes counts for both.
    """
    rows, cols, labels = (np.asarray(v) for v in (rows, cols, labels))
    grid = _raster(rows, cols, labels)
    lookup = {**{int(k): v for k, v in label_class.items()}, -1: 'outside'}
    code = {c: i for i, c in enumerate(classes)}
    want = np.isin(labels, list(units))
    r, c, own = rows[want] + 1, cols[want] + 1, labels[want]
    hits = np.zeros((len(own), len(classes)), bool)
    is_edge = np.zeros(len(own), bool)
    for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        nb = grid[r + dr, c + dc]
        diff = nb != own
        is_edge |= diff
        names = pd.Series(nb).map(lookup).to_numpy()
        for name, i in code.items():
            hits[:, i] |= diff & (names == name)
    frame = pd.DataFrame(hits[is_edge], columns=list(classes))
    frame['unit'] = own[is_edge]
    share = frame.groupby('unit').mean()
    return share.reindex(list(units)).fillna(0.0)


def class_markers(ref_cp10k, pt_cp10k, *, min_cp10k=1.0, ratio=20.0, others=(), other_ratio=5.0):
    """Genes with ref >= ``min_cp10k``, ref / PT interior >= ``ratio`` and ref / each other profile >= ``other_ratio``.

    Profiles are pandas Series of CP10k on the same gene index.
    """
    ok = (ref_cp10k >= min_cp10k) & (ref_cp10k >= ratio * pt_cp10k)
    for other in others:
        ok &= ref_cp10k >= other_ratio * other
    return ref_cp10k.index[ok.to_numpy()]


def nonnative_fraction(unit_share, pt_share, ref_share):
    """(unit marker share - PT interior share) / (reference share - PT interior share)."""
    return (np.asarray(unit_share, float) - pt_share) / (ref_share - pt_share)


def log_cp10k(counts, pseudo=.5):
    """log2 CP10k of the pooled rows of ``counts`` (genes vector), with ``pseudo`` counts added."""
    total = np.asarray(counts.sum(axis=0)).ravel().astype(float)
    return np.log2((total + pseudo) / max(total.sum(), 1.0) * 1e4)


def end_contrast(counts, x, *, lo=.3, hi=.7, pseudo=.5):
    """log2 CP10k of the pooled units with x >= ``hi`` minus that of the units with x <= ``lo``."""
    x = np.asarray(x, float)
    m = sparse.csr_matrix(counts)
    return log_cp10k(m[x >= hi], pseudo) - log_cp10k(m[x <= lo], pseudo)
