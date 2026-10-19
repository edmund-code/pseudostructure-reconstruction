import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from pseudospace import contamination as ct


def _two_squares():
    """A 9 x 9 square (label 0) beside a 3 x 3 square (label 1), the rest unassigned (-1)."""
    rr, cc = np.meshgrid(np.arange(11), np.arange(16), indexing='ij')
    rows, cols = rr.ravel(), cc.ravel()
    labels = np.full(rows.size, -1)
    labels[(rows >= 1) & (rows <= 9) & (cols >= 1) & (cols <= 9)] = 0
    labels[(rows >= 4) & (rows <= 6) & (cols >= 10) & (cols <= 12)] = 1
    return rows, cols, labels


def test_bin_grid_and_depth():
    r, c = ct.bin_grid(['s_002um_01022_01191-1', 's_002um_00000_00007-1'])
    assert r.tolist() == [1022, 0] and c.tolist() == [1191, 7]
    with pytest.raises(ValueError):
        ct.bin_grid(['bad'])
    rows, cols, labels = _two_squares()
    depth = ct.label_depth(rows, cols, labels)
    at = lambda i, j: depth[(rows == i) & (cols == j)][0]
    assert at(1, 1) == 0 and at(5, 5) == 4 and at(2, 5) == 1   # centre of the 9 x 9 square is 4 bins in
    assert at(5, 9) == 0 and at(5, 10) == 0 and at(5, 11) == 1  # the shared side is an edge on both squares


def test_aggregate_and_size_matched_sample():
    rows, cols, labels = _two_squares()
    counts = sparse.csr_matrix(np.ones((rows.size, 2)))
    depth = ct.label_depth(rows, cols, labels)
    full = ct.unit_aggregate(counts, labels, 2).toarray()
    core = ct.unit_aggregate(counts, labels, 2, keep=depth >= 2).toarray()
    assert full[:, 0].tolist() == [81, 9] and core[:, 0].tolist() == [25, 0]
    keep = depth >= 2
    draw = ct.size_matched_sample(labels, keep, np.random.default_rng(0))
    assert np.bincount(labels[draw & (labels >= 0)], minlength=2).tolist() == [25, 0]
    assert not draw[labels < 0].any()
    other = ct.size_matched_sample(labels, keep, np.random.default_rng(1))
    assert (draw != other).any()


def test_contact_markers_fraction_and_contrast():
    rows, cols, labels = _two_squares()
    share = ct.contact_share(rows, cols, labels, {0: 'PT', 1: 'AL'}, [0, 1], ['AL', 'outside'])
    assert share.loc[0, 'outside'] == pytest.approx(29 / 32)     # all edge bins but the 3 facing the small square
    assert share.loc[0, 'AL'] == pytest.approx(3 / 32)           # 3 of its 32 edge bins touch the small square
    genes = pd.Index(['Lrp2', 'Umod', 'Both'])
    ref = pd.Series([1., 400., 50.], genes)
    pt = pd.Series([300., 2., 10.], genes)
    assert ct.class_markers(ref, pt).tolist() == ['Umod']
    assert ct.class_markers(ref, pt, others=[pd.Series([0., 100., 0.], genes)]).tolist() == []
    assert ct.nonnative_fraction([.01, .11], .01, 1.01).tolist() == pytest.approx([0, .1])
    m = sparse.csr_matrix(np.array([[10, 0], [10, 0], [0, 10], [5, 5]]))
    with np.errstate(divide='ignore'):
        e = ct.end_contrast(m, [0., .1, .9, 1.], pseudo=0)
    assert e[0] == pytest.approx(-2) and np.isposinf(e[1])         # late pool: 5 of 20 counts; gene 1 absent early


def test_peritubular_ring():
    rows, cols, labels = _two_squares()
    ring = ct.peritubular_ring(rows, cols, labels, max_dist=2)
    at = lambda i, j: ring[(rows == i) & (cols == j)][0]
    assert at(0, 5) == 0 and at(5, 0) == 0                 # outside bins one step from the big square
    assert at(5, 14) == 1 and at(5, 15) == -1              # 2 bins and 3 bins from the small square
    assert (ring[labels >= 0] == -1).all()                 # bins inside polygons are not ring bins
