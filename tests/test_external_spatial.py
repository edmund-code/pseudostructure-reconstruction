import numpy as np
import pandas as pd
import pytest

from pseudospace import external_spatial as es


def test_grid_and_edge_distance():
    x, y = es.visium_v1_um([0, 1], [0, 1])
    assert x.tolist() == [0, 50] and y[1] == pytest.approx(86.6, abs=.1)
    rows = np.repeat(np.arange(6), 12)
    cols = np.tile(np.arange(12), 6) * 2 + np.repeat(np.arange(6) % 2, 12)   # hexagonal v1 layout
    tissue = (rows >= 1) & (rows <= 4) & (cols >= 4) & (cols <= 18)
    d = es.edge_distance(rows[tissue], cols[tissue], rows, cols)
    assert d.min() == pytest.approx(100, abs=1e-6) and d.max() > 100   # interior spots lie farther from the edge
    with pytest.raises(ValueError):
        es.edge_distance(rows, cols, rows, cols)


def test_markers_nnls_rank_slope():
    genes = [f'g{i}' for i in range(6)]
    prof = pd.DataFrame({'A': [10, 8, 0, 0, 1, 1], 'B': [0, 0, 9, 7, 1, 1]}, index=genes, dtype=float)
    m = es.type_markers(prof, n=2)
    assert m['A'] == ['g0', 'g1'] and m['B'] == ['g2', 'g3']
    spots = np.array([0.25 * prof.A + 0.75 * prof.B, prof.A, np.zeros(6)])
    w = es.nnls_weights(spots, prof)
    assert w[0] == pytest.approx([.25, .75]) and w[1] == pytest.approx([1, 0]) and w[2].tolist() == [0, 0]
    pos = es.rank_position([3., 1., 2., 10., 20.], ['a', 'a', 'a', 'b', 'b'])
    assert pos.tolist() == pytest.approx([5 / 6, 1 / 6, 1 / 2, 1 / 4, 3 / 4])
    assert es.origin_slope([2., 4.], [1., 2.]) == pytest.approx(2)
