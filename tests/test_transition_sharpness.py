import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from pseudospace import transition_sharpness as ts


def _sig(x, m, w):
    return 1 / (1 + np.exp(-(x - m) * np.log(81) / w))


def test_equal_count_bins_and_binned_log_rates():
    rng = np.random.default_rng(0)
    pos = rng.uniform(0, 1, 400)
    bins = ts.equal_count_bins(pos, 8)
    assert np.bincount(bins).tolist() == [50] * 8
    assert (np.diff([pos[bins == k].max() for k in range(8)]) > 0).all()
    counts = sparse.csr_matrix(np.column_stack([np.full(400, 3), (pos > np.median(pos)) * 8]))
    rates = ts.binned_log_rates(counts, np.full(400, 10.), bins, 8)
    assert rates.shape == (2, 8)
    assert rates[0] == pytest.approx(np.log2((150 + .5) / 500))
    assert (rates[1, :4] < rates[1, 4:] - 5).all()


def test_fit_edges_recovers_a_sharp_step_a_gradient_a_hump_and_flat():
    rng = np.random.default_rng(1)
    x = np.linspace(.1, 2.9, 40)
    profiles = np.array([
        1.2 * _sig(x, 2.0, .15),                           # sharp step at the S2/S3 boundary
        -1.0 * _sig(x, 1.3, 2.2),                          # gradual fall
        _sig(x, 1.0, .3) - _sig(x, 2.0, .3),               # hump over S2
        np.zeros_like(x),                                  # flat
    ]) + rng.normal(0, .03, (4, 40))
    models, edges = ts.fit_edges(x, profiles)
    assert models.model.tolist() == ['M1', 'M1', 'M2', 'flat']
    e = edges.set_index(['profile', 'edge'])
    assert e.loc[(0, 0), 'position'] == pytest.approx(2.0, abs=.05) and e.loc[(0, 0), 'width'] <= .3
    assert e.loc[(1, 0), 'direction'] == -1 and e.loc[(1, 0), 'width'] >= 1.5
    assert e.loc[(2, 0), 'direction'] == 1 and e.loc[(2, 1), 'direction'] == -1
    assert e.loc[(2, 0), 'position'] == pytest.approx(1.0, abs=.08) and e.loc[(2, 1), 'position'] == pytest.approx(2.0, abs=.08)
    assert (e.xs(2, level='profile').width < .5).all()


def test_class_location_and_agreement():
    assert [ts.edge_class(w) for w in (.5, .7, 1., np.nan)] == ['sharp', 'graded', 'gradual', '']
    assert [ts.edge_location(p) for p in (.9, 1.1, 1.5, 2.15, 2.7)] == ['at S1/S2', 'at S1/S2', 'inside S2', 'at S2/S3', 'inside S3']
    a = pd.DataFrame({'position': [1.0, 2.0], 'width': [.3, .4], 'amplitude': [1., -1.], 'direction': [1, -1]})
    b = pd.DataFrame({'position': [1.1, 2.5, 2.1], 'width': [1.4, .4, 2.], 'amplitude': [1., -1., 1.], 'direction': [1, -1, 1]})
    assert ts.match_edges(a, b) == [(0, 0)]                # the falling 2.5 is too far; the rising 2.1 has the wrong sign
    agreed = ts.agreed_edges(a, b)
    assert len(agreed) == 1 and agreed.loc[0, 'class'] == 'mice disagree on width'
    assert agreed.loc[0, 'location'] == 'at S1/S2' and agreed.loc[0, 'position'] == pytest.approx(1.05)


def test_spline_step_response_is_wider_than_the_binned_fit():
    pytest.importorskip('patsy')
    rng = np.random.default_rng(2)
    pos = np.sort(rng.uniform(0, 1, 3000))
    knots = np.quantile(pos, [.25, .5, .75])
    grid = np.linspace(pos.min(), pos.max(), 101)
    curve = ts.spline_step_response(pos, .6, knots, grid)
    spline_w = ts.crossing_width(grid, curve)
    bins = ts.equal_count_bins(pos, 40)
    frac = np.bincount(bins, weights=(pos >= .6).astype(float)) / np.bincount(bins)
    x = np.array([np.median(pos[bins == k]) for k in range(40)])
    _, edges = ts.fit_edges(x, frac[None], flat=.1)
    assert spline_w > .1 and edges.width.iloc[0] < spline_w / 3


def test_edge_composition_labels():
    genes = pd.DataFrame({'gene': list('abcdef'), 'position': [1., 1.3, 1.6, 1.9, 1.2, 3.],
                          'direction': [1, 1, 1, 1, -1, 1], 'class': ['sharp'] * 4 + ['gradual', 'sharp']})
    out = ts.edge_composition({'position': 1.5, 'direction': 1}, genes, list('abcdef'))
    assert out['genes'] == 4 and out['composition'] == 'staggered gene switches'
    genes['class'] = 'gradual'
    assert ts.edge_composition({'position': 1.5, 'direction': 1}, genes, list('abcdef'))['composition'] == 'gene gradient'
    assert ts.edge_composition({'position': 1.5, 'direction': 1}, genes, ['a', 'b'])['composition'] == 'not evaluable'
