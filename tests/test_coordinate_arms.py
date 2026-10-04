import numpy as np
import pytest
from scipy import sparse
from scipy.stats import spearmanr

from pseudospace.coordinate_arms import (heldout_prediction, marker_tip_score, registration_gap,
                                         scanpy_lognorm, within_segment_agreement)


def _curve(n=400, seed=0):
    rng = np.random.default_rng(seed)
    t = rng.uniform(size=n)
    X = np.column_stack([np.cos(np.pi * t), np.sin(np.pi * t), .5 * t,
                         rng.normal(0, .02, n), rng.normal(0, .02, n)])
    return t, X + rng.normal(0, .02, X.shape)


def test_scanpy_lognorm_divides_and_sorts():
    counts = sparse.csr_matrix(([3., 1., 4.], ([0, 0, 1], [1, 0, 1])), shape=(2, 2))
    counts.has_sorted_indices = False
    out = scanpy_lognorm(counts)
    np.testing.assert_allclose(out.toarray(), np.log1p(np.array([[1, 3], [0, 4]]) / np.array([[4], [4]]) * 1e4))
    assert out.has_sorted_indices


def test_marker_tip_score_is_scaled_within_species():
    values = np.array([[1., 0.], [3., 0.], [10., 5.], [30., 5.]])
    species = np.array(['m', 'm', 'h', 'h'])
    score = marker_tip_score(values, species, [0], [1])
    np.testing.assert_allclose(score, [-1, 1, -1, 1])


def test_scfates_curve_dpt_and_loso_follow_a_synthetic_path():
    pytest.importorskip('scFates')
    from pseudospace.coordinate_arms import fit_dpt, fit_scfates_curve, loso_positions
    t, X = _curve()
    groups = np.where(np.arange(len(t)) % 2 == 0, 'a', 'b')
    fit = fit_scfates_curve(X, 1 - t, nodes=15, seed=0, agree_groups=groups)
    assert spearmanr(fit['z'], t).statistic > .95
    dpt = fit_dpt(X, fit['root_cell'], n_neighbors=30)
    assert spearmanr(dpt, t, nan_policy='omit').statistic > .9
    specimen = np.array(list('wxyz'))[np.arange(len(t)) % 4]
    projected = loso_positions(X, 1 - t, specimen, 'w', nodes=15, seed=0)
    assert spearmanr(projected, t[specimen == 'w']).statistic > .9


def test_automatic_root_must_be_a_tip():
    pytest.importorskip('scFates')
    from pseudospace.coordinate_arms import fit_scfates_curve
    t, X = _curve(seed=2)
    with pytest.raises(ValueError, match='not a tip'):
        fit_scfates_curve(X, -t - 1, nodes=15, seed=0)


def test_curve_root_requires_group_agreement():
    pytest.importorskip('scFates')
    from pseudospace.coordinate_arms import fit_scfates_curve
    t, X = _curve(seed=1)
    groups = np.where(t < .5, 'a', 'b')  # each group sees one end only, with opposite score signs
    score = np.where(groups == 'a', -t, t)
    with pytest.raises(ValueError, match='disagree'):
        fit_scfates_curve(X, score, nodes=15, seed=0, agree_groups=groups)


def test_registration_gap_and_within_segment_agreement():
    z = np.array([.1, .2, .5, .6, .3, .4, .5, .6])
    species = np.array(['mouse'] * 4 + ['human'] * 4)
    specimen = np.array(['m1', 'm1', 'm2', 'm2', 'h1', 'h1', 'h2', 'h2'])
    segment = np.array(['S1', 'S1', 'S1', 'S1', 'S1', 'S1', 'S1', 'S1'])
    gap, worst = registration_gap(z, species, specimen, segment)
    assert gap['S1'] == pytest.approx(.45 - .35) and worst == pytest.approx(.1)
    cells = np.repeat(['a', 'b'], 5)
    table = within_segment_agreement(np.arange(10.0), np.arange(10.0) ** 2, cells, np.repeat('S1', 10))
    assert np.allclose(table.rho, 1) and len(table) == 2


def test_heldout_prediction_prefers_the_coordinate_for_a_smooth_gene():
    rng = np.random.default_rng(3)
    z_train, z_test = rng.uniform(size=300), rng.uniform(size=300)
    seg = lambda z: np.where(z < .33, 'S1', np.where(z < .66, 'S2', 'S3'))
    gene = lambda z: np.column_stack([np.sin(3 * z), np.where(z < .33, 0., 1.)])
    out = heldout_prediction(gene(z_train), z_train, seg(z_train), gene(z_test), z_test, seg(z_test))
    assert out['coordinate_mse'][0] < .2 * out['segment_mse'][0]
    assert out['segment_mse'][1] == pytest.approx(0)
