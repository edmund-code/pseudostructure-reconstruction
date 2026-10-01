"""Small synthetic checks for inductive repeated-structure baselines."""
import numpy as np
import pytest

pytest.importorskip('scipy')
pytest.importorskip('sklearn')

from pseudospace.repeated_baselines import (
    anatomy_initialize,
    baseline_coordinates,
    curve_agreement,
    topology_coordinate,
)


def test_pca1_is_training_oriented_and_holdout_uses_frozen_training_scale():
    rng = np.random.default_rng(42)
    z = np.linspace(0, 1, 90)
    X = np.column_stack([z, 2*z, -z]) + rng.normal(0, .01, (90, 3))
    anatomy = np.repeat(['S1', 'S2', 'S3'], 30)
    train = np.arange(0, 90, 2)
    test = np.arange(1, 90, 2)
    a, b = baseline_coordinates(X[train], X[test], anatomy[train], 'pca1')
    assert spearman(a, z[train]) > .98
    assert np.all((a >= 0) & (a <= 1)) and np.all((b >= 0) & (b <= 1))
    # No test sample can change the training coordinate transform.
    _, b2 = baseline_coordinates(X[train], X[test] + 100, anatomy[train], 'pca1')
    np.testing.assert_allclose(a, baseline_coordinates(X[train], X[test] + 100, anatomy[train])[0])
    assert not np.allclose(b, b2)


def test_panel_dpt_fits_training_graph_and_projects_heldout_samples():
    pytest.importorskip('scanpy')
    pytest.importorskip('anndata')
    rng = np.random.default_rng(231)
    z = np.linspace(0, 1, 72)
    X = np.column_stack([z, np.sin(2*np.pi*z), np.cos(2*np.pi*z)])
    X += rng.normal(0, .025, X.shape)
    anatomy = np.repeat(['S1', 'S2', 'S3'], 24)
    ztrain, ztest = baseline_coordinates(X[::2], X[1::2], anatomy[::2], 'panel_dpt', seed=4)
    assert ztrain.shape == (36,) and ztest.shape == (36,)
    assert np.isfinite(ztrain).all() and np.isfinite(ztest).all()
    assert np.all((ztrain >= 0) & (ztrain <= 1))
    assert np.all((ztest >= 0) & (ztest <= 1))
    assert spearman(ztrain, z[::2]) > .4


def spearman(a, b):
    from scipy.stats import spearmanr
    return spearmanr(a, b).statistic


def test_anatomy_initialization_keeps_coarse_intervals_and_has_within_interval_ranks():
    rng = np.random.default_rng(1)
    anatomy = np.repeat(['S1', 'S2', 'S3'], 20)
    X = rng.normal(size=(60, 5))
    z = anatomy_initialize(X, anatomy)
    assert np.all((z[:20] >= 0) & (z[:20] < 1/3))
    assert np.all((z[20:40] >= 1/3) & (z[20:40] < 2/3))
    assert np.all((z[40:] >= 2/3) & (z[40:] <= 1))
    assert len(np.unique(z[:20])) == 20


def test_curve_agreement_recovers_shared_smooth_gene_shape():
    z = np.tile(np.linspace(0, 1, 40), 3)
    specimens = np.repeat(['a', 'b', 'c'], 40)
    y = np.column_stack([np.sin(np.pi*z), np.cos(np.pi*z)])
    out = curve_agreement(y, z, specimens, ['sin', 'cos'])
    assert len(out) == 6
    assert out.curve_pearson.min() > .95
    assert out.monotonic_sign_agreement.notna().all()


def test_topology_coordinate_does_not_bridge_component_without_both_endpoints():
    rng = np.random.default_rng(9)
    # Three adjacent clouds make an S1-S2-S3 path; a remote S2 component has
    # no endpoint and must remain unmapped instead of being bridged.
    x = np.r_[np.linspace(0, 1, 36), np.linspace(10, 10.1, 8)][:, None]
    x += rng.normal(0, .001, x.shape)
    anatomy = np.array(['S1']*12 + ['S2']*12 + ['S3']*12 + ['S2']*8)
    z, diagnostics = topology_coordinate(x, anatomy, k=5)
    assert diagnostics['n_components'] >= 2
    assert np.isnan(z[-8:]).all()
    assert np.all(z[:12] == 0) and np.all(z[24:36] == 1)
    assert np.isfinite(z[12:24]).all()
    assert np.all((z[np.isfinite(z)] >= 0) & (z[np.isfinite(z)] <= 1))
