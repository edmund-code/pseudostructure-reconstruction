import numpy as np
import pytest
pytest.importorskip("scipy")
pytest.importorskip("sklearn")

from pseudospace.trajectory_benchmark import (
    fit_dpt_axis, fit_scfates_axis, project_neighbor_axis, _trajectory_anatomy, _orient,
)


def test_neighbor_projection_is_frozen_and_subset_invariant():
    train = np.array([[0.], [1.], [2.], [3.]])
    z = np.array([0., .2, .8, 1.])
    query = np.array([[.5], [2.5]])
    first = project_neighbor_axis(train, z, query, k=2)
    altered = project_neighbor_axis(train, z, np.vstack([query, [[100.]]]), k=2)
    np.testing.assert_allclose(first, altered[:2])
    assert np.all((first >= 0) & (first <= 1))


def test_neighbor_projection_averages_exact_duplicate_neighbors_only():
    train = np.array([[0.], [0.], [1.], [2.]])
    z = np.array([.1, .9, .4, .8])
    result = project_neighbor_axis(train, z, np.array([[0.]]), k=3)
    assert result[0] == pytest.approx(.5)


def test_neighbor_assignments_preserve_projection_and_zero_distance_rule():
    train = np.array([[0.], [0.], [1.], [2.]])
    z = np.array([.1, .9, .4, .8])
    query = np.array([[0.], [1.5]])
    result = project_neighbor_axis(train, z, query, k=3, return_distribution=True)
    np.testing.assert_array_equal(result['z_mean'], project_neighbor_axis(train, z, query, k=3))
    np.testing.assert_allclose(result['weights'].sum(axis=1), 1.)
    np.testing.assert_allclose((result['weights'] * result['support']).sum(axis=1), result['z_mean'])
    assert np.count_nonzero(result['weights'][0]) == 2
    assert result['weights'][0].max() == .5
    # Averaging nonlinear responses uses the original assignments, not their mean.
    assert (result['weights'][0] * result['support'][0] ** 2).sum() > result['z_mean'][0] ** 2


def test_neighbor_assignments_allow_empty_queries_and_validate_flag():
    result = project_neighbor_axis(np.array([[0.], [1.]]), np.array([0., 1.]),
                                   np.empty((0, 1)), k=15, return_distribution=True)
    assert result['support'].shape == result['weights'].shape == (0, 2)
    assert result['z_mean'].shape == (0,)
    with pytest.raises(ValueError, match='return_distribution'):
        project_neighbor_axis(np.array([[0.], [1.]]), np.array([0., 1.]),
                              np.array([[.5]]), return_distribution=1)


@pytest.mark.parametrize('kwargs', [
    {'k': 0}, {'k': 1.5},
])
def test_neighbor_projection_rejects_invalid_k(kwargs):
    with pytest.raises(ValueError):
        project_neighbor_axis(np.array([[0.], [1.]]), np.array([0., 1.]),
                              np.array([[.5]]), **kwargs)


def test_neighbor_projection_rejects_misalignment_and_nonfinite():
    with pytest.raises(ValueError):
        project_neighbor_axis(np.array([[0.], [1.]]), np.array([0.]), np.array([[.5]]))
    with pytest.raises(ValueError):
        project_neighbor_axis(np.array([[0.], [1.]]), np.array([0., 1.]), np.array([[np.nan]]))


def test_trajectory_anchor_mask_ignores_unlabeled_anatomy_values():
    labels = np.array([0, 0, 1, 1, 2, 2, -1, -1], dtype=float)
    mask = np.array([1, 1, 1, 1, 1, 1, 0, 0], dtype=bool)
    validated, returned_mask = _trajectory_anatomy(labels, len(labels), mask)
    altered = labels.copy()
    altered[~mask] = [900, np.nan]
    validated_altered, _ = _trajectory_anatomy(altered, len(altered), mask)
    np.testing.assert_array_equal(validated, [0, 0, 1, 1, 2, 2, -1, -1])
    np.testing.assert_array_equal(validated_altered, validated)
    np.testing.assert_array_equal(returned_mask, mask)
    raw = np.array([.1, .2, .3, .4, .7, .8, .5, .6])
    z, flipped = _orient(raw, validated)
    z_altered, flipped_altered = _orient(raw, validated_altered)
    np.testing.assert_array_equal(z_altered, z)
    assert flipped_altered == flipped


@pytest.mark.parametrize('mask', [
    np.array([True]),
    np.array([[True, False]]),
    np.array([1, 0]),
])
def test_trajectory_anchor_mask_requires_matching_boolean_vector(mask):
    with pytest.raises(ValueError, match='anchor_mask'):
        _trajectory_anatomy(np.array([0, 1]), 2, mask)


def test_trajectory_anchor_mask_requires_all_anatomy_segments():
    with pytest.raises(ValueError, match='all labels'):
        _trajectory_anatomy(np.array([0, 0, 1, 1, -1]), 5,
                            np.array([True, True, True, True, False]))


def test_trajectory_anchor_default_remains_strict():
    with pytest.raises(ValueError):
        _trajectory_anatomy(np.array([0, 1, 2, -1]), 4)


@pytest.mark.parametrize('use_mask', [False, True])
def test_scanpy_dpt_smoke_when_installed(use_mask):
    pytest.importorskip('anndata')
    pytest.importorskip('scanpy')
    t = np.linspace(0, 1, 60)
    Y = np.column_stack([t, t**2, np.sin(t * np.pi)])
    anatomy = np.where(t < 1/3, 0, np.where(t < 2/3, 1, 2))
    anchor_mask = np.arange(len(Y)) % 4 != 0 if use_mask else None
    if anchor_mask is not None:
        anatomy[~anchor_mask] = -1
    z, meta = fit_dpt_axis(Y, anatomy, anchor_mask=anchor_mask,
                           n_neighbors=10, seed=2)
    assert np.isfinite(z).all() and np.ptp(z) > 0
    assert meta['root_index'] >= 0
    assert meta['anchor_count'] == (int(anchor_mask.sum()) if use_mask else len(Y))
    assert meta['unlabeled_count'] == (int((~anchor_mask).sum()) if use_mask else 0)
    if use_mask:
        assert anchor_mask[meta['root_index']]


@pytest.mark.parametrize('use_mask', [False, True])
def test_scfates_smoke_when_installed(use_mask):
    pytest.importorskip('anndata')
    pytest.importorskip('scFates')
    t = np.linspace(0, 1, 60)
    Y = np.column_stack([t, t**2, np.sin(t * np.pi)])
    anatomy = np.where(t < 1/3, 0, np.where(t < 2/3, 1, 2))
    anchor_mask = np.arange(len(Y)) % 4 != 0 if use_mask else None
    if anchor_mask is not None:
        anatomy[~anchor_mask] = -1
    z, meta = fit_scfates_axis(Y, anatomy, anchor_mask=anchor_mask,
                               nodes=10, seed=2)
    assert np.isfinite(z).all() and np.ptp(z) > 0
    assert meta['tips'] == 2 and meta['forks'] == 0
    assert meta['anchor_count'] == (int(anchor_mask.sum()) if use_mask else len(Y))
    assert meta['unlabeled_count'] == (int((~anchor_mask).sum()) if use_mask else 0)



def test_path_guard_rejects_disconnected_cycle_with_two_other_tips():
    from types import SimpleNamespace
    from pseudospace.trajectory_benchmark import _scfates_graph_guard
    # A separate 3-cycle plus a two-node path has two tips but is not one path.
    adjacency=np.array([[0,1,1,0,0],[1,0,1,0,0],[1,1,0,0,0],
                        [0,0,0,0,1],[0,0,0,1,0]],float)
    obj=SimpleNamespace(uns={'graph':{'B':adjacency}})
    with pytest.raises(ValueError,match='unbranched'):
        _scfates_graph_guard(obj)
