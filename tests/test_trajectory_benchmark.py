import numpy as np
import pytest
pytest.importorskip("scipy")
pytest.importorskip("sklearn")

from pseudospace.trajectory_benchmark import (
    fit_dpt_axis, fit_scfates_axis, fit_unoriented_scfates_axis,
    project_neighbor_axis, _trajectory_anatomy, _orient,
    _anatomy_curve_initialization,
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


def test_anatomy_curve_initialization_uses_ordered_anchor_only_centroids():
    Y = np.array([[0., 0.], [2., 0.], [10., 1.], [12., 1.],
                  [20., 2.], [22., 2.], [1000., -1000.]])
    labels = np.array([0, 0, 1, 1, 2, 2, -1])
    anchors = np.array([True, True, True, True, True, True, False])
    first = _anatomy_curve_initialization(Y, labels, anchors)
    np.testing.assert_array_equal(first['InitNodePositions'],
                                  [[1., 0.], [11., 1.], [21., 2.]])
    np.testing.assert_array_equal(first['InitEdges'], [[0, 1], [1, 2]])
    assert first['InitNodes'] == 3
    Y[-1] = [-10000., 10000.]
    labels[-1] = 0
    second = _anatomy_curve_initialization(Y, labels, anchors)
    np.testing.assert_array_equal(second['InitNodePositions'], first['InitNodePositions'])


def test_anatomy_curve_initialization_rejects_duplicate_centroids():
    Y = np.array([[0.], [1.], [0.], [1.], [2.], [3.]])
    labels = np.repeat([0, 1, 2], 2)
    with pytest.raises(ValueError, match='distinct'):
        _anatomy_curve_initialization(Y, labels, np.ones(6, dtype=bool))


@pytest.mark.parametrize('initialization', ['invalid', None, 1])
def test_scfates_rejects_invalid_initialization_before_import(initialization):
    with pytest.raises(ValueError, match='initialization'):
        fit_scfates_axis(np.ones((3, 1)), np.array([0, 1, 2]),
                         initialization=initialization)


@pytest.mark.parametrize('Y', [np.empty((0, 2)), np.ones(3), np.array([[np.nan]])])
def test_unoriented_scfates_rejects_invalid_matrix(Y):
    with pytest.raises(ValueError):
        fit_unoriented_scfates_axis(Y)


@pytest.mark.parametrize('seed', [True, 1.5, -1])
def test_unoriented_scfates_rejects_invalid_seed(seed):
    with pytest.raises(ValueError, match='seed'):
        fit_unoriented_scfates_axis(np.ones((5, 2)), seed=seed)


@pytest.mark.parametrize('nodes', [0, 2, 3.5, True])
def test_unoriented_scfates_rejects_invalid_nodes(nodes):
    with pytest.raises(ValueError, match='nodes'):
        fit_unoriented_scfates_axis(np.ones((5, 2)), nodes=nodes)


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
@pytest.mark.parametrize('initialization', ['default', 'anatomy'])
def test_scfates_smoke_when_installed(use_mask, initialization):
    pytest.importorskip('anndata')
    pytest.importorskip('scFates')
    t = np.linspace(0, 1, 60)
    Y = np.column_stack([t, t**2, np.sin(t * np.pi)])
    anatomy = np.where(t < 1/3, 0, np.where(t < 2/3, 1, 2))
    anchor_mask = np.arange(len(Y)) % 4 != 0 if use_mask else None
    if anchor_mask is not None:
        anatomy[~anchor_mask] = -1
    z, meta = fit_scfates_axis(Y, anatomy, anchor_mask=anchor_mask,
                               nodes=10, seed=2, initialization=initialization)
    assert np.isfinite(z).all() and np.ptp(z) > 0
    assert meta['tips'] == 2 and meta['forks'] == 0
    assert meta['anchor_count'] == (int(anchor_mask.sum()) if use_mask else len(Y))
    assert meta['unlabeled_count'] == (int((~anchor_mask).sum()) if use_mask else 0)
    assert meta['initialization'] == initialization


def test_unoriented_scfates_smoke_when_installed():
    pytest.importorskip('anndata')
    pytest.importorskip('scFates')
    t = np.linspace(0, 1, 60)
    Y = np.column_stack([t, t**2, np.sin(t * np.pi)])
    z, meta = fit_unoriented_scfates_axis(Y, nodes=10, seed=2)
    assert np.isfinite(z).all() and np.ptp(z) == pytest.approx(1.)
    assert meta['tips'] == 2 and meta['forks'] == 0
    assert meta['root_tip'] >= 0
    assert 'arbitrary orientation' in meta['root_rule']



def test_path_guard_rejects_disconnected_cycle_with_two_other_tips():
    from types import SimpleNamespace
    from pseudospace.trajectory_benchmark import _scfates_graph_guard
    # A separate 3-cycle plus a two-node path has two tips but is not one path.
    adjacency=np.array([[0,1,1,0,0],[1,0,1,0,0],[1,1,0,0,0],
                        [0,0,0,0,1],[0,0,0,1,0]],float)
    obj=SimpleNamespace(uns={'graph':{'B':adjacency}})
    with pytest.raises(ValueError,match='unbranched'):
        _scfates_graph_guard(obj)
