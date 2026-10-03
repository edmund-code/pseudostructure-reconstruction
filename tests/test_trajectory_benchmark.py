import numpy as np
import pytest
pytest.importorskip("scipy")
pytest.importorskip("sklearn")

from pseudospace.trajectory_benchmark import (
    fit_dpt_axis, fit_scfates_axis, project_neighbor_axis,
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


def test_scanpy_dpt_smoke_when_installed():
    pytest.importorskip('anndata')
    pytest.importorskip('scanpy')
    t = np.linspace(0, 1, 60)
    Y = np.column_stack([t, t**2, np.sin(t * np.pi)])
    anatomy = np.where(t < 1/3, 0, np.where(t < 2/3, 1, 2))
    z, meta = fit_dpt_axis(Y, anatomy, n_neighbors=10, seed=2)
    assert np.isfinite(z).all() and np.ptp(z) > 0
    assert meta['root_index'] >= 0


def test_scfates_smoke_when_installed():
    pytest.importorskip('anndata')
    pytest.importorskip('scFates')
    t = np.linspace(0, 1, 60)
    Y = np.column_stack([t, t**2, np.sin(t * np.pi)])
    anatomy = np.where(t < 1/3, 0, np.where(t < 2/3, 1, 2))
    z, meta = fit_scfates_axis(Y, anatomy, nodes=10, seed=2)
    assert np.isfinite(z).all() and np.ptp(z) > 0
    assert meta['tips'] == 2 and meta['forks'] == 0



def test_path_guard_rejects_disconnected_cycle_with_two_other_tips():
    from types import SimpleNamespace
    from pseudospace.trajectory_benchmark import _scfates_graph_guard
    # A separate 3-cycle plus a two-node path has two tips but is not one path.
    adjacency=np.array([[0,1,1,0,0],[1,0,1,0,0],[1,1,0,0,0],
                        [0,0,0,0,1],[0,0,0,1,0]],float)
    obj=SimpleNamespace(uns={'graph':{'B':adjacency}})
    with pytest.raises(ValueError,match='unbranched'):
        _scfates_graph_guard(obj)
