"""Synthetic checks for specimen offset fitting and frozen query projection."""
import numpy as np
import pytest

pytest.importorskip('scipy')
pytest.importorskip('sklearn')

from pseudospace.specimen_registration import fit_specimen_offset
from pseudospace.trajectory_benchmark import project_neighbor_axis


def _atlas(mean):
    return {"grid": np.linspace(0., 1., len(mean)), "mean": np.asarray(mean, dtype=float)}


def test_restricted_offset_is_orthogonal_to_retained_curve_space():
    x = np.linspace(0., 1., 101)
    mean = np.column_stack([x, np.zeros_like(x), np.zeros_like(x)])
    y = mean[::5] + np.array([0.7, -0.5, 0.2])
    offset, meta = fit_specimen_offset(_atlas(mean), y, penalty=0)

    assert meta["rank"] == 1
    np.testing.assert_allclose(offset, [0., -0.5, 0.2], atol=1e-10)
    assert meta["retained_fraction"] == pytest.approx(1.)
    assert meta["calibration_n"] == len(y)
    assert len(meta["history"]) == 5


def test_unrestricted_absorbs_tangent_while_restricted_does_not():
    x = np.linspace(0., 1., 101)
    mean = np.column_stack([x, np.zeros_like(x)])
    y = mean[::5] + [0.5, 0.25]
    restricted, _ = fit_specimen_offset(_atlas(mean), y, penalty=0)
    unrestricted, _ = fit_specimen_offset(_atlas(mean), y, restricted=False, penalty=0)
    np.testing.assert_allclose(restricted, [0., .25], atol=1e-10)
    assert unrestricted[0] > .25
    assert unrestricted[1] == pytest.approx(.25)


def test_full_rank_restriction_returns_exact_zero():
    t = np.linspace(-1., 1., 41)
    rotation = np.array([[np.cos(.37), -np.sin(.37)],
                         [np.sin(.37), np.cos(.37)]])
    mean = np.column_stack([t, t ** 2]) @ rotation
    offset, meta = fit_specimen_offset(_atlas(mean), mean[::4] + [.4, -.2],
                                       variance_fraction=1.)
    assert meta["rank"] == 2
    np.testing.assert_array_equal(offset, np.zeros(2))
    assert meta["retained_fraction"] == 1.


def test_projection_uses_frozen_training_population_after_offset_fit():
    z = np.linspace(0., 1., 101)
    train = np.column_stack([z, z ** 2])
    atlas = _atlas(train)
    calibration = train[::10] + [0., .3]
    offset, _ = fit_specimen_offset(atlas, calibration, restricted=False, penalty=0)
    query = train[[20, 70]] + [0., .3]
    projected = project_neighbor_axis(train, z, query - offset, k=3)
    projected_with_unrelated_query = project_neighbor_axis(
        train, z, np.vstack([query - offset, [[100., -100.]]]), k=3
    )
    np.testing.assert_allclose(projected_with_unrelated_query[:2], projected)
    assert np.isfinite(projected).all()


def test_invalid_dimensions_and_parameters_are_rejected():
    atlas = _atlas(np.column_stack([np.linspace(0., 1., 5), np.zeros(5)]))
    with pytest.raises(ValueError, match="feature count"):
        fit_specimen_offset(atlas, np.ones((3, 1)))
    with pytest.raises(ValueError, match="finite"):
        fit_specimen_offset(atlas, [[0., np.nan]])
    for kwargs in ({"variance_fraction": 0}, {"variance_fraction": 1.1},
                   {"penalty": -1}, {"penalty": np.inf}, {"restricted": 1},
                   {"iterations": 0}):
        with pytest.raises(ValueError):
            fit_specimen_offset(atlas, np.ones((3, 2)), **kwargs)
    with pytest.raises(ValueError, match="grid"):
        fit_specimen_offset({"grid": [0., 1.], "mean": np.ones((3, 2))}, np.ones((2, 2)))
    for grid in ([.1, 1.], [0., .4, .3, .7, 1.]):
        malformed = {"grid": grid, "mean": np.ones((len(grid), 2))}
        with pytest.raises(ValueError, match="grid"):
            fit_specimen_offset(malformed, np.ones((2, 2)))
