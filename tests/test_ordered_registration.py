from itertools import combinations_with_replacement

import numpy as np
import pytest

from pseudospace.ordered_registration import calibration_support, local_curve_means, register_ordered_means


def test_partial_reverse_path_does_not_stretch_to_reference_endpoints():
    reference = np.arange(5.)[:, None]
    result = register_ordered_means(np.array([[3.], [2.], [1.]]), reference, np.linspace(0, 1, 5))
    np.testing.assert_array_equal(result['z'], [.75, .5, .25])
    assert result['orientation'] == -1 and result['cost'] == 0
    assert result['reference_span'] == .5


def test_exact_ordered_cost_matches_exhaustive_paths_in_both_orientations():
    reference = np.array([[0., 0.], [1., 2.], [2., 1.], [3., 3.]])
    query = np.array([[2.6, 2.8], [.8, 1.9], [2.1, 1.1]])
    cost = np.mean((query[:, None] - reference[None]) ** 2, axis=2)
    paths = list(combinations_with_replacement(range(len(reference)), len(query)))
    possible = [float(cost[np.arange(len(query)), path].mean())
                for p in paths for path in [p, p[::-1]]]
    result = register_ordered_means(query, reference, np.linspace(0, 1, len(reference)))
    assert result['cost'] == pytest.approx(min(possible))
    assert result['cost'] >= result['independent_cost']
    assert np.all(np.diff(result['reference_index']) * result['orientation'] >= 0)


def test_flat_path_ties_are_deterministic_and_collapse_is_reported():
    result = register_ordered_means(np.ones((3, 2)), np.ones((4, 2)), np.linspace(0, 1, 4))
    np.testing.assert_array_equal(result['z'], np.zeros(3))
    assert result['orientation'] == 1 and result['orientation_gap'] == 0
    assert result['plateau_fraction'] == 1 and result['reference_span'] == 0


def test_local_means_use_normalized_weights_and_report_effective_support():
    Y = np.array([[0., 1.], [1., 2.], [2., 3.]])
    z = np.array([0., .5, 1.])
    means, support = local_curve_means(Y, z, z, bandwidth=.1)
    weights = np.exp(-.5 * ((z[:, None] - z) / .1) ** 2)
    weights /= weights.sum(axis=1)[:, None]
    np.testing.assert_allclose(means, weights @ Y)
    np.testing.assert_allclose(support, 1 / (weights ** 2).sum(axis=1))
    assert np.all((support >= 1) & (support <= len(Y)))


@pytest.mark.parametrize('grid', [[0., 0., 1.], [-.1, .5, 1.], [0., np.nan, 1.], [0., 1.]])
def test_registration_rejects_invalid_reference_grid(grid):
    with pytest.raises(ValueError, match='grid'):
        register_ordered_means(np.ones((3, 1)), np.ones((3, 1)), grid)


@pytest.mark.parametrize('query', [np.empty((0, 1)), np.ones(3), [[np.inf]]])
def test_registration_rejects_invalid_means(query):
    with pytest.raises(ValueError):
        register_ordered_means(query, np.ones((3, 1)), np.linspace(0, 1, 3))


def test_registration_rejects_feature_mismatch():
    with pytest.raises(ValueError, match='features'):
        register_ordered_means(np.ones((3, 2)), np.ones((3, 1)), np.linspace(0, 1, 3))


@pytest.mark.parametrize('bandwidth', [0., -1., np.nan, np.inf, True])
def test_local_means_rejects_invalid_bandwidth(bandwidth):
    with pytest.raises(ValueError, match='bandwidth'):
        local_curve_means(np.ones((3, 1)), np.linspace(0, 1, 3), np.linspace(0, 1, 3), bandwidth=bandwidth)


def test_support_radius_excludes_self_and_preserves_uncovered_reference_region():
    result = calibration_support(np.array([[0.], [1.], [2.]]),
                                 np.array([[.5], [10.]]), k=1)
    assert result['radius'] == 1.
    np.testing.assert_array_equal(result['supported'], [True, False])
    empty = calibration_support(np.array([[0.], [1.]]), np.empty((0, 1)), k=1)
    assert empty['supported'].shape == (0,)


@pytest.mark.parametrize('k', [0, True, 3, 1.5])
def test_support_rejects_invalid_neighbor_count(k):
    with pytest.raises(ValueError, match='k must'):
        calibration_support(np.ones((3, 1)), np.ones((2, 1)), k=k)
