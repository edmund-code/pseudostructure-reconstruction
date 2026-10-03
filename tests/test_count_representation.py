"""Synthetic checks for frozen count-exposure residual representations."""
import numpy as np
import pytest
pytest.importorskip('scipy')
pytest.importorskip('sklearn')
from scipy.stats import spearmanr

from pseudospace.count_representation import (
    fit_count_representation, transform_count_representation,
)


def test_neutral_poisson_profiles_have_small_library_association_after_residualization():
    rng = np.random.default_rng(61)
    n, g = 800, 60
    library = np.exp(rng.uniform(np.log(100), np.log(1000), n))
    probability = rng.uniform(.001, .0025, g)
    counts = rng.poisson(library[:, None] * probability[None, :])
    specimens = np.repeat(['A', 'B', 'C', 'D'], n // 4)
    embedding, _ = fit_count_representation(counts, library, specimens,
                                            n_components=8, theta=100.)
    corr = spearmanr(embedding[:, 0], library).statistic
    low, high = library < np.median(library), library >= np.median(library)
    assert abs(corr) < .15
    assert np.linalg.norm(embedding[low].mean(axis=0) - embedding[high].mean(axis=0)) < .3


def test_gene_probabilities_balance_specimens_and_ignore_specimen_duplication():
    counts = np.vstack([np.tile([80, 10], (4, 1)), np.tile([20, 10], (40, 1))])
    library = np.full(len(counts), 100.)
    specimens = np.array(['A']*4 + ['B']*40)
    _, transform = fit_count_representation(counts, library, specimens, n_components=1)
    np.testing.assert_allclose(transform['probability'], [.5, .1])
    duplicate = np.flatnonzero(specimens == 'A')
    _, transform_dup = fit_count_representation(
        np.vstack([counts, counts[duplicate]]), np.r_[library, library[duplicate]],
        np.r_[specimens, specimens[duplicate]], n_components=1)
    np.testing.assert_allclose(transform_dup['probability'], transform['probability'])


def test_query_projection_uses_frozen_probabilities_and_zero_training_gene_is_uninformative():
    rng = np.random.default_rng(7)
    library = np.full(80, 1000.)
    gene0 = rng.poisson(120, 80)
    training_counts = np.column_stack([gene0, np.zeros(80, dtype=int)])
    specimens = np.repeat(['A', 'B'], 40)
    _, transform = fit_count_representation(training_counts, library, specimens,
                                            n_components=1)
    query_a = np.array([[125, 0], [110, 0]])
    query_b = np.array([[125, 500], [110, 300]])
    projected_a = transform_count_representation(query_a, [1000, 1000], transform)
    projected_b = transform_count_representation(query_b, [1000, 1000], transform)
    np.testing.assert_allclose(projected_a, projected_b)
    np.testing.assert_array_equal(transform['probability'][1:], [0.])
    original_probability = transform['probability'].copy()
    _ = transform_count_representation(np.array([[100, 0]]), [5000], transform)
    np.testing.assert_array_equal(transform['probability'], original_probability)


def test_invalid_counts_exposures_and_component_counts_are_rejected():
    counts = np.array([[2, 1], [1, 2], [2, 2]])
    library = np.array([10., 10., 10.])
    specimens = ['a', 'a', 'b']
    with pytest.raises(ValueError, match='integer-valued'):
        fit_count_representation(counts + .1, library, specimens, n_components=1)
    with pytest.raises(ValueError, match='at least'):
        fit_count_representation(counts, [2., 10., 10.], specimens, n_components=1)
    with pytest.raises(ValueError, match='positive integer'):
        fit_count_representation(counts, library, specimens, n_components=1.5)
