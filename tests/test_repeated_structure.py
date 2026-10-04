import numpy as np
import pytest

from pseudospace.repeated_structure import (
    evaluate_prediction,
    fit_atlas,
    map_atlas,
    refine_coordinate,
    select_panel,
)


def toy():
    rng = np.random.default_rng(7)
    z = np.tile(np.linspace(0.03, 0.97, 40), 3)
    specimen = np.repeat(["a", "b", "c"], 40)
    x = np.column_stack([z + rng.normal(0, .03, len(z)), np.sin(4*z) + rng.normal(0, .04, len(z)), rng.normal(size=len(z))])
    anatomy = np.minimum((z * 3).astype(int), 2)
    return x, z, specimen, anatomy


def test_atlas_scaling_is_frozen_from_training_and_anatomy_bounds_mapping():
    x, z, s, a = toy()
    atlas = fit_atlas(x[:80], z[:80], s[:80])
    frozen_mean = atlas.mean.copy()
    held = x[80:]
    held[0] += 1000
    mapped = map_atlas(held, atlas, anatomy=a[80:], initial=z[80:], radius=.2)
    assert np.array_equal(atlas.mean, frozen_mean)
    assert np.all(mapped >= a[80:] / 3)
    assert np.all(mapped <= (a[80:] + 1) / 3 + 1e-12)


def test_gene_panel_selection_is_reproducible_and_prior_is_optional():
    x, z, s, a = toy()
    i1, score1 = select_panel(x, z, s, a, 2)
    i2, score2 = select_panel(x, z, s, a, 2)
    assert np.array_equal(i1, i2)
    assert np.allclose(score1, score2)
    assert 0 in i1
    _, prior_score = select_panel(x, z, s, a, 3, prior=[0, 0, 1])
    assert prior_score[2] > score1[2]


def test_refinement_returns_cross_fitted_updates_and_history():
    x, z, s, a = toy()
    updated, history = refine_coordinate(x, z, s, a, iterations=1, radius=.08)
    assert len(history) == 2
    assert np.array_equal(history[0], z)
    assert updated.shape == z.shape
    assert np.all(np.abs(updated - z) <= .08 + 1e-12)
    assert np.all(updated >= a / 3)
    assert np.all(updated <= (a + 1) / 3 + 1e-12)


def test_heldout_prediction_reports_baseline_comparison():
    x, z, s, _ = toy()
    result = evaluate_prediction(x[:80], z[:80], s[:80], x[80:], z[80:], ["gradient", "peak", "noise"])
    assert result["gene_metrics"]["gene"].tolist() == ["gradient", "peak", "noise"]
    assert result["mean_mse"] < result["mean_baseline_mse"]




def test_invalid_inputs_fail_early():
    x, z, s, a = toy()
    single_specimen_atlas = fit_atlas(x[:40], z[:40], np.repeat("one", 40))
    assert single_specimen_atlas.beta.shape[0] == x.shape[1]
    with pytest.raises(ValueError):
        map_atlas(x, fit_atlas(x, z, s), anatomy=np.full(len(z), 3))
    with pytest.raises(ValueError):
        map_atlas(x, fit_atlas(x, z, s), anatomy=np.full(len(z), .5))


def test_atlas_preserves_constant_gene_and_constant_panel_scores_tie():
    x, z, s, a = toy()
    x = np.column_stack([x, np.full(len(x), 7.0)])
    atlas = fit_atlas(x, z, s)
    assert np.allclose(atlas.predict(np.linspace(0, 1, 21))[:, 3], 7.0, atol=1e-5)
    constant = np.full((120, 3), 4.0)
    zc = np.tile(np.linspace(0, 1, 40), 3)
    sc = np.repeat(["a", "b", "c"], 40)
    ac = np.minimum((zc * 3).astype(int), 2)
    _, scores = select_panel(constant, zc, sc, ac, 3)
    assert np.allclose(scores, scores[0])


def test_panel_signed_agreement_downgrades_reversed_replicate_curve():
    z = np.tile(np.linspace(.01, .99, 120), 2)
    s = np.repeat(["a", "b"], 120)
    anatomy = np.minimum((z * 3).astype(int), 2)
    good = z + np.random.default_rng(4).normal(0, .01, len(z))
    reversed_in_b = good.copy()
    reversed_in_b[120:] = 1 - z[120:]
    x = np.column_stack([good, reversed_in_b])
    _, score = select_panel(x, z, s, anatomy, 2)
    assert score[0] > score[1]






def test_mapping_rejects_anatomy_interval_absent_from_custom_grid():
    x, z, s, _ = toy()
    grid = np.array([0, .01, .02, .03, .04, .05, .06, 1.0])
    atlas = fit_atlas(x, z, s, grid=grid)
    with pytest.raises(ValueError, match="no points in anatomy interval"):
        map_atlas(x[:1], atlas, anatomy=[1])


def test_explicit_original_specimen_weights_survive_duplicate_rows():
    z = np.tile(np.linspace(.05, .95, 30), 2)
    specimen = np.repeat(["a", "b"], 30)
    x = np.r_[np.zeros(30), np.ones(30)][:, None]
    default = fit_atlas(x, z, specimen)
    weighted = fit_atlas(x, z, specimen, specimen_weights={"a": 1, "b": 2})
    assert default.mean[0] == pytest.approx(.5)
    assert weighted.mean[0] == pytest.approx(2 / 3)
    with pytest.raises(ValueError):
        fit_atlas(x, z, specimen, specimen_weights={"a": 0, "b": 1})
    with pytest.raises(ValueError):
        fit_atlas(x, z, specimen, specimen_weights={"a": 1})


def test_spline_distribution_prediction_averages_nonlinear_decoder_values():
    z = np.tile(np.linspace(0., 1., 100), 2)
    specimen = np.repeat(['a', 'b'], len(z) // 2)
    x = np.column_stack([z ** 2, np.sin(2 * np.pi * z)])
    atlas = fit_atlas(x, z, specimen, ridge=.001)
    support = np.array([.1, .9])
    weights = np.array([[.25, .75], [.6, .4]])

    actual = atlas.predict_distribution(support, weights)
    expected = weights @ atlas.predict(support)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
    assert not np.allclose(actual[0], atlas.predict([weights[0] @ support])[0], atol=1e-3)


def test_spline_distribution_prediction_supports_row_specific_grids_and_empty_queries():
    x, z, specimen, _ = toy()
    atlas = fit_atlas(x, z, specimen)
    support = np.array([[.05, .4, .8], [.2, .6, .95]])
    weights = np.array([[.2, .5, .3], [.6, .1, .3]])
    actual = atlas.predict_distribution(support, weights)
    expected = np.vstack([w @ atlas.predict(row) for row, w in zip(support, weights)])
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
    assert atlas.predict_distribution(np.array([.1, .9]), np.empty((0, 2))).shape == (0, x.shape[1])
    assert atlas.predict_distribution(np.empty((0, 3)), np.empty((0, 3))).shape == (0, x.shape[1])


@pytest.mark.parametrize('support,weights', [
    ([.1, .9], [[.5]]),
    ([[.1, .9]], [[.5, .4]]),
    ([.1, 1.1], [[.5, .5]]),
    ([[.1, np.nan]], [[.5, .5]]),
    ([.1, .9], [[.5, -.5]]),
    ([.1, .9], [[np.nan, np.nan]]),
    ([.1, .9], [[.2, .2]]),
    ([.1, .9], [0., 1.]),
    ([], np.empty((1, 0))),
])
def test_spline_distribution_prediction_rejects_invalid_support_or_weights(support, weights):
    x, z, specimen, _ = toy()
    atlas = fit_atlas(x, z, specimen)
    with pytest.raises(ValueError):
        atlas.predict_distribution(support, weights)
