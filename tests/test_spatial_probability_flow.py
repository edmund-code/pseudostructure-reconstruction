"""Synthetic checks for the finite-step spatial probability-flow pilot."""
import numpy as np
import pytest
import pseudospace.spatial_probability_flow as spf

pytest.importorskip('scipy')
pytest.importorskip('sklearn')

from pseudospace.spatial_probability_flow import (
    advect, energy_distance, fit_density_atlas, fit_flow, posterior,
    scaffold_centroid_rows, window_clouds,
)


def test_flow_learns_translation_and_predicts_independent_cloud():
    rng = np.random.default_rng(8)
    base = rng.normal(size=(40, 2))
    shift = np.array([.8, -.35])
    clouds = [base + j * shift for j in range(3)]
    model = fit_flow(clouds, [0., .5, 1.], ridge=.01, epsilon=.5)
    independent = rng.normal(size=(35, 2)) + np.array([.2, .1])
    target = independent + shift
    predicted = advect(model, independent, 0., .5, steps=12)
    assert energy_distance(predicted, target) < energy_distance(independent, target)
    assert np.dot(np.mean(predicted - independent, axis=0), shift) > 0


def test_representation_freezes_train_transform_and_accepts_empty_test():
    rng = np.random.default_rng(31)
    train = rng.normal(size=(24, 5))
    test = rng.normal(size=(6, 5))
    emb_a, held_a, fit_a = spf.representation(train, test, n_components=4)
    emb_b, held_b, fit_b = spf.representation(train, test + 1e6, n_components=4)
    np.testing.assert_allclose(emb_a, emb_b)
    assert not np.allclose(held_a, held_b)
    np.testing.assert_allclose(fit_a['scaler'].mean_, fit_b['scaler'].mean_)
    _, empty, _ = spf.representation(train, np.empty((0, 5)), n_components=4)
    assert empty.shape == (0, 4)


def test_benchmark_local_regression_differs_from_centroid_on_curved_mean():
    rng = np.random.default_rng(44)
    ztrain = np.linspace(0, 1, 240)
    # The local conditional mean changes nonlinearly, while broad symmetric
    # window centroids average over a different span of the curve.
    ytrain = np.column_stack([ztrain**2, np.sin(2*np.pi*ztrain)])
    ytrain += rng.normal(0, .015, ytrain.shape)
    ztest = np.linspace(0, 1, 240)
    ytest = np.column_stack([ztest**2, np.sin(2*np.pi*ztest)])
    ytest += rng.normal(0, .015, ytest.shape)
    result = spf.benchmark(ytrain, ztrain, ytest, ztest, [.2, .5, .8], .16,
                           min_count=10, max_points=100, seed=4)
    scores = result.pivot(index='transition', columns='method', values='energy_distance')
    assert np.any(np.abs(scores['dpt_local'] - scores['centroid']) > 1e-5)


def test_benchmark_holdout_changes_scores_but_not_fitted_flow(monkeypatch):
    rng = np.random.default_rng(45)
    ztrain = np.linspace(0, 1, 180)
    ytrain = np.column_stack([ztrain, ztrain**2]) + rng.normal(0, .03, (180, 2))
    ztest = np.linspace(0, 1, 180)
    ytest = np.column_stack([ztest, ztest**2]) + rng.normal(0, .03, (180, 2))
    original = spf.fit_flow
    fitted = []

    def spy(clouds, centers, *args, **kwargs):
        model = original(clouds, centers, *args, **kwargs)
        fitted.append(model['coef'].copy())
        return model

    monkeypatch.setattr(spf, 'fit_flow', spy)
    kwargs = dict(centers=[.2, .5, .8], halfwidth=.15, min_count=10,
                  max_points=80, seed=12)
    a = spf.benchmark(ytrain, ztrain, ytest, ztest, **kwargs)
    b = spf.benchmark(ytrain, ztrain, ytest + 4, ztest, **kwargs)
    np.testing.assert_allclose(fitted[0], fitted[1])
    assert not np.allclose(a.energy_distance, b.energy_distance)


def test_fitted_flow_is_training_only_when_external_target_changes():
    rng = np.random.default_rng(2)
    source = rng.normal(size=(28, 2))
    train = [source, source + [0.4, 0.2], source + [0.8, 0.4]]
    model_a = fit_flow(train, [0., .5, 1.], epsilon=.4)
    # Evaluation targets are deliberately different; neither enters fit_flow.
    target_a = rng.normal(size=(20, 2))
    target_b = target_a + 50
    assert not np.allclose(target_a, target_b)
    model_b = fit_flow(train, [0., .5, 1.], epsilon=.4)
    np.testing.assert_allclose(model_a['coef'], model_b['coef'])


def test_sinkhorn_handles_separated_outlier_points_with_finite_barycenters():
    rng = np.random.default_rng(22)
    source = np.vstack([rng.normal(size=(30, 2)), [[24., 24.], [-24., 24.]]])
    target = source + np.array([.3, -.2])
    model = fit_flow([source, target], [0., .5], ridge=.1, epsilon=.5)
    moved = advect(model, source, 0., .5)
    assert moved.shape == source.shape
    assert np.isfinite(moved).all()


def test_posterior_is_normalized_and_unrelated_labels_do_not_affect_unlabeled_projection():
    rng = np.random.default_rng(11)
    z = np.linspace(0, 1, 120)
    Y = np.column_stack([z, np.sin(2*np.pi*z)]) + rng.normal(0, .035, (120, 2))
    labels = rng.choice(['S1', 'S2', 'S3'], size=len(z))  # unrelated to molecular order
    atlas = fit_density_atlas(Y, z, labels, np.linspace(.1, .9, 5), .14,
                              anatomy_strength=0., min_count=8)
    mean_z, probs, entropy = posterior(atlas, Y[:9])
    mean_z2, probs2, _ = posterior(atlas, Y[:9], labels=np.repeat('S1', 9))
    assert probs.shape == (9, 5)
    np.testing.assert_allclose(probs.sum(axis=1), 1.)
    np.testing.assert_allclose(probs, probs2)
    np.testing.assert_allclose(mean_z, mean_z2)
    assert np.isfinite(entropy).all()


def test_window_clouds_refuses_sparse_support_and_is_deterministic():
    Y = np.arange(80, dtype=float).reshape(40, 2)
    z = np.linspace(0, 1, 40)
    with pytest.raises(ValueError, match='requires'):
        window_clouds(Y, z, [.1, .9], .03, min_count=8)
    a = window_clouds(Y, z, [.4, .6], .3, min_count=8, max_points=12, seed=9)
    b = window_clouds(Y, z, [.4, .6], .3, min_count=8, max_points=12, seed=9)
    assert all(len(x) == 12 for x in a)
    for x, y in zip(a, b):
        np.testing.assert_array_equal(x, y)


def test_window_clouds_balances_specimens_and_refuses_sparse_specimen_window():
    z = np.full(90, .5)
    specimens = np.array(['A'] * 30 + ['B'] * 60)
    Y = np.r_[np.zeros((30, 1)), np.full((60, 1), 10.)]
    clouds = window_clouds(Y, z, [.4, .6], .2, min_count=10,
                           max_points=128, seed=3, specimens=specimens)
    assert [len(x) for x in clouds] == [60, 60]
    np.testing.assert_allclose([x.mean() for x in clouds], [5., 5.])
    sparse_specimens = np.array(['A'] * 5 + ['B'] * 20)
    with pytest.raises(ValueError, match='per-specimen'):
        window_clouds(np.zeros((25, 1)), np.full(25, .5), [.4, .6], .2,
                      min_count=6, specimens=sparse_specimens)


def test_endpoint_only_anatomy_does_not_use_s2_projection_label():
    rng = np.random.default_rng(17)
    z = np.linspace(0, 1, 150)
    Y = np.column_stack([z, z**2]) + rng.normal(0, .04, (150, 2))
    labels = np.where(z < .35, 'S1', np.where(z > .65, 'S3', 'S2'))
    atlas = fit_density_atlas(Y, z, labels, np.linspace(.1, .9, 5), .15,
                              endpoint_only=True, min_count=8)
    probe = Y[70:72]
    _, p_unlabeled, _ = posterior(atlas, probe)
    _, p_s2, _ = posterior(atlas, probe, labels=['S2', 'S2'])
    np.testing.assert_allclose(p_unlabeled, p_s2)


def test_specimen_balanced_latent_fit_returns_normalized_unlabeled_projection():
    rng = np.random.default_rng(28)
    z = np.tile(np.linspace(.02, .98, 100), 2)
    specimens = np.repeat(['A', 'B'], 100)
    Y = np.column_stack([z, np.sin(2*np.pi*z)])
    Y += np.where((specimens == 'B')[:, None], .05, -.05)
    Y += rng.normal(0, .025, Y.shape)
    labels = np.where(z < .33, 'S1', np.where(z < .67, 'S2', 'S3'))
    result = spf.latent_fit(Y, z, labels, np.linspace(.15, .85, 4), .18,
                            iterations=2, anatomy_strength=0., min_count=10,
                            specimens=specimens)
    projected_z, probabilities, entropy = posterior(result['atlas'], Y[:12])
    assert np.isfinite(projected_z).all() and np.isfinite(entropy).all()
    np.testing.assert_allclose(probabilities.sum(axis=1), 1.)


def test_scaffold_centroids_match_nearest_within_specimen_and_mark_absent():
    query_xy = np.array([[.2, .1], [8., 8.], [.1, .2]])
    query_samples = ['a', 'a', 'b']
    reference_xy = np.array([[0., 0.], [0., .3], [10., 10.]])
    reference_samples = ['a', 'b', 'b']
    rows = scaffold_centroid_rows(query_xy, query_samples, reference_xy,
                                  reference_samples, tolerance=.5)
    np.testing.assert_array_equal(rows, [0, -1, 1])


@pytest.mark.parametrize('query_xy, reference_xy, match', [
    ([[0., 0.]], [[-.1, 0.], [.1, 0.]], 'ambiguous'),
    ([[0., 0.], [.01, 0.]], [[0., 0.]], 'reuse'),
])
def test_scaffold_centroids_refuses_ambiguous_or_reused_matches(query_xy, reference_xy, match):
    with pytest.raises(ValueError, match=match):
        scaffold_centroid_rows(query_xy, ['s'] * len(query_xy), reference_xy,
                              ['s'] * len(reference_xy), tolerance=.5)


def test_scaffold_centroids_handles_zero_or_one_reference_row():
    empty = scaffold_centroid_rows([[0., 0.]], ['s'], np.empty((0, 2)), [], tolerance=1.)
    np.testing.assert_array_equal(empty, [-1])
    one = scaffold_centroid_rows([[.5, 0.]], ['s'], [[0., 0.]], ['s'], tolerance=1.)
    np.testing.assert_array_equal(one, [0])
    with pytest.raises(ValueError, match='tolerance'):
        scaffold_centroid_rows([[0., 0.]], ['s'], [[0., 0.]], ['s'], tolerance=0.)
