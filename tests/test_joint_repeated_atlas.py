"""Synthetic checks for the joint latent repeated-specimen atlas."""
import numpy as np
import pytest

pytest.importorskip('scipy')
pytest.importorskip('sklearn')

from pseudospace.joint_repeated_atlas import (
    fit_balanced_representation, fit_joint_atlas, project_joint_atlas,
    transform_balanced_representation,
)


def _synthetic(seed=5, n_per=90):
    rng = np.random.default_rng(seed)
    z0 = np.linspace(.01, .99, n_per)
    z = np.tile(z0, 2)
    specimens = np.repeat(['A', 'B'], n_per)
    mean = np.column_stack([1.8*z, np.sin(np.pi*z), np.cos(2*np.pi*z)])
    offsets = np.where((specimens == 'A')[:, None], [0., -.22, .12], [0., .22, -.12])
    Y = mean + offsets + rng.normal(0, .08, mean.shape)
    anatomy = np.where(z < .2, 0, np.where(z < .8, 1, 2))
    return Y, specimens, anatomy, z


def test_joint_atlas_recovers_order_learns_cutpoints_and_fixes_offset_gauge():
    Y, specimens, anatomy, z = _synthetic()
    atlas = fit_joint_atlas(Y, specimens, anatomy, z, max_iter=35)
    from scipy.stats import spearmanr
    assert spearmanr(atlas['z_mean'], z).statistic > .9
    assert atlas['cutpoints'][0] < atlas['cutpoints'][1]
    assert not np.allclose(atlas['cutpoints'], [.33, .67], atol=.04)
    np.testing.assert_allclose(np.mean(list(atlas['offsets'].values()), axis=0), 0., atol=1e-10)
    objective = np.array([record['objective'] for record in atlas['history']])
    assert np.isfinite(objective).all()
    assert np.all(np.diff(objective) <= 1e-7)
    assert all({'iteration', 'objective', 'sigma2', 'objective_change'} <= set(row)
               for row in atlas['history'])


def test_expression_only_projection_has_full_support_and_normalized_posterior():
    Y, specimens, anatomy, z = _synthetic(seed=9, n_per=60)
    atlas = fit_joint_atlas(Y, specimens, anatomy, z, max_iter=20)
    held = Y[::11] + np.array([.03, -.02, .01])
    projected = project_joint_atlas(atlas, held, block_size=3)
    p = projected['posterior']
    assert p.shape == (len(held), len(atlas['grid']))
    assert np.all(p > 0)
    np.testing.assert_allclose(p.sum(axis=1), 1.)
    for key in ('z_mean', 'z_MAP', 'entropy', 'residual_norm'):
        assert np.isfinite(projected[key]).all()


def test_balanced_weights_are_equal_by_specimen_even_after_unequal_duplication():
    Y, specimens, anatomy, z = _synthetic(n_per=40)
    # Duplicate only specimen A observations; each specimen retains half the total mass.
    dup = np.flatnonzero(specimens == 'A')
    Y2 = np.vstack([Y, Y[dup]])
    s2 = np.r_[specimens, specimens[dup]]
    a2 = np.r_[anatomy, anatomy[dup]]
    z2 = np.r_[z, z[dup]]
    atlas = fit_joint_atlas(Y2, s2, a2, z2, max_iter=5)
    masses = [atlas['weights'][s2 == sample].sum() for sample in ['A', 'B']]
    np.testing.assert_allclose(masses, [.5, .5])


def test_balanced_representation_uses_frozen_training_transform():
    rng = np.random.default_rng(17)
    X = rng.normal(size=(45, 8))
    specimens = np.repeat(['small', 'large'], [15, 30])
    Y, transform = fit_balanced_representation(X, specimens, n_components=5)
    held = rng.normal(size=(4, 8))
    projected = transform_balanced_representation(held, transform)
    np.testing.assert_allclose(Y, transform_balanced_representation(X, transform))
    assert projected.shape == (4, 5)
    with pytest.raises(ValueError, match='feature count'):
        transform_balanced_representation(np.zeros((2, 7)), transform)


def test_invalid_anatomy_and_initial_coordinates_are_rejected():
    Y, specimens, anatomy, z = _synthetic(n_per=20)
    with pytest.raises(ValueError, match='anatomy labels'):
        fit_joint_atlas(Y, specimens, np.full(len(Y), 4), z, max_iter=2)
    with pytest.raises(ValueError, match='initial_z'):
        fit_joint_atlas(Y, specimens, anatomy, np.full(len(Y), 1.5), max_iter=2)
    with pytest.raises(ValueError, match='one label per row'):
        fit_joint_atlas(Y, specimens, anatomy[:, None], z, max_iter=2)


def test_unlabeled_anatomy_is_opt_in_and_ignored_in_objective():
    from scipy.special import expit, logsumexp
    from pseudospace.joint_repeated_atlas import _basis

    Y, specimens, anatomy, z = _synthetic(seed=31, n_per=24)
    labels = anatomy.copy()
    labels[::3] = -1
    with pytest.raises(ValueError, match='anatomy labels'):
        fit_joint_atlas(Y, specimens, labels, z, max_iter=3)
    atlas = fit_joint_atlas(Y, specimens, labels, z, max_iter=4,
                            allow_unlabeled_anatomy=True)
    assert np.isfinite(atlas['posterior']).all()
    np.testing.assert_allclose(atlas['posterior'].sum(axis=1), 1.)
    assert np.isfinite([r['objective'] for r in atlas['history']]).all()

    # Reconstruct the recorded observed objective. Unlabeled rows contribute
    # only their expression likelihood, while known rows also receive ordinal
    # log probability; both groups retain the normalizing constants/penalties.
    grid, mean = atlas['grid'], atlas['mean']
    samples = atlas['specimens']
    offsets = np.vstack([atlas['offsets'][sample] for sample in samples])
    sample_index = np.searchsorted(samples, specimens.astype(str))
    residual = Y[:, None, :] - mean[None, :, :] - offsets[sample_index, None, :]
    sigma2, d = atlas['sigma2'], Y.shape[1]
    logp = -np.sum(residual ** 2, axis=2) / (2 * sigma2) - .5 * d * np.log(2 * np.pi * sigma2)
    c1, c2 = atlas['cutpoints']
    p1, p2 = expit((grid-c1)/.1), expit((grid-c2)/.1)
    ordinal = np.column_stack([1-p1, p1-p2, p2])
    known = labels >= 0
    logp[known] += np.log(np.maximum(ordinal[:, labels[known]].T, 1e-12))
    logp -= np.log(len(grid))
    weights = atlas['weights']
    objective = -float(weights @ logsumexp(logp, axis=1))
    B = _basis(grid, grid)
    d2 = np.diff(np.eye(B.shape[1]), n=2, axis=0)
    objective += .5 * .001 * np.sum((d2 @ atlas['beta']) ** 2)
    objective += .5 * 4. * np.mean(np.sum(offsets ** 2, axis=1))
    assert atlas['history'][-1]['objective'] == pytest.approx(objective, rel=1e-10, abs=1e-10)


def test_all_unknown_anatomy_has_no_ordinal_effect_or_cutpoint_fit():
    Y, specimens, _, z = _synthetic(seed=22, n_per=18)
    unknown = np.full(len(Y), -1)
    zero = fit_joint_atlas(Y, specimens, unknown, z, max_iter=3,
                           anatomy_strength=0., allow_unlabeled_anatomy=True)
    positive = fit_joint_atlas(Y, specimens, unknown, z, max_iter=3,
                               anatomy_strength=2., allow_unlabeled_anatomy=True)
    np.testing.assert_array_equal(zero['cutpoints'], [.33, .67])
    np.testing.assert_array_equal(positive['cutpoints'], [.33, .67])
    np.testing.assert_allclose(zero['posterior'], positive['posterior'])
    np.testing.assert_allclose([x['objective'] for x in zero['history']],
                               [x['objective'] for x in positive['history']])


@pytest.mark.parametrize('bad', [-2, 3, 1.5])
def test_opt_in_still_rejects_invalid_anatomy_labels(bad):
    Y, specimens, anatomy, z = _synthetic(n_per=20)
    anatomy = anatomy.astype(float)
    anatomy[0] = bad
    with pytest.raises(ValueError, match='anatomy labels'):
        fit_joint_atlas(Y, specimens, anatomy, z, max_iter=2,
                        allow_unlabeled_anatomy=True)


def test_iteration_limit_is_reported_as_not_converged():
    Y, specimens, anatomy, z = _synthetic(seed=4, n_per=35)
    atlas = fit_joint_atlas(Y, specimens, anatomy, z, max_iter=1)
    assert len(atlas['history']) == 1
    assert atlas['converged'] is False


@pytest.mark.parametrize('scale', [.2, 7.])
def test_prior_scale_makes_joint_fit_unit_invariant_and_restores_raw_units(scale):
    Y, specimens, anatomy, z = _synthetic(seed=43, n_per=24)
    params = dict(max_iter=12, tol=1e-7)
    raw = fit_joint_atlas(Y, specimens, anatomy, z, **params)
    scaled = fit_joint_atlas(Y * scale, specimens, anatomy, z,
                             prior_scale=scale, **params)

    assert raw['converged'] == scaled['converged']
    assert len(raw['history']) == len(scaled['history'])
    for field in ('objective', 'objective_change'):
        np.testing.assert_allclose([item[field] for item in raw['history']],
                                   [item[field] for item in scaled['history']],
                                   rtol=1e-9, atol=1e-9, equal_nan=True)
    for item, scaled_item in zip(raw['history'], scaled['history']):
        assert scaled_item['sigma2'] == pytest.approx(item['sigma2'] * scale ** 2,
                                                       rel=1e-9, abs=1e-12)
    np.testing.assert_allclose(raw['posterior'], scaled['posterior'], rtol=1e-6, atol=5e-9)
    np.testing.assert_allclose(raw['z_mean'], scaled['z_mean'], rtol=1e-7, atol=1e-8)
    np.testing.assert_array_equal(raw['z_MAP'], scaled['z_MAP'])
    np.testing.assert_allclose(scaled['mean'], raw['mean'] * scale, rtol=1e-7, atol=1e-8)
    np.testing.assert_allclose(scaled['beta'], raw['beta'] * scale, rtol=1e-7, atol=1e-8)
    np.testing.assert_allclose(scaled['sigma2'], raw['sigma2'] * scale ** 2,
                               rtol=1e-8, atol=1e-12)
    for sample in raw['specimens']:
        np.testing.assert_allclose(scaled['offsets'][sample], raw['offsets'][sample] * scale,
                                   rtol=1e-7, atol=1e-8)

    held = Y[::7] + [.03, -.02, .01]
    offset = raw['offsets']['B']
    p_raw = project_joint_atlas(raw, held, offset=offset)
    p_scaled = project_joint_atlas(scaled, held * scale, offset=offset * scale)
    np.testing.assert_allclose(p_scaled['posterior'], p_raw['posterior'], rtol=1e-6, atol=5e-9)
    np.testing.assert_allclose(p_scaled['z_mean'], p_raw['z_mean'], rtol=1e-7, atol=1e-8)
    np.testing.assert_array_equal(p_scaled['z_MAP'], p_raw['z_MAP'])


def test_default_prior_scale_matches_explicit_one_exactly():
    Y, specimens, anatomy, z = _synthetic(seed=44, n_per=20)
    implicit = fit_joint_atlas(Y, specimens, anatomy, z, max_iter=5)
    explicit = fit_joint_atlas(Y, specimens, anatomy, z, max_iter=5, prior_scale=1.)
    for key in ('mean', 'beta', 'sigma2', 'posterior', 'z_mean', 'z_MAP', 'cutpoints'):
        np.testing.assert_array_equal(implicit[key], explicit[key])
    assert implicit['history'] == explicit['history']
    assert implicit['prior_scale'] == explicit['prior_scale'] == 1.
    for sample in implicit['specimens']:
        np.testing.assert_array_equal(implicit['offsets'][sample], explicit['offsets'][sample])


@pytest.mark.parametrize('scale', [0., -1., np.nan, np.inf, True, False,
                                   [1.], np.array([1.]), 1e-200, 1e200])
def test_invalid_prior_scale_is_rejected(scale):
    Y, specimens, anatomy, z = _synthetic(seed=45, n_per=12)
    with pytest.raises(ValueError, match='prior_scale'):
        fit_joint_atlas(Y, specimens, anatomy, z, max_iter=2, prior_scale=scale)


def test_prior_scale_rejects_underflowed_restored_variance():
    with pytest.raises(ValueError, match='prior_scale restoration'):
        fit_joint_atlas(np.zeros((6, 2)), np.repeat(['A', 'B'], 3),
                        np.tile([0, 1, 2], 2), np.tile([.1, .5, .9], 2),
                        max_iter=2, prior_scale=1e-160)
