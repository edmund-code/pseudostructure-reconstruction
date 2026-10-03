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


def test_iteration_limit_is_reported_as_not_converged():
    Y, specimens, anatomy, z = _synthetic(seed=4, n_per=35)
    atlas = fit_joint_atlas(Y, specimens, anatomy, z, max_iter=1)
    assert len(atlas['history']) == 1
    assert atlas['converged'] is False
