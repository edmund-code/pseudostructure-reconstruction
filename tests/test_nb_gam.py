import numpy as np
import pytest
from scipy import sparse

from pseudospace.nb_gam import (A_MIN, GENE_MODEL, _contrast_se, _orthobasis, _outer, nb_dispersion, nb_fit,
                                nb_group_difference_curves, nb_species_trajectories, nested_nb_lr)
from pseudospace.pathway_calibration import contrast_designs, nested_partial_f
from pseudospace.pathway_remodeling import fit_nested_trajectories
from pseudospace.stats_gam import gam_internal_knots, make_gam_design


def _cohort(seed, sizes=(60, 60, 60, 60)):
    rng = np.random.default_rng(seed)
    specimen = np.repeat(['m1', 'm2', 'h1', 'h2'], sizes)
    human = np.isin(specimen, ['h1', 'h2']).astype(float)
    position = rng.uniform(0, 1, len(specimen))
    library = rng.lognormal(np.log(2000), .4, len(specimen))
    return rng, specimen, human, position, library


def _nb_counts(rng, log_mu, alpha):
    """Gamma-Poisson draws (structures x genes) with per-gene dispersion ``alpha``."""
    return rng.poisson(rng.gamma(1 / alpha, np.exp(log_mu) * alpha)).astype(float)


@pytest.mark.filterwarnings('ignore:The design matrix is rank-deficient')
def test_irls_at_fixed_dispersion_matches_statsmodels_with_a_redundant_column():
    sm = pytest.importorskip('statsmodels.api')
    rng, specimen, human, position, library = _cohort(1, (70, 50, 60, 40))
    designs, w = contrast_designs(position, human, specimen, gam_internal_knots(position, 4))
    x = np.column_stack([designs['full'], designs['full'][:, 1] + designs['full'][:, 2]])
    off, alpha = np.log(library), np.array([.05, .3, 1.])
    y = _nb_counts(rng, off[:, None] + np.log(3e-3) + np.sin(3 * position)[:, None] + .7 * human[:, None], alpha)
    fit = nb_fit(sparse.csr_matrix(y), x, w, off, alpha)
    assert fit['converged'].all() and fit['beta'].shape == (3, x.shape[1])
    for g in range(3):
        ref = sm.GLM(y[:, g], x, family=sm.families.NegativeBinomial(alpha=alpha[g]), offset=off,
                     var_weights=w).fit(tol=1e-12, maxiter=300)
        assert fit['deviance'][g] == pytest.approx(ref.deviance, rel=1e-6)
        # R glm's criterion (relative deviance change < 1e-8) pins the deviance far more tightly than the
        # means: Fisher scoring converges linearly here, and at a = 1 the means are good to about 1e-4.
        np.testing.assert_allclose(np.exp(x @ fit['beta'][g] + off), ref.mu, rtol=5e-4)


def test_dispersion_is_the_weighted_ml_and_flags_a_gene_without_overdispersion():
    from scipy.optimize import minimize
    from scipy.stats import nbinom
    rng, specimen, human, position, library = _cohort(2, (80, 50, 70, 40))
    designs, w = contrast_designs(position, human, specimen, gam_internal_knots(position, 4))
    x, off = designs['level'], np.log(library)
    log_mu = off[:, None] + np.log(4e-3) + .8 * position[:, None] - .5 * human[:, None]
    y = _nb_counts(rng, np.repeat(log_mu, 2, axis=1), np.array([.1, .5]))
    # Binomial counts proportional to library: variance below the mean, so the ML sits on the Poisson bound.
    y = np.column_stack([y, rng.binomial(np.round(library / 250).astype(int), .5)])
    disp = nb_dispersion(y, x, w, off)
    assert disp['gene_model'] == GENE_MODEL and disp['converged'].all()
    assert list(disp['at_lower']) == [False, False, True] and disp['alpha'][2] == pytest.approx(A_MIN)
    beta = nb_fit(y, x, w, off, disp)['beta']

    def nll(par, col):
        mu, theta = np.exp(x @ par[:-1] + off), np.exp(-par[-1])
        return -(w * nbinom.logpmf(y[:, col], theta, theta / (theta + mu))).sum()

    for col in range(2):
        start = np.r_[np.linalg.lstsq(x, np.log(y[:, col] + .5) - off, rcond=None)[0], 0.]
        ref = minimize(nll, start, args=(col,), method='L-BFGS-B',
                       options={'maxiter': 5000, 'ftol': 1e-15, 'gtol': 1e-10})
        assert disp['alpha'][col] == pytest.approx(np.exp(ref.x[-1]), rel=1e-3)
        assert nll(np.r_[beta[col], np.log(disp['alpha'][col])], col) <= ref.fun + 1e-6


def test_nested_lr_is_the_statsmodels_deviance_difference_with_partial_f_df():
    sm = pytest.importorskip('statsmodels.api')
    from scipy.linalg import orth
    rng, specimen, human, position, library = _cohort(3)
    labels = np.where(position < .4, 'S1', np.where(position < .7, 'S2', 'S3'))
    designs, w = contrast_designs(position, human, specimen, gam_internal_knots(position, 4), segments=labels)
    off, alpha = np.log(library), np.array([.1, .4])
    y = _nb_counts(rng, off[:, None] + np.log(5e-3) + np.cos(2 * position)[:, None] * [1, -1]
                   + (human * position)[:, None], alpha)
    # 'union' carries redundant columns.
    comparisons = {'T_spatial': ('level', 'full'), 'T_position_given_segments': ('discrete', 'union')}
    got = nested_nb_lr(y, designs, w, off, alpha, comparisons)
    assert got['converged'].all()
    assert got['df'] == nested_partial_f(np.log1p(y), designs, w, comparisons)['df']
    for g in range(2):
        # statsmodels' IRLS stalls on the eight redundant 'union' columns; a deviance depends only on the column
        # space, so the reference fits use scipy's orthonormal basis of each design.
        dev = {k: sm.GLM(y[:, g], orth(designs[k]), family=sm.families.NegativeBinomial(alpha=alpha[g]),
                         offset=off, var_weights=w).fit(tol=1e-12, maxiter=300).deviance
               for k in ('level', 'full', 'discrete', 'union')}
        for name, (reduced, full) in comparisons.items():
            assert got[name][g] == pytest.approx(dev[reduced] - dev[full], rel=1e-6, abs=1e-6)
    with pytest.raises(ValueError, match='not nested'):
        nested_nb_lr(y, designs, w, off, alpha, {'bad': ('full', 'discrete')})


def test_block_size_and_workers_do_not_change_results():
    rng, specimen, human, position, library = _cohort(4, (50, 40, 50, 30))
    designs, w = contrast_designs(position, human, specimen, gam_internal_knots(position, 4))
    off = np.log(library)
    y = _nb_counts(rng, off[:, None] + np.log(4e-3) + np.sin(3 * position)[:, None] * rng.normal(size=7)
                   + human[:, None] * rng.normal(size=7), rng.uniform(.05, .5, 7))
    settings = ((3, 1), (50, 1), (3, 2))
    disp = [nb_dispersion(y, designs['full'], w, off, block_size=b, n_jobs=j) for b, j in settings]
    lr = [nested_nb_lr(y, designs, w, off, d, block_size=b, n_jobs=j, quasi=True) for d, (b, j) in zip(disp, settings)]
    for (one, wide, parallel), keys in ((disp, ('alpha', 'deviance')), (lr, ('T_level', 'T_spatial_ql', 'T_total'))):
        np.testing.assert_array_equal(one['converged'], wide['converged'])
        np.testing.assert_array_equal(one['converged'], parallel['converged'])
        for key in keys:
            np.testing.assert_array_equal(one[key], parallel[key])
            # A wider block runs other BLAS kernels, so only rounding may differ.
            np.testing.assert_allclose(one[key], wide[key], rtol=1e-10)


def test_lr_stays_calibrated_where_the_log_cp10k_partial_f_inflates():
    """Pure species offset, 20-fold, at low mouse counts: no position-dependent species difference."""
    rng, specimen, human, position, library = _cohort(0, (150, 150, 150, 150))
    designs, w = contrast_designs(position, human, specimen, gam_internal_knots(position, 6))
    off = np.log(library)
    log_mu = off + np.log(.2 / 2000) + 1.5 * np.sin(2 * np.pi * position) + 3 * human
    y = _nb_counts(rng, np.repeat(log_mu[:, None], 60, axis=1), np.full(60, .3))
    disp = nb_dispersion(y, designs['full'], w, off)
    lr = nested_nb_lr(y, designs, w, off, disp, {'T_spatial': ('level', 'full')})
    assert lr['converged'].all() and lr['df']['T_spatial'][0] == 6
    assert .7 < lr['T_spatial'].mean() / 6 < 1.4
    gauss = nested_partial_f(np.log1p(y * 1e4 / library[:, None]), designs, w, {'T_spatial': ('level', 'full')})
    assert gauss['T_spatial'].mean() > 5


def test_difference_curve_and_sandwich_errors():
    rng, specimen, human, position, library = _cohort(6, (100, 100, 100, 100))
    knots, grid = gam_internal_knots(position, 4), np.linspace(.1, .9, 9)
    designs, w = contrast_designs(position, human, specimen, knots)
    off, x = np.log(library), designs['full']
    truth = .4 + .6 * (position - .5)
    y = _nb_counts(rng, (off + np.log(1e-2) + np.sin(3 * position) + human * truth)[:, None] + rng.normal(0, .2, 5),
                   np.full(5, .1))
    disp = nb_dispersion(y, x, w, off)
    curves = nb_group_difference_curves(y, library, position, human, specimen, knots, grid, disp)
    assert curves['converged'].all()
    np.testing.assert_allclose(curves['delta'].mean(axis=0), .4 + .6 * (grid - .5), atol=.1)
    assert np.all(np.abs(curves['delta'] - (.4 + .6 * (grid - .5))) < 4 * curves['se'])
    # Independent per-gene loop on the original columns [intercept, basis, group, group x basis, specimens].
    base = make_gam_design(grid, knots)
    p_b = base.shape[1]
    rows = np.zeros((len(grid), x.shape[1]))
    rows[:, p_b], rows[:, p_b + 1:2 * p_b] = 1, base[:, 1:]
    fit = nb_fit(y, x, w, off, disp, contrasts={'delta': rows})
    # Different IRLS starts (warm from the level model here) agree to the deviance criterion's precision.
    np.testing.assert_allclose(fit['contrasts']['delta']['value'], curves['delta'], rtol=1e-4)
    np.testing.assert_allclose(fit['contrasts']['delta']['se_hc3'], curves['se'], rtol=1e-4)
    for g in range(5):
        mu, a = np.exp(x @ fit['beta'][g] + off), disp['alpha'][g]
        info = w * mu / (1 + a * mu)
        inv = np.linalg.inv(x.T @ (info[:, None] * x))
        lev = info * np.einsum('ij,jk,ik->i', x, inv, x)
        score = w * (y[:, g] - mu) / (1 + a * mu) / (1 - lev)
        for key, middle in (('se_hc3', score ** 2), ('se_model', w * info)):
            cov = inv @ (x.T @ (middle[:, None] * x)) @ inv
            expected = np.sqrt(np.einsum('ij,jk,ik->i', rows, cov, rows))
            np.testing.assert_allclose(fit['contrasts']['delta'][key][g], expected, rtol=1e-8)


def test_hc0_form_equals_statsmodels_sandwich_with_expected_information():
    sm = pytest.importorskip('statsmodels.api')
    from statsmodels.stats.sandwich_covariance import cov_white_simple
    rng, specimen, human, position, library = _cohort(5)
    designs, w = contrast_designs(position, human, specimen, gam_internal_knots(position, 4))
    assert np.allclose(w, 1)
    x, off, alpha = designs['full'], np.log(library), np.array([.2, .6])
    y = _nb_counts(rng, (off + np.log(5e-3) + np.sin(3 * position) + .5 * human)[:, None], alpha)
    beta = nb_fit(y, x, w, off, alpha)['beta']
    q = _orthobasis(x)
    mu = np.exp(beta @ x.T + off).T
    hc0 = _contrast_se(y, mu, q, _outer(q), w, alpha, np.linalg.pinv(x) @ q, hc3=False)[0]
    for g in range(2):
        model = sm.GLM(y[:, g], x, family=sm.families.NegativeBinomial(alpha=alpha[g]), offset=off)
        # statsmodels' IRLS default uses the observed Hessian; the IRLS bread here is the expected information.
        cov = cov_white_simple((model.score_obs(beta[g]), np.linalg.inv(model.hessian(beta[g], observed=False))),
                               use_correction=False)
        np.testing.assert_allclose(hc0[g], np.sqrt(np.diag(cov)), rtol=1e-8)


def test_a_gene_without_counts_is_filled_but_a_separated_gene_keeps_its_limiting_lr():
    rng, specimen, human, position, library = _cohort(7)
    designs, w = contrast_designs(position, human, specimen, gam_internal_knots(position, 4))
    off = np.log(library)
    y = _nb_counts(rng, (off + np.log(4e-3) + np.sin(3 * position))[:, None] + rng.normal(0, .3, 5), np.full(5, .2))
    y[:, 0] = 0                      # never detected: no information
    y[human == 1, 1] = 0             # absent from one species: boundary MLE, the species level runs to -inf
    disp = nb_dispersion(y, designs['full'], w, off)
    assert list(disp['converged']) == [False, True, True, True, True]
    assert list(disp['separated']) == [True, True, False, False, False]
    for dispersion in (disp, np.full(5, .2)):   # a plain array: the gene without counts still fails
        got = nested_nb_lr(y, designs, w, off, dispersion)
        assert list(got['converged']) == [False, True, True, True, True]
        assert list(got['separated']) == [True, True, False, False, False]
        for name in ('T_level', 'T_spatial', 'T_total'):
            assert got[name][0] == np.median(got[name][1:])
        # No positional information in the empty species, but a real offset.
        assert got['T_spatial'][1] < 1e-3 and got['T_level'][1] > 10 * got['T_level'][2:].max()
    raw = nested_nb_lr(y, designs, w, off, disp, fill=None)
    assert np.isnan(raw['T_level'][0]) and np.isfinite(raw['T_level'][1:]).all()
    assert raw['T_level'][1] == nested_nb_lr(y, designs, w, off, disp)['T_level'][1]


def test_species_trajectories_return_the_gaussian_keys_and_more():
    rng, specimen, human, position, library = _cohort(8)
    grid = np.linspace(.1, .9, 7)
    y = _nb_counts(rng, (np.log(library) + np.log(5e-3) + np.sin(3 * position) + .6 * human)[:, None]
                   + rng.normal(0, .3, 4), np.full(4, .2))
    gauss = fit_nested_trajectories(np.log1p(y * 1e4 / library[:, None]), position, human, specimen, grid,
                                    return_residuals=True)
    designs, w = contrast_designs(position, human, specimen, gauss['knots'])
    disp = nb_dispersion(y, designs['full'], w, np.log(library))
    nb = nb_species_trajectories(sparse.csr_matrix(y), library, position, human, specimen, grid, disp,
                                 return_residuals=True)
    assert set(gauss) < set(nb) and {'se_model', 'converged', 'separated', 'gene_model'} <= set(nb)
    for key in gauss:
        assert np.shape(nb[key]) == np.shape(gauss[key])
    np.testing.assert_array_equal(nb['design_df'], gauss['design_df'])
    np.testing.assert_allclose(nb['delta'], nb['human'] - nb['mouse'], atol=1e-10)
    assert nb['converged'].all() and not nb['separated'].any() and np.all(nb['beta_human_level'] > 0)
