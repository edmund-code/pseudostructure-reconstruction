import numpy as np
import pytest
from scipy import sparse

from pseudospace import gene_model_audit as ga


def _four_specimens(rng, n=400, genes=12):
    specimen = np.repeat(['h1', 'h2', 'm1', 'm2'], n)
    human = np.isin(specimen, ['h1', 'h2'])
    position = rng.uniform(0, 1, 4 * n)
    library = rng.uniform(500, 1500, 4 * n)
    base = np.exp(rng.uniform(np.log(3e-3), np.log(3e-2), genes))
    shape = np.exp(np.outer(np.sin(3 * position), rng.normal(0, .8, genes)) + np.outer(human * position, rng.normal(0, .6, genes)))
    mu = library[:, None] * base[None, :] * shape
    counts = rng.negative_binomial(5, 5 / (5 + mu))
    return counts, library, position, human.astype(float), specimen


def test_positional_lrs_reproduce_the_standard_model():
    pytest.importorskip('patsy')
    from pseudospace.pathway_pipelines import specimen_shape_lr, standard_gene_model
    from pseudospace.stats_gam import gam_internal_knots
    rng = np.random.default_rng(0)
    counts, library, position, human, specimen = _four_specimens(rng)
    model = standard_gene_model(counts, library, position, human, specimen)
    shape = specimen_shape_lr(counts, library, position, human, specimen, model.dispersion)
    out = ga.positional_lrs(counts, library, position, human, specimen, gam_internal_knots(position, 6))
    np.testing.assert_allclose(out['T_spatial'], model.statistics.T_spatial, rtol=1e-6)
    np.testing.assert_allclose(out['T_total'], model.statistics.T_total, rtol=1e-6)
    np.testing.assert_allclose(out['LR_spec'], shape['LR_spec'], rtol=1e-6)
    assert out['df1'] == 6 and out['df2'] == shape['df2'] == 12


def test_zero_check_and_residuals_accept_nb_and_flag_zero_inflation():
    rng = np.random.default_rng(1)
    mu = np.tile(rng.uniform(.2, 5, (3000, 1)), (1, 2))
    y = rng.negative_binomial(2, 2 / (2 + mu))
    y[:, 1] = np.where(rng.uniform(size=3000) < .2, 0, y[:, 1])          # 20% extra zeros in gene 2
    z = ga.zero_check(y, mu, np.array([.5, .5]))
    assert z.p[0] > .01 and z.p[1] < 1e-6 and z.observed[1] - z.nb[1] > .1
    r = ga.quantile_residuals(y[:, :1], mu[:, :1], np.array([.5]), rng)
    assert ga.residual_ks(r).D[0] < .03
    disp = ga.pearson_dispersion(y[:, :1], mu[:, :1], np.array([.5]), 1)
    assert disp.nb[0] == pytest.approx(1, abs=.1) and disp.poisson[0] > 1.5


def test_grouped_alpha_recovers_different_dispersions():
    rng = np.random.default_rng(2)
    mu = np.full((6000, 1), 4.)
    groups = np.repeat(['a', 'b'], 3000)
    size = np.where(groups == 'a', 10., 1.)[:, None]
    y = rng.negative_binomial(size, size / (size + mu))
    per, test = ga.grouped_alpha(y, mu, np.array([.3]), groups)
    assert per.loc[0, 'a'] == pytest.approx(.1, abs=.04) and per.loc[0, 'b'] == pytest.approx(1., abs=.15)
    assert test.p[0] < 1e-10


def test_poscounts_trend_and_shrinkage():
    rng = np.random.default_rng(3)
    sf = np.exp(rng.normal(0, .4, 300))
    rates = np.exp(rng.uniform(0, 3, 200))
    counts = rng.poisson(np.outer(sf, rates))
    est = ga.poscounts_size_factors(sparse.csr_matrix(counts))
    assert np.corrcoef(np.log(est), np.log(sf))[0, 1] > .97
    mean = np.geomspace(.1, 100, 500)
    alpha = (.05 + .5 / mean) * np.exp(rng.normal(0, .1, 500))
    a, b, trend = ga.dispersion_trend(alpha, mean)
    assert a == pytest.approx(.05, rel=.3) and b == pytest.approx(.5, rel=.3)
    shrunk, weight = ga.shrunken_alpha(alpha, trend, 10000, 7)
    assert weight > .99 and np.allclose(shrunk, alpha, rtol=.02)


def test_heldout_deviance_prefers_the_curve_over_a_flat_line():
    pytest.importorskip('patsy')
    from pseudospace.stats_gam import make_gam_design
    rng = np.random.default_rng(4)
    position = rng.uniform(0, 1, 2000)
    specimen = np.repeat(['a', 'b'], 1000)
    mu = 3 * np.exp(1.5 * np.sin(4 * position))
    y = rng.negative_binomial(5, 5 / (5 + mu))[:, None]
    folds = ga.cv_folds(specimen, 5, rng)
    assert np.bincount(folds).tolist() == [400] * 5
    flat, _ = ga.heldout_deviance(y, np.ones((2000, 1)), np.ones(2000), np.zeros(2000), folds)
    curve, ok = ga.heldout_deviance(y, make_gam_design(position, [.25, .5, .75]), np.ones(2000), np.zeros(2000), folds)
    assert ok.all() and curve[0] < flat[0]
