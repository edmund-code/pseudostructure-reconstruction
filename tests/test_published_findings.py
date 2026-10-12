import numpy as np
import pandas as pd
import pytest

from pseudospace import published_findings as pf


def test_sage_counts_use_the_smallest_value_as_one_tag():
    values = pd.DataFrame({'S1': [0., .2, 1., 2.6], 'S3': [.1, .1, 0., 5.]})
    counts, library = pf.sage_counts(values)
    assert counts.S1.tolist() == [0, 1, 5, 13] and counts.S3.tolist() == [1, 1, 0, 50]
    assert library.S1 == pytest.approx(5e4) and library.S3 == pytest.approx(1e5)
    with pytest.raises(ValueError):
        pf.sage_counts(pd.DataFrame({'x': [0., 0.]}))


def test_sage_contrast_and_zonated():
    d, se = pf.sage_contrast([40, 10], [10, 10], 1e5, 1e5)
    assert d[0] == pytest.approx(np.log2(40.5 / 10.5)) and d[1] == pytest.approx(0)
    assert se[0] == pytest.approx(np.sqrt(1 / 40.5 + 1 / 10.5) / np.log(2))
    assert pf.zonated(d, se).tolist() == [True, False]


def test_strength_agreement_and_trend():
    rng = np.random.default_rng(1)
    a = rng.uniform(-3, 3, 4000)
    # agreement rises with |a|: weak genes are noise, strong genes share the sign
    b = np.where(np.abs(a) >= 2, a, np.where(np.abs(a) >= 1, a * rng.choice([1, 1, -1], 4000), rng.normal(size=4000)))
    frame, (z, p) = pf.strength_agreement(a, b, np.abs(a) >= .5)
    assert frame.share.iloc[-1] == 1 and frame.share.iloc[0] < .65
    assert z > 3 and p < 1e-3
    z0, p0 = pf.cochran_armitage([50, 50, 50], [100, 100, 100])
    assert z0 == pytest.approx(0) and p0 == pytest.approx(.5)


def test_slope_through_origin():
    x = np.linspace(-2, 2, 50)
    est, lo, hi = pf.slope_through_origin(x, .3 * x, n_boot=200)
    assert est == pytest.approx(.3) and lo == pytest.approx(.3) and hi == pytest.approx(.3)


def test_permutation_spearman_partial_removes_a_shared_driver():
    rng = np.random.default_rng(2)
    m = rng.normal(size=600)
    x = -m + .3 * rng.normal(size=600)          # both interactions follow the mouse contrast (flattening)
    y = -m + .3 * rng.normal(size=600)
    rho, p = pf.permutation_spearman(x, y, n_perm=500, seed=0)
    assert rho > .8 and p < .01
    rho_p, p_p = pf.permutation_spearman(x, y, covariates=m[:, None], n_perm=500, seed=0)
    assert abs(rho_p) < .15 and p_p > .01
    shared = rng.normal(size=600)                # a real interaction beyond the driver
    rho_s, p_s = pf.permutation_spearman(x + shared, y + shared, covariates=m[:, None], n_perm=500, seed=0)
    assert rho_s > .5 and p_s < .01


def test_library_class_rules():
    assert pf.library_class([.5] * 9 + [-.1])['class'] == 'zonated'
    assert pf.library_class([.5] * 8 + [-.1] * 2)['class'] == 'attenuated'      # sign test p 0.055
    assert pf.library_class([-.6] * 9 + [.1])['class'] == 'reversed'
    assert pf.library_class([.1, -.1, .2, -.2, 0, .05])['class'] == 'absent'
    assert pf.library_class([.4, .5, -.3, -.2, .6, .3])['class'] == 'attenuated'
    assert pf.library_class([.4, .5, np.nan])['class'] == 'not testable'


def test_donor_class_rules():
    assert pf.donor_class([.5, .6, .4, .7, .5, .45])['class'] == 'zonated'
    assert pf.donor_class([-.5, -.6, -.4, -.7, -.5, -.45])['class'] == 'reversed'
    assert pf.donor_class([.1, -.1, .05, 0, .2, -.15])['class'] == 'absent'
    assert pf.donor_class([.5, .1, .4, .2, .6, .3])['class'] == 'attenuated'      # 4 of 6 donors ≥ 0.25
    assert pf.donor_class([.5, .6])['class'] == 'not testable'


def test_segment_contrasts_and_donor_rates():
    rates = pd.DataFrame({'S1': [1., 0.], 'S2': [0., 0.], 'S3': [-1., 3.]})
    c = pf.segment_contrasts(rates)
    assert c.loc[0].tolist() == [1.5, 0., -1.5] and c.loc[1, 'S3'] == 3
    meta = pd.DataFrame({'donor': ['a'] * 3 + ['b'] * 3, 'segment': ['S1', 'S2', 'S3'] * 2,
                         'cells': [60, 60, 60, 60, 10, 60]})
    sums = np.array([[10, 0], [20, 0], [30, 0], [1, 1], [1, 1], [1, 1]], float)
    out = pf.donor_segment_rates(sums, [100, 100, 100, 10, 10, 10], meta, ['g1', 'g2'])
    assert list(out) == ['a']                    # donor b has too few S2 cells
    assert out['a'].loc['g1', 'S3'] == pytest.approx(np.log2(.3e6 + 1))


def test_gene_class_order():
    pairs = lambda q: np.tile(q, (4, 1))        # noqa: E731
    assert pf.gene_class([1, 1, 1], pairs([1, 1, 1]), 0) == 'present in all segments'
    assert pf.gene_class([1, .1, .5], pairs([1, .1, .5]), 0) == 'restricted to the listed segment'
    assert pf.gene_class([.1, .8, .1], pairs([.1, .8, .1]), 0) == 'elsewhere only'
    assert pf.gene_class([.1, -.5, .1], pairs([.1, -.5, .1]), 0) == 'opposite somewhere'
    assert pf.gene_class([.3, .1, .1], np.array([[.3, .1, .1]] * 3 + [[-.1, .1, .1]]), 0) == 'not reproduced'
    assert pf.gene_class([np.nan, 1, 1], pairs([1, 1, 1]), 0) == 'not measured'


def test_mantel_haenszel_matches_a_hand_count_and_stratifies():
    exposed = np.array([1] * 20 + [0] * 20, bool)
    outcome = np.array([1] * 10 + [0] * 10 + [1] * 2 + [0] * 18, bool)
    res = pf.mantel_haenszel(exposed, outcome, np.zeros(40))
    assert res['or_mh'] == pytest.approx(10 * 18 / (10 * 2)) and res['p'] < .01 and res['a'] == 10
    # confounded by strata: within each stratum no association
    strata = np.repeat([0, 1], 20)
    e = np.r_[np.ones(16), np.zeros(4), np.ones(4), np.zeros(16)].astype(bool)
    o = np.r_[np.ones(8), np.zeros(8), np.ones(2), np.zeros(2), np.zeros(4), np.zeros(16)].astype(bool)
    assert pf.mantel_haenszel(e, o, strata)['or_mh'] == pytest.approx(1)


def test_program_curve_centres_and_takes_the_median():
    curves = np.log(np.array([[1, 2, 4], [10, 20, 40], [5, 5, 5], [1, 4, 16], [2, 4, 8]], float))
    prof, n = pf.program_curve(curves, [True] * 5)
    assert n == 5 and prof[2] - prof[0] == pytest.approx(2)
    assert np.isnan(pf.program_curve(curves, [True, True, False, False, False])[0]).all()


def test_reference_curve_removes_the_shape_every_gene_shares():
    rng = np.random.default_rng(3)
    grid = np.linspace(0, 1, 50)
    shared = np.sin(3 * grid)                                   # e.g. a composition trend along the axis
    genes = shared + rng.normal(0, .05, (200, 50)) + rng.normal(0, 1, (200, 1))
    genes[:20] += 2 * (grid > .6)                               # a real program switching at 0.6
    ref = pf.reference_curve(genes, np.ones(200, bool))
    random_prog, _ = pf.program_curve(genes[100:120], np.ones(20, bool), reference=ref, log=1.)
    real_prog, _ = pf.program_curve(genes[:20], np.ones(20, bool), reference=ref, log=1.)
    assert np.ptp(random_prog) < .2
    assert real_prog[grid > .7].mean() - real_prog[grid < .5].mean() > 1.8


def test_shape_descriptors_classes():
    u = np.linspace(0, 3, 301)
    sharp = 1 / (1 + np.exp(-(u - 1.6) / .05))
    d = pf.shape_descriptors(sharp, u)
    assert d['shape'] == 'rising' and d['switch'] == 'sharp switch' and d['position'] == pytest.approx(1.6, abs=.02)
    ramp = u / 3
    g = pf.shape_descriptors(-ramp, u)
    assert g['shape'] == 'falling' and g['switch'] == 'gradient' and g['position'] == pytest.approx(1.5, abs=.02)
    peak = np.exp(-((u - 1.5) / .2) ** 2)
    p = pf.shape_descriptors(peak, u)
    assert p['shape'] == 'interior peak' and p['position'] == pytest.approx(1.5, abs=.02)
    assert pf.shape_descriptors(.1 * ramp, u)['shape'] == 'flat'
    assert pf.shape_descriptors(np.sin(4 * u), u)['shape'] in ('mixed', 'interior peak', 'interior trough')


def test_species_statement_and_cross_species():
    a = {'shape': 'rising', 'position': 1.5, 'width': .3, 'switch': 'sharp switch'}
    b = {'shape': 'rising', 'position': 1.7, 'width': .4, 'switch': 'sharp switch'}
    s = pf.species_statement(a, b)
    assert s['agree'] and s['position'] == pytest.approx(1.6) and s['flag'] == 'human S2: weak cell'
    assert pf.species_statement(a, {**b, 'position': 2.0})['shape'] == 'rising, positions disagree'
    assert pf.species_statement(a, {**b, 'shape': 'flat'})['shape'] == 'specimens disagree'
    m = pf.species_statement(a, b, weak=None)
    h = pf.species_statement({**a, 'position': 2.1}, {**b, 'position': 2.2}, weak=None)
    assert pf.cross_species(m, h) == 'later in human by 0.55 segment units'
    assert pf.cross_species(m, m) == 'same place within registration error'
    assert pf.segment_position(1.62) == 'S2 0.62'


def test_pathway_label_takes_the_largest_overlap_then_library_order():
    libs = {'Hallmark': {'A': ['g1', 'g2', 'g3'], 'B': ['g1', 'g2', 'g3', 'g4']},
            'KEGG': {'C': ['g1', 'g2', 'g3', 'g4']}}
    assert pf.pathway_label(['g1', 'g2', 'g3', 'g4', 'g9'], libs) == ('Hallmark', 'B', ['g1', 'g2', 'g3', 'g4'])
    assert pf.pathway_label(['g1', 'g2'], libs) is None
