import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.stats import false_discovery_control, hypergeom

from pseudospace.discrete_comparators import (family_ora, membership, segment_f_test, segment_step_designs,
                                              segment_step_lr, swap_rule)

SPECIMENS = ['m1', 'm2', 'h1', 'h2']
SEGMENTS = ['S1', 'S2', 'S3']
STEP_GENE = 0                       # a gene whose human - mouse difference changes between segments
STEPS = np.array([.6, 0., -.6])     # its true per-segment log ratio, human - mouse


def _cohort(seed=11, n_per=300, n_genes=120, alpha=.05, specimen_steps=0.):
    """NB counts: segment shape and species offsets everywhere; STEP_GENE adds human steps; optional specimen noise."""
    rng = np.random.default_rng(seed)
    specimen = np.repeat(SPECIMENS, n_per)
    segment = np.tile(np.repeat(SEGMENTS, n_per // 3), len(SPECIMENS))
    human = np.isin(specimen, ['h1', 'h2']).astype(float)
    library = rng.lognormal(np.log(4000), .3, len(specimen))
    seg_index = pd.Series(range(3), index=SEGMENTS)[segment].to_numpy()
    shape = rng.normal(0, .4, (3, n_genes))[seg_index]                     # shared segment shape per gene
    offset = rng.normal(0, .3, n_genes)[None, :] * human[:, None]            # species offset per gene
    noise = rng.normal(0, specimen_steps, (len(SPECIMENS), 3, n_genes))      # specimen x segment deviations
    spec_index = pd.Series(range(4), index=SPECIMENS)[specimen].to_numpy()
    rate = np.log(3e-3) + shape + offset + noise[spec_index, seg_index]
    rate[:, STEP_GENE] += human * STEPS[seg_index]
    mu = library[:, None] * np.exp(rate)
    counts = rng.negative_binomial(1 / alpha, 1 / (1 + alpha * mu))
    return {'counts': sparse.csr_matrix(counts.astype(float)), 'library': library, 'specimen': specimen,
            'segment': segment, 'human': human, 'alpha': np.full(n_genes, alpha)}


@pytest.fixture(scope='module')
def cohort():
    return _cohort()


def _rank(x):
    return np.linalg.matrix_rank(x)


def test_designs_are_nested_with_the_expected_degrees_of_freedom(cohort):
    c = cohort
    designs, weights, labels = segment_step_designs(c['segment'], c['human'], c['specimen'])
    assert list(labels) == SEGMENTS
    ranks = {k: _rank(v) for k, v in designs.items()}
    assert ranks['steps'] - ranks['level'] == 2 and ranks['spec'] - ranks['steps'] == 4   # condition split
    for small, big in (('level', 'steps'), ('steps', 'spec')):
        assert _rank(np.column_stack([designs[small], designs[big]])) == ranks[big]
    # each group half the weight, each specimen an equal share
    for name in SPECIMENS:
        assert np.isclose(weights[c['specimen'] == name].sum(), len(weights) / 4)
    relabel = np.isin(c['specimen'], ['m1', 'h1']).astype(float)
    designs_r, _, _ = segment_step_designs(c['segment'], relabel, c['specimen'], nuisance=c['human'])
    ranks_r = {k: _rank(v) for k, v in designs_r.items()}
    assert ranks_r['steps'] - ranks_r['level'] == 2 and ranks_r['spec'] - ranks_r['steps'] == 2   # relabeling


def test_group_must_be_constant_within_specimen(cohort):
    c = cohort
    bad = c['human'].copy()
    bad[0] = 1 - bad[0]
    with pytest.raises(ValueError):
        segment_step_designs(c['segment'], bad, c['specimen'])


def test_step_lr_finds_the_step_gene_and_recovers_its_segment_ratios(cohort):
    c = cohort
    out = segment_step_lr(c['counts'], c['library'], c['segment'], c['human'], c['specimen'], c['alpha'])
    assert (out['df1'], out['df2']) == (2, 4)
    assert out['segments'] == SEGMENTS
    null = np.delete(out['T_step'], STEP_GENE)
    assert out['converged'].all() and not out['separated'].any()
    assert out['T_step'][STEP_GENE] > 50 and np.median(null) < 3        # chi2(2) median is 1.39
    # the per-segment difference: STEPS plus the gene's species offset, which is common to all segments
    delta = out['delta'][STEP_GENE]
    np.testing.assert_allclose(delta - delta.mean(), STEPS - STEPS.mean(), atol=.12)
    assert (out['se_hc3'][STEP_GENE] > 0).all() and (out['se_hc3'][STEP_GENE] < .2).all()
    # a null gene's steps are equal across segments
    spread = np.ptp(np.delete(out['delta'], STEP_GENE, axis=0), axis=1)
    assert np.median(spread) < .15


def test_segment_f_test_is_calibrated_under_specimen_step_noise():
    pytest.importorskip('rpy2')
    from rpy2.robjects.packages import isinstalled
    if not isinstalled('limma'):
        pytest.skip('R package limma is not installed')
    c = _cohort(seed=5, n_genes=200, specimen_steps=.15)
    frame = segment_f_test(c['counts'], c['library'], c['segment'], c['human'], c['specimen'], c['alpha'],
                           genes=[f'g{i}' for i in range(200)],
                           covariate=np.log(np.asarray(c['counts'].mean(axis=0)).ravel()))
    assert {'F', 'p', 'z', 'q', 'T_step', 'LR_spec', 'df2', 'delta_S1', 'se_S3'} <= set(frame.columns)
    assert (frame.df2 == 4).all()
    null = frame.drop(index='g0')
    # specimen step noise inflates T_step, but the specimen-step denominator absorbs it
    assert np.mean(null.p <= .05) < .12
    relabel = np.isin(c['specimen'], ['m1', 'h1']).astype(float)
    swap = segment_f_test(c['counts'], c['library'], c['segment'], relabel, c['specimen'], c['alpha'],
                          nuisance=c['human'], covariate=np.log(np.asarray(c['counts'].mean(axis=0)).ravel()))
    assert (swap.df2 == 2).all() and np.mean(swap.p <= .05) < .12


def test_family_ora_matches_hypergeometric_and_family_bh():
    universe = [f'g{i}' for i in range(200)]
    sets = {'a': universe[:20], 'b': universe[10:40], 'c': universe[100:130] + ['not_in_universe']}
    member = membership(sets, universe)
    assert member.shape == (3, 200) and member.sum() == 80
    up = np.zeros(200, bool)
    up[:15] = True
    down = np.zeros(200, bool)
    down[150:160] = True
    out = family_ora([up, down, np.zeros(200, bool)], member)
    expected_up = hypergeom.sf(np.array([15, 5, 0]) - 1, 200, np.array([20, 30, 30]), 15)
    np.testing.assert_allclose(out.p_0, expected_up)
    np.testing.assert_allclose(out.p_2, 1.)                              # an empty list tests nothing
    p = out[['p_0', 'p_1', 'p_2']].to_numpy()
    np.testing.assert_allclose(out[['q_0', 'q_1', 'q_2']].to_numpy(), false_discovery_control(p.ravel()).reshape(p.shape))
    assert out.call.tolist() == [True, False, False]
    np.testing.assert_allclose(out.calibration_p, np.minimum(1, 3 * p.min(axis=1)))
    # lists from different backgrounds: one member matrix per list, one BH over the family
    small = universe[:100]
    mixed = family_ora([up, up[:100]], [member, membership(sets, small)])
    np.testing.assert_allclose(mixed.p_1, hypergeom.sf(np.array([15, 5, 0]) - 1, 100, np.array([20, 30, 0]), 15))
    with pytest.raises(ValueError):
        family_ora([up, down], [member])


def test_swap_rule_reproduces_notebook_71_verdicts():
    rng = np.random.default_rng(0)

    def pvalues(share, n=3026):
        p = rng.uniform(.06, 1, n)
        p[:int(round(share * n))] = .01
        return p

    hvm = swap_rule(83, [1, 0], 1513, pvalues(.071))                     # notebook 71, human vs mouse: passes
    assert hvm['pass'] == 'pass' and np.isclose(hvm['NCDR'], .5 / 83)
    aki = swap_rule(103, [13, 29], 1469, pvalues(.120, 2938))             # AKI: relabeling passes, calibration fails
    assert (aki['relabeling rule'], aki['calibration'], aki['pass']) == ('pass', 'fail', 'fail')
    assert swap_rule(0, [0, 0], 1513, pvalues(0.))['pass'] == 'fail'      # no calls never passes
    assert swap_rule(100, [10, 80], 1513, pvalues(0.))['relabeling rule'] == 'fail'   # one relabeling >= 5%
    two = swap_rule([97, 152], [0, 0], 1513, pvalues(.002))               # several condition runs: their mean
    assert two['calls'] == 124.5 and two['pass'] == 'pass'
