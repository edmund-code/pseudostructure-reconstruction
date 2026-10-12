"""program_localization on synthetic data: pseudobulks, closed-form contrasts, moderated t, set tests, the decision
rules of docs/results/pt-known-programs.md (section 5), gene spans and the program-file checks."""
import gzip

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from pseudospace.program_localization import (added_information, bin_rates, centred_contrast, gene_span, informative,
                                              localization_class, partition_effects, program_test, read_programs,
                                              replication_class, segment_contrast, segment_pseudobulk, segment_summary,
                                              stratified_rank_normal, transfer_class, file_sha256, moderated_t)


def test_segment_pseudobulk_sums_counts_and_floors_zeros():
    counts = sparse.csr_matrix(np.array([[2, 0], [3, 0], [5, 1], [0, 0]], float))
    library = np.array([10., 10., 20., 20.])
    out = segment_pseudobulk(counts, library, ['A', 'A', 'B', 'B'], ['S1', 'S3', 'S1', 'S1'],
                             {'S1': ['S1'], 'PT': None}, genes=['g1', 'g2'])
    assert out['S1'].loc['g1', 'A'] == pytest.approx(np.log2(2.5 / 10))
    assert out['S1'].loc['g2', 'A'] == pytest.approx(np.log2(.5 / 10))          # zero count floored at 0.5
    assert out['PT'].loc['g1', 'A'] == pytest.approx(np.log2(5.5 / 20))
    assert out['S1'].loc['g1', 'B'] == pytest.approx(np.log2(5.5 / 40))
    with pytest.raises(ValueError):
        segment_pseudobulk(counts, library, ['A', 'A', 'B', 'B'], ['S1', 'S3', 'S1', 'S1'], {'S3': ['S3']})


def test_partition_effects_equal_ols_on_four_samples():
    rng = np.random.default_rng(1)
    genes = [f'g{i}' for i in range(200)]
    rates = pd.DataFrame(rng.normal(size=(200, 4)), index=genes, columns=['c1', 'c2', 'r1', 'r2'])
    rates.iloc[:, :2] += .2                                                     # a shared case shift (centred away)
    designs = {'condition': np.c_[np.ones(4), [1, 1, 0, 0]],
               'relabeling 1': np.c_[np.ones(4), [1, 1, 0, 0], [1, 0, 1, 0]],
               'relabeling 2': np.c_[np.ones(4), [1, 1, 0, 0], [1, 0, 0, 1]]}
    for label, X in designs.items():
        out = partition_effects(rates, ('c1', 'c2'), ('r1', 'r2'), label)
        y = rates.to_numpy().T
        beta, rss, *_ = np.linalg.lstsq(X, y, rcond=None)
        df = 4 - X.shape[1]
        cov = np.linalg.inv(X.T @ X)[-1, -1]
        se2 = rss / df * cov
        effect = beta[-1] - np.median(beta[-1])                                 # centring by the median gene
        np.testing.assert_allclose(out.effect, effect, atol=1e-10)
        assert (out.df == df).all()
        # with every contrast's median gene at 0 (centring is then a no-op) s2 is exactly OLS's se^2
        one = pd.DataFrame(np.zeros((5, 4)), columns=rates.columns, index=list('abcde'))
        one.loc['a'] = rates.iloc[0]
        exact = partition_effects(one, ('c1', 'c2'), ('r1', 'r2'), label)
        b1, r1, *_ = np.linalg.lstsq(X, one.loc['a'].to_numpy(), rcond=None)
        assert exact.loc['a', 's2'] == pytest.approx(float(r1[0]) / df * cov)
        assert exact.loc['a', 'effect'] == pytest.approx(b1[-1])
    plain = pd.DataFrame({'c1': [1., 0, 0], 'c2': [3., 0, 0], 'r1': [0., 0, 0], 'r2': [2., 0, 0]}, index=['a', 'b', 'c'])
    cond = partition_effects(plain, ('c1', 'c2'), ('r1', 'r2'))
    assert cond.loc['a', 'effect'] == pytest.approx(1.)                        # (1 + 3)/2 - (0 + 2)/2, median 0
    assert cond.loc['a', 's2'] == pytest.approx(((-2.) ** 2 + (-2.) ** 2) / 4)


def test_segment_contrast_and_centring():
    genes = ['a', 'b', 'c']
    rates = {s: pd.DataFrame({'m1': v, 'm2': v}, index=genes) for s, v in
             {'S1': [0., 0, 0], 'S2': [0., 0, 0], 'S3': [3., 1, 1]}.items()}
    out = segment_contrast(rates, ['m1', 'm2'], 'S3')
    np.testing.assert_allclose(out.effect, [2., 0, 0])                         # median gene (b, c) centred to 0
    np.testing.assert_allclose(out.s2, 0)
    with pytest.raises(ValueError):
        segment_contrast(rates, ['m1'], 'S3')
    assert centred_contrast(rates['S3'], {'m1': 1.}).tolist() == [2., 0., 0.]


def test_moderated_t_shrinks_variances_towards_the_prior():
    pytest.importorskip('rpy2')
    rng = np.random.default_rng(2)
    sigma2 = 8 / rng.chisquare(8, 2000)                                         # inverse-gamma gene variances, prior df 8
    s2 = sigma2 * rng.chisquare(2, 2000) / 2
    effect = rng.normal(size=2000)
    out = moderated_t(effect, s2, 2.)
    assert np.isfinite(out.d0).all() and (out.d0 > 0).all()
    small, large = s2 < np.quantile(s2, .1), s2 > np.quantile(s2, .9)
    assert (out.s2_post[small] > s2[small]).all() and (out.s2_post[large] < s2[large]).all()
    assert ((out.p > 0) & (out.p <= 1)).all()


def test_program_test_direction_handling():
    pytest.importorskip('rpy2')
    rng = np.random.default_rng(3)
    genes = [f'g{i}' for i in range(500)]
    stat = pd.Series(rng.normal(size=500), index=genes)
    stat.iloc[:20] += 3                                                         # a planted case-high program
    programs = {'up': genes[:20], 'other': genes[100:130]}
    out = program_test(stat, programs, {'up': 1, 'other': -1})
    assert out.loc['up', 'p_claimed'] < 1e-6 and out.loc['up', 'p_opposite'] > .99
    claimed_down = program_test(stat, {'up': genes[:20]}, {'up': -1})
    assert claimed_down.loc['up', 'p_claimed'] > .99


def test_localization_classes():
    pairs = lambda q: np.tile(q, (4, 1))                                        # noqa: E731
    s3 = np.array([0., .1, 1.])
    out = localization_class(.01, s3, pairs(s3), [[0, 0, .1], [0, 0, -.1]])
    assert out['class'] == 'localized to S3' and out['present'] == ['S3'] and out['unresolved'] == ['S1', 'S2']
    flat = np.array([.6, .65, .7])
    assert localization_class(.4, flat, pairs(flat), np.zeros((2, 3)))['class'] == 'uniform within resolution'
    mixed = pairs(s3)
    mixed[0] = [1.2, .1, 1.]                                                    # one pair peaks in S1
    assert localization_class(.01, s3, mixed, np.zeros((2, 3)))['class'] == 'graded'
    assert localization_class(.6, np.zeros(3), pairs(np.zeros(3)) + .01, np.zeros((2, 3)))['class'] == 'no difference'
    blocked = localization_class(.01, s3, pairs(s3), [[0, 0, .3], [0, 0, 0]])  # relabeled ratio 0.3
    assert not blocked['specific'] and blocked['class'] == 'not resolved'
    assert localization_class(.3, s3, pairs(s3), np.zeros((2, 3)))['class'] == 'not resolved'


def test_replication_and_transfer_classes():
    pairs = np.array([.3, .4, .5, .6])
    assert replication_class(.01, True, .4, pairs, .9) == 'replicated'
    assert replication_class(.01, False, .4, pairs, .9) == 'direction only'
    assert replication_class(.2, True, .4, pairs, .9) == 'direction only'
    assert replication_class(.01, True, .2, pairs, .9) == 'direction only'     # below the 0.25 effect
    assert replication_class(.9, True, -.4, -pairs, .01) == 'opposite'
    assert replication_class(.9, True, .1, np.array([.1, -.1, .2, .2]), .5) == 'not replicated'
    assert transfer_class(.01, [.3, .5]) == 'zonated'
    assert transfer_class(.2, [.3, .5]) == 'attenuated'
    assert transfer_class(.01, [.1, .2]) == 'absent'
    assert transfer_class(.01, [.3, -.5]) == 'absent'
    assert transfer_class(.9, [-.3, -.5]) == 'reversed'
    assert informative('replicated', 'localized to S3') == 'position adds information'
    assert informative('not replicated', 'not resolved') == 'not informative'


def test_segment_summary_min_members_and_stat():
    frame = pd.DataFrame({'S1': [1., 2, 3, 4, 100], 'S2': [1., np.nan, np.nan, np.nan, 1]}, index=list('abcde'))
    out = segment_summary({'condition': frame}, {'p': list('abcde')}, min_members=5).set_index('scope')
    assert out.loc['S1', 'value'] == 3 and np.isnan(out.loc['S2', 'value'])
    mean = segment_summary({'condition': frame}, {'p': list('abcde')}, min_members=5, stat='mean').set_index('scope')
    assert mean.loc['S1', 'value'] == pytest.approx(22)


def test_read_programs_checks(tmp_path):
    path = tmp_path / 'programs.csv'
    pd.DataFrame({'program_id': ['A'], 'genes': ['g1, g2']}).to_csv(path, index=False)
    programs, _ = read_programs(path, ['g1', 'g2', 'g3'], sha256=file_sha256(path))
    assert programs == {'A': ['g1', 'g2']}
    with pytest.raises(ValueError):
        read_programs(path, ['g1', 'g2'], sha256='0' * 64)
    with pytest.raises(ValueError):
        read_programs(path, ['g1'])
    with pytest.raises(ValueError):
        read_programs(path, ['g1', 'g2'], exclude=['g2'])


def test_gene_span_and_stratified_rank_normal(tmp_path):
    fasta = tmp_path / 'x.fa.gz'
    with gzip.open(fasta, 'wt') as handle:
        handle.write('>T1 cdna chromosome:GRCh38:1:100:200:1 gene:G1 gene_biotype:x transcript_biotype:x gene_symbol:AAA1 d\n'
                     'ACGT\n'
                     '>T2 cdna chromosome:GRCh38:1:50:150:1 gene:G1 gene_biotype:x transcript_biotype:x gene_symbol:AAA1 d\n'
                     '>T3 cdna scaffold:GRCh38:KI1.1:10:19:-1 gene:G2 gene_biotype:x transcript_biotype:x gene_symbol:Bbb2 d\n')
    span = gene_span(fasta)
    assert span['AAA1'] == 151 and span['BBB2'] == 10
    assert gene_span(fasta, ['aaa1', 'zzz']).tolist()[0] == 151
    stat = pd.Series([1., 2, 3, 10, 20, 30], index=list('abcdef'))
    strata = pd.Series([0, 0, 0, 1, 1, 1], index=stat.index)
    out = stratified_rank_normal(stat, strata)
    np.testing.assert_allclose(out[list('abc')].to_numpy(), out[list('def')].to_numpy())


def test_bin_rates_and_added_information():
    counts = sparse.csr_matrix(np.arange(40, dtype=float).reshape(20, 2))
    rates, majority, names = bin_rates(counts, np.full(20, 10.), np.linspace(0, 1, 20), ['A'] * 10 + ['B'] * 10,
                                       ['S1'] * 5 + ['S3'] * 15, [0], n_bins=2)
    assert rates.shape == (2, 2, 1) and names == ['A', 'B'] and majority.loc['A', 0] == 'S1'
    line = added_information({'source': 'X', 'claim': 'up in S3', 'replication': 'replicated', 'scope': 'S3',
                              'value': .5, 'pairs': 4, 'q': .01, 'localization': 'localized to S3', 'peak': 'S3',
                              'unresolved': ['S1'], 'axis': 'agrees', 'flags': ['probe-sensitive']})
    assert 'S1: no difference resolved' in line and 'Axis agrees' in line and 'probe-sensitive' in line
