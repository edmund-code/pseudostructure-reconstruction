import numpy as np
import pandas as pd
import pytest

from pseudospace import gene_atlas as ga


def test_join_and_convergence_respect_the_testable_mask():
    genes = ['a', 'b', 'c', 'd']
    calls, testable = ga.join_lines(genes, {
        'L1': (['a', 'b'], genes),
        'L2': (['a', 'b'], ['a', 'b']),       # c and d were never looked at by L2
        'L3': (['a', 'c'], genes),
    })
    assert calls.loc['a'].tolist() == [True, True, True]
    assert testable.loc['c'].tolist() == [True, False, True]
    counts = ga.convergence_counts(calls, testable)
    assert counts[1] == 3 and counts[2] == 2 and counts[3] == 1       # a:3, b:2, c:1, d:0


def test_stratified_null_keeps_line_sizes_and_kills_agreement():
    rng = np.random.default_rng(0)
    n = 600
    genes = [f'g{i}' for i in range(n)]
    strata = pd.Series(np.repeat([0, 1], n // 2), index=genes)
    shared = rng.choice(n, 40, replace=False)                          # the same genes in three lines
    lines = {}
    for name in ('L1', 'L2', 'L3'):
        hit = set(shared) | set(rng.choice(n, 40, replace=False))
        lines[name] = ([genes[i] for i in sorted(hit)], genes)
    calls, testable = ga.join_lines(genes, lines)
    observed = ga.convergence_counts(calls, testable)
    null = ga.stratified_line_null(calls, testable, strata, draws=40, seed=1)
    assert observed[3] >= 40                                           # the planted core
    assert null[3].max() < observed[3]                                 # chance alone cannot make it
    assert null[2].max() < observed[2]
    total = int(calls.sum().sum())
    assert null[1].max() <= total and null[1].mean() > observed[1]      # calls preserved, spread wider


def test_similarity_and_yield_profile():
    genes = ['a', 'b', 'c', 'd', 'e']
    calls, testable = ga.join_lines(genes, {
        'M1': (['a', 'b', 'c'], genes),
        'M2': (['a', 'b'], genes),
        'X': (['a', 'd'], genes),
    })
    jac, phi = ga.line_similarity(calls, testable)
    assert jac.loc['M1', 'M2'] == pytest.approx(2 / 3) and jac.loc['M1', 'M1'] == 1
    assert phi.loc['M1', 'M2'] > phi.loc['M1', 'X']
    flag = lambda names: pd.Series(np.isin(genes, names), index=genes)
    prof = ga.yield_profile(calls, testable, ['M1', 'M2'],
                            clean=flag(['a', 'b', 'c', 'd']), supported=flag(['a', 'c', 'd']),
                            anchored=flag(['a']), min_lines=2).set_index('method')
    assert prof.loc['M1', 'called'] == 3 and prof.loc['M1', 'yield'] == 2      # a and c
    assert prof.loc['M1', 'novel yield'] == 0                                  # c has one line only
    assert prof.loc['M2', 'yield'] == 1 and prof.loc['M1', 'unique share'] == pytest.approx(1 / 3)


def test_evidence_model_weights_the_informative_line_and_ranks_held_out_genes():
    rng = np.random.default_rng(3)
    n = 1200
    informative = rng.random(n) < .3
    noise = rng.random(n) < .3
    truth = rng.random(n) < np.where(informative, .6, .1)            # only the first line carries the outcome
    X = pd.DataFrame({'informative': informative.astype(float), 'noise': noise.astype(float),
                      'constant': np.ones(n)})
    train, test = slice(0, 800), slice(800, None)
    coef, score = ga.evidence_model(X.iloc[train], truth[train], seed=0)
    assert coef['informative'] > 1.0 > abs(coef['noise'])
    assert coef['constant'] == 0                                      # constant column is dropped, not fitted
    held = score(X.iloc[test])
    assert held[truth[test]].mean() > held[~truth[test]].mean()
