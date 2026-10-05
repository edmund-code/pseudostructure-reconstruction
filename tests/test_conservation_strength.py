import numpy as np
import pandas as pd
import pytest

from pseudospace.conservation_strength import (binned_statistics, block_draws, coexpression_modules, donor_folds,
                                               gene_draws, percentile_interval, percentile_rank, resampled_gradient,
                                               strength_bins, strength_score, subset_statistics)
from pseudospace.zonation_reliability import conservation_index, paired_gradient


def _replicates(truth, noise, n, rng, distortion=0.):
    systematic = distortion * rng.normal(0, np.std(truth), len(truth))
    return truth[None, :] + systematic[None, :] + rng.normal(0, noise, (n, len(truth)))


def _gradient(reps, genes):
    mean, var = paired_gradient(pd.DataFrame(reps, columns=genes))
    return mean, var


def test_estimator_reproduces_notebook43_index_and_interval():
    rng = np.random.default_rng(0)
    genes = pd.Index([f'g{i}' for i in range(800)])
    shared = rng.standard_t(3, 800)
    series = {k: _gradient(_replicates(shared * w, .6, 3, rng, .5), genes)
              for k, w in (('h', .5), ('m', 1.), ('h2', .6), ('m2', 1.1))}
    reference = conservation_index(series['h'], series['m'], series['h2'], series['m2'], genes, n_boot=40, seed=7, q=.01)
    data = {key: (series[k][0].to_numpy(), series[k][1].to_numpy())
            for key, k in (('human', 'h'), ('mouse', 'm'), ('human_ref', 'h2'), ('mouse_ref', 'm2'))}
    point = subset_statistics(data, np.arange(800), q=.01, match_reliability=False)
    draws = [subset_statistics(data, idx, q=.01, match_reliability=False)['index'] for idx in gene_draws(800, 40, 7)]
    low, high, _ = percentile_interval(np.array(draws)[:, None])
    assert point['index'] == pytest.approx(reference['index'], abs=1e-12)
    assert point['rstar_ceiling_human'] == pytest.approx(reference['ceiling_human'], abs=1e-12)
    assert low[0] == pytest.approx(reference['ci_low'], abs=1e-12) and high[0] == pytest.approx(reference['ci_high'], abs=1e-12)
    # Matched reliability differs only through the winsorised variance; with q = 0 the conventions agree.
    assert subset_statistics(data, np.arange(800), q=0)['index'] == pytest.approx(
        subset_statistics(data, np.arange(800), q=0, match_reliability=False)['index'])


def test_strength_bins_separate_shared_strong_genes_from_species_specific_weak_genes():
    """Strong genes share their gradient; weak genes are reproducible within species but differ."""
    rng = np.random.default_rng(1)
    n = 4000
    strong = np.arange(n) >= n - 400
    shared = np.where(strong, rng.normal(0, 3, n), 0)
    human_own, mouse_own = rng.normal(0, .5, n), rng.normal(0, .5, n)
    truth_h, truth_m = shared + np.where(strong, 0, human_own), shared + np.where(strong, 0, mouse_own)
    data = {k: tuple(a.to_numpy() for a in _gradient(_replicates(t, .3, 4, rng), range(n)))
            for k, t in (('human', truth_h), ('mouse', truth_m), ('human_ref', truth_h), ('mouse_ref', truth_m))}
    score = strength_score([percentile_rank(truth_h), percentile_rank(truth_m)])
    bins = strength_bins(score, (0, .9, 1))
    out = binned_statistics([(data, bins)], np.arange(n), {'weak': [0], 'strong': [1]}, q=0, min_reliability=.2)
    weak, top = out
    from pseudospace.conservation_strength import STATS
    k = {s: i for i, s in enumerate(STATS)}
    assert top[k['index']] > .9 and abs(weak[k['index']]) < .15
    assert weak[k['rstar_ceiling_human']] > .9 and weak[k['rstar_ceiling_mouse']] > .9


def test_min_reliability_marks_noise_only_subsets_not_estimable():
    rng = np.random.default_rng(2)
    n = 500
    data = {k: tuple(a.to_numpy() for a in _gradient(rng.normal(0, 1, (3, n)), range(n))) for k in ('human', 'mouse', 'human_ref', 'mouse_ref')}
    out = subset_statistics(data, np.arange(n), q=0, min_reliability=.2)
    assert np.isnan(out['rstar_cross']) and np.isnan(out['index'])


def test_bins_folds_blocks_and_donor_resampling():
    bins = strength_bins(np.linspace(0, 1, 1000))
    assert np.bincount(bins).tolist() == [500, 300, 100, 50, 40, 10]
    folds = donor_folds([f'd{i}' for i in range(7)])
    assert sorted(np.bincount(folds).tolist()) == [3, 4]
    assert (donor_folds([f'd{i}' for i in range(7)][::-1]) == folds[::-1]).all()
    blocks = np.repeat(np.arange(5), [1, 2, 3, 4, 5])
    for draw in block_draws(blocks, 20, 0):
        counts = pd.Series(blocks[draw]).value_counts()
        assert all(counts[b] % (b + 1) == 0 for b in counts.index)
    per_donor = np.arange(12.).reshape(4, 3)
    mean, noise = resampled_gradient(per_donor, np.random.default_rng(0))
    assert mean.shape == (3,) and (noise >= 0).all()


def test_coexpression_modules_recover_planted_programs():
    rng = np.random.default_rng(3)
    cells, per = 600, 20
    factors = rng.normal(0, 1, (cells, 3))
    matrix = np.hstack([factors[:, [j]] + .3 * rng.normal(0, 1, (cells, per)) for j in range(3)])
    specimen = np.repeat(['a', 'b', 'c', 'd'], cells // 4)
    species = np.where(np.isin(specimen, ['a', 'b']), 'mouse', 'human')
    labels = coexpression_modules(matrix, specimen, species, n_components=5, n_modules=3, seed=0)
    assert all(len(set(labels[j * per:(j + 1) * per])) == 1 for j in range(3))
    assert len(set(labels)) == 3
