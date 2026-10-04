import numpy as np
import pandas as pd
import pytest

from pseudospace.zonation_reliability import (confirm_flat, conservation_index, disattenuated, paired_gradient,
                                              reliability, unpaired_gradient)


def _dataset(truth, noise, n, rng, distortion=0.):
    """Replicates of a gradient: truth plus a dataset-specific distortion plus replicate noise."""
    genes = truth.index
    systematic = distortion * rng.normal(0, truth.std(), len(truth))
    reps = truth.to_numpy()[None, :] + systematic[None, :] + rng.normal(0, noise, (n, len(truth)))
    return paired_gradient(pd.DataFrame(reps, columns=genes))


def test_reliability_and_disattenuation_recover_true_correlation():
    rng = np.random.default_rng(0)
    truth = pd.Series(rng.normal(0, 1, 5000), index=[f'g{i}' for i in range(5000)])
    a, b = _dataset(truth, 1.5, 3, rng), _dataset(truth, 1.5, 3, rng)
    out = disattenuated(a, b, truth.index, q=0)
    assert out['r'] < .7 and out['r_corrected'] == pytest.approx(1, abs=.08)
    assert reliability(*a) == pytest.approx(1 / (1 + 1.5 ** 2 / 3), abs=.05)


def test_conservation_index_separates_ceiling_from_divergence():
    rng = np.random.default_rng(1)
    genes = pd.Index([f'g{i}' for i in range(3000)])
    shared = pd.Series(rng.normal(0, 1, 3000), index=genes)
    human_specific = pd.Series(rng.normal(0, 1, 3000), index=genes)
    human_truth = .5 * shared + np.sqrt(.75) * human_specific        # true cross-species r = 0.5
    data = {k: _dataset(t, .5, 3, rng, distortion=.6) for k, t in
            (('h', human_truth), ('h2', human_truth), ('m', shared), ('m2', shared))}
    out = conservation_index(data['h'], data['m'], data['h2'], data['m2'], genes, n_boot=50, seed=2, q=0)
    assert out['ceiling_human'] < .9 and out['index'] == pytest.approx(.5, abs=.1)
    assert out['ci_low'] < out['index'] < out['ci_high']


def test_unpaired_gradient_variance():
    late = pd.DataFrame([[2., 0.], [4., 0.]], columns=['a', 'b'])
    early = pd.DataFrame([[1., 0.], [1., 0.], [1., 0.]], columns=['a', 'b'])
    mean, noise = unpaired_gradient(late, early)
    assert mean['a'] == 2 and noise['a'] == pytest.approx(1.)


def test_confirm_flat_drops_species_only_calls_the_reference_contradicts():
    classes = pd.Series(['mouse-only', 'mouse-only', 'human-only', 'conserved'], index=list('abcd'))
    calls_mouse = pd.Series(['up', 'up', 'flat', 'up'], index=list('abcd'))
    calls_human = pd.Series(['flat', 'flat', 'down', 'up'], index=list('abcd'))
    ref_human = pd.DataFrame({'mean': [.1, .8, 0, 0], 't': [.5, 6., 0, 0], 'df': [10.] * 4}, index=list('abcd'))
    ref_mouse = pd.DataFrame({'mean': [0, 0, -.9, 0], 't': [0, 0, -5., 0], 'df': [10.] * 4}, index=list('abcd'))
    out = confirm_flat(classes, calls_mouse, calls_human, ref_mouse, ref_human)
    assert out.tolist() == ['mouse-only', 'indeterminate', 'indeterminate', 'conserved']
