import numpy as np
import pandas as pd

from pseudospace.zonation_classes import (confirm, cross_species_class, moderated_contrast, union_class,
                                          zonation_calls)


def _replicates(true, noise, n=2, seed=0):
    rng = np.random.default_rng(seed)
    genes = [f'g{i}' for i in range(len(true))]
    return pd.DataFrame(np.asarray(true)[None, :] + rng.normal(0, noise, (n, len(true))), columns=genes)


def test_moderated_contrast_is_calibrated_under_null():
    reps = _replicates(np.zeros(4000), .3, n=3)
    out = moderated_contrast(reps, pd.Series(np.arange(4000), index=reps.columns))
    assert abs((out.p <= .05).mean() - .05) < .015


def test_calls_distinguish_flat_from_indeterminate_by_noise():
    true = np.r_[np.full(200, 1.5), np.zeros(1800)]
    quiet = moderated_contrast(_replicates(true, .05, n=3, seed=1), pd.Series(np.zeros(2000), index=[f'g{i}' for i in range(2000)]))
    noisy = moderated_contrast(_replicates(true, 1.2, n=3, seed=2), pd.Series(np.zeros(2000), index=[f'g{i}' for i in range(2000)]))
    quiet_calls, _ = zonation_calls(quiet)
    noisy_calls, _ = zonation_calls(noisy)
    assert (quiet_calls.iloc[:200] == 'up').mean() > .95 and (quiet_calls.iloc[200:] == 'flat').mean() > .9
    assert (noisy_calls.iloc[200:] == 'flat').mean() < .2          # a noisy null is not "flat"
    assert (noisy_calls.iloc[200:] == 'indeterminate').mean() > .7


def test_confirmation_drops_unsupported_calls():
    calls = pd.Series(['up', 'down', 'up', 'flat'], index=list('abcd'))
    reference = pd.DataFrame({'mean': [1, 1, 0, 0], 't': [5., 5., 0., 0.], 'df': [10.] * 4}, index=list('abcd'))
    out = confirm(calls, reference)
    assert out.tolist() == ['up', 'indeterminate', 'indeterminate', 'flat']


def test_cross_species_classes_and_union():
    assert cross_species_class('up', 'up') == 'conserved'
    assert cross_species_class('up', 'down') == 'reversal'
    assert cross_species_class('down', 'flat') == 'mouse-only'
    assert cross_species_class('down', 'indeterminate') == 'indeterminate'
    assert cross_species_class('flat', 'up') == 'human-only'
    assert cross_species_class('flat', 'indeterminate') == 'neither'
    assert cross_species_class('indeterminate', 'indeterminate') == 'indeterminate'
    table = pd.DataFrame({'a': ['mouse-only', 'neither', np.nan], 'b': ['reversal', 'human-only', np.nan]})
    assert union_class(table).iloc[:2].tolist() == ['reversal', 'human-only'] and pd.isna(union_class(table).iloc[2])
