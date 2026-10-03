"""Small synthetic contract checks for repeated-section simulation."""
import numpy as np
import pytest

from pseudospace.repeated_structure_simulation import simulate_repeated_sections


def test_simulation_is_reproducible_balanced_and_returns_unordered_partial_axis():
    kwargs = dict(seed=26, n_per_specimen=40, specimens=('A', 'B', 'C'),
                  interval=(.2, .8), exposure_multiplier=1.7,
                  nuisance_mode='coupled')
    first = simulate_repeated_sections(**kwargs)
    second = simulate_repeated_sections(**kwargs)
    for key in ('counts', 'library', 'specimens', 'anatomy', 'gene_names',
                'structure_ids', 'true_z', 'nuisance_q'):
        np.testing.assert_array_equal(first[key], second[key])
    assert first['counts'].shape == (120, 40)
    assert np.all(first['true_z'] >= .2) and np.all(first['true_z'] <= .8)
    np.testing.assert_array_equal(first['counts'].sum(axis=1), first['library'])
    assert all(np.sum(first['specimens'] == sample) == 40 for sample in ('A', 'B', 'C'))
    assert len(np.unique(first['structure_ids'])) == len(first['structure_ids'])
    assert np.corrcoef(first['true_z'], first['nuisance_q'])[0, 1] == pytest.approx(1., abs=1e-12)


def test_independent_and_reversed_nuisance_modes():
    independent = simulate_repeated_sections(seed=3, n_per_specimen=150,
                                             nuisance_mode='independent')
    assert abs(np.corrcoef(independent['true_z'], independent['nuisance_q'])[0, 1]) < .15
    reversed_data = simulate_repeated_sections(seed=3, n_per_specimen=80,
                                               nuisance_mode='reversed')
    assert np.corrcoef(reversed_data['true_z'], reversed_data['nuisance_q'])[0, 1] < -.999


@pytest.mark.parametrize('kwargs', [
    {'n_per_specimen': 19},
    {'specimens': ()},
    {'specimens': ('A', 'A')},
    {'interval': (-.1, .8)},
    {'interval': (.8, .2)},
    {'spatial_strength': 0.},
    {'nuisance_strength': -1.},
    {'exposure_multiplier': 0.},
    {'nuisance_mode': 'unknown'},
])
def test_simulator_rejects_invalid_configuration(kwargs):
    with pytest.raises(ValueError):
        simulate_repeated_sections(seed=1, **kwargs)
