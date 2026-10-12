import numpy as np
import pandas as pd
import pytest

from pseudospace import species_calls as sc


def _classes():
    columns = ['S2-S1', 'S2-S1 mouse call', 'S2-S1 human call', 'S3c-early', 'S3c-early mouse call', 'S3c-early human call']
    rows = {
        'm1': ['mouse-only', 'up', 'flat', 'indeterminate', 'up', 'indeterminate'],
        'm2': ['indeterminate', 'down', 'indeterminate', np.nan, np.nan, np.nan],
        'm3': ['mouse-only', 'up', 'flat', 'mouse-only', 'up', 'flat'],
        'h1': ['human-only', 'flat', 'up', np.nan, np.nan, np.nan],
        's1': ['mouse-only', 'up', 'flat', 'conserved', 'up', 'up'],
        'x1': ['mouse-only', 'up', 'flat', 'human-only', 'flat', 'down'],
        'n1': ['flat', 'flat', 'flat', 'flat', 'flat', 'flat'],
    }
    return pd.DataFrame.from_dict(rows, orient='index', columns=columns)


def test_member_category_rules():
    classes = _classes()
    get = lambda g: sc.member_category(classes.loc[g])  # noqa: E731
    assert [get(g) for g in ('m1', 'm2', 'h1', 's1', 'x1', 'n1')] == [
        'mouse-led', 'mouse-led', 'human-led', 'shared', 'mixed', 'not zonated']
    assert sc.member_category(None) == 'not classifiable'


def test_call_categories_labels_and_unresolved():
    classes = _classes()
    sets = {'A': ['m1', 'm2', 'm3', 'h1', 'z'], 'B': ['m1', 'h1', 'n1', 'z'], 'C': ['s1', 'x1', 'm1', 'h1']}
    out = sc.call_categories(sets, list(sets), contributing={'m1', 'm2', 'm3', 'h1', 's1', 'x1', 'n1'}, classes=classes)
    assert out.loc['A', 'category'] == 'summary of mouse zonation'      # 3 of 4 zonated members mouse-led
    assert out.loc['B', 'category'] == 'unresolved'                     # only 2 zonated contributing members
    assert out.loc['C', 'category'] == 'mixed'                          # 1 shared, 1 mixed, 1 mouse, 1 human
    assert out.loc['A', 'n_not classifiable'] == 0 and out.loc['A', 'contributing_members'] == 4


def test_call_peaks_recovers_a_planted_segment_and_direction():
    rng = np.random.default_rng(0)
    genes = [f'g{i}' for i in range(200)]
    z = pd.DataFrame(rng.normal(size=(200, 3)), index=genes)
    z.iloc[:10, 1] = -6.0                               # members differ most in S2, mouse higher
    z.iloc[10, 0] = np.nan                              # a failed gene leaves the ranking
    out = sc.call_peaks(z, {'P': genes[:11]}, ['P'], np.array([0., 1., 2.]))
    assert out.loc['P', 'peak_position'] == 1.0 and not out.loc['P', 'human_high']
    assert out.loc['P', 'members_used'] == 10


def test_transitions_and_segment_units():
    t = pd.DataFrame({'specimen': ['a', 'a', 'b', 'b'], 'transition': ['S1→S2', 'S2→S3'] * 2,
                      'position': [.3, .6, .5, .4]})
    assert sc.mean_transitions(t, ['a']) == (.3, .6)
    assert sc.mean_transitions(t, ['b']) is None                       # crossed: undefined
    units = sc.to_segment_units([.15, .45, .8], (.3, .6))
    np.testing.assert_allclose(units, [.5, 1.5, 2.5])
    assert np.isnan(sc.to_segment_units([.2], None)).all()
    np.testing.assert_array_equal(sc.segment_index([.5, 1.99, 3.0, np.nan])[:3], [0, 1, 2])


def _fake_camera(statistic, sets):
    rows = []
    for p, members in sets.items():
        effect = sc.set_effect(statistic, members)
        rows.append({'pathway_id': p, 'direction': 'Up' if effect > 0 else 'Down', 'fdr': .01 if effect > 1 else .5})
    return pd.DataFrame(rows)


def test_subset_retention():
    genes = [f'g{i}' for i in range(100)]
    z = pd.Series(0., index=genes)
    z.iloc[:12] = 3.
    sets = {'A': genes[:12], 'B': genes[50:70]}
    full = sc.subset_retention(z, sets, ['A', 'B'], genes, camera=_fake_camera)
    assert full.loc['A', 'retained'] and full.loc['A', 'effect_ratio'] == pytest.approx(1)
    assert not full.loc['B', 'retained']
    dropped = sc.subset_retention(z, sets, ['A'], genes[5:], camera=_fake_camera)   # 7 members left: untestable
    assert not dropped.loc['A', 'testable'] and not dropped.loc['A', 'retained']


def test_ours_test_arms():
    claim = {'mouse': '-', 'human': '0', 'interaction': '+'}
    values = pd.Series({'M1': -2., 'M2': -1.5, 'H1': -.1, 'H2': .1})
    assert sc.ours_test(claim, values, ['H1', 'H2'], ['M1', 'M2']) is True
    values['M2'] = .2                                   # one mouse specimen breaks the signed arm
    assert sc.ours_test(claim, values, ['H1', 'H2'], ['M1', 'M2']) is False
    assert sc.ours_test({'mouse': '', 'human': '', 'interaction': ''}, values, ['H1', 'H2'], ['M1', 'M2']) is None


def test_replication_decision_and_per_call_counts():
    row = {'p_random, calls': .001, 'p calls vs noise': .01, 'p programs (label permutation)': .2}
    assert sc.replication_decision(row, 'calls', 'noise') == 'better than random only'
    row['p programs (label permutation)'] = .01
    assert sc.replication_decision(row, 'calls', 'noise') == 'specific'
    row['p_random, calls'] = .5
    assert sc.replication_decision(row, 'calls', 'noise') == 'not replicated'
    out = sc.per_call_replication(pd.Series({'a': 4., 'b': -3., 'c': .1}), ['a', 'b', 'c', 'missing'])
    assert out.replicated.tolist() == [True, False, False] and list(out.index) == ['a', 'b', 'c']


def test_replication_inputs_and_collection_run_on_synthetic_donors():
    rng = np.random.default_rng(1)
    genes = [f'g{i}' for i in range(300)]
    signal = np.r_[np.full(30, 2.), np.zeros(270)]
    ours = pd.DataFrame({'M1': signal + rng.normal(0, .2, 300), 'M2': signal + rng.normal(0, .2, 300),
                         'H1': rng.normal(0, .2, 300), 'H2': rng.normal(0, .2, 300)}, index=genes)
    cov = pd.DataFrame({'detection_mouse': .5, 'detection_human': .5, 'coverage': rng.uniform(.5, 1, 300)}, index=genes)
    base = pd.Series(rng.uniform(3, 8, 300), index=genes)
    donors = lambda n, shift: pd.DataFrame([base + shift + rng.normal(0, .1, 300) for _ in range(n)],  # noqa: E731
                                           index=[f'd{shift}{k}' for k in range(n)])
    h_late, h_early = donors(4, 0), donors(4, 0)
    m_early = donors(4, 0)
    m_late = m_early + signal + rng.normal(0, .1, (4, 300))
    sets = {'zonated': genes[:20], 'random': genes[100:130]}
    ref = sc.replication_inputs(ours, ['H1', 'H2'], ['M1', 'M2'], cov, sets, h_late, h_early, m_late, m_early, m_late.index)
    assert len(ref['genes']) == 300 and ref['component']['slope of interaction on mouse gradient'] < -.5
    run = sc.collection_run(ref, 'signed', ref['sets'], {'zonated': ['zonated'], 'random': ['random']}, ref['strata'],
                            n_perm=49, seed=0)
    z = pd.Series(run['set_z'], index=run['set_ids'])
    assert z['zonated'] > z['random']


def test_census_donors_keeps_complete_adult_donors():
    meta = pd.DataFrame({'development_stage': ['adult'] * 4 + ['embryo'] * 2,
                         'cell_type': ['early', 'late', 'early', 'late', 'early', 'late'],
                         'n_cells': [60, 60, 60, 10, 60, 60], 'donor_id': ['a', 'a', 'b', 'b', 'c', 'c'],
                         'sex': ['male'] * 6, 'total_counts_all_genes': [1e6] * 6}, index=list('uvwxyz'))
    table = pd.DataFrame({'G1': [10, 20, 10, 20, 10, 20]}, index=list('uvwxyz'))
    late, early, whole, info = sc.census_donors(table, meta, {'early': 'early', 'late': 'S3'}, ['early'], ['G1'])
    assert list(late.index) == ['a'] and list(info.index) == ['a']
    assert late.loc['a', 'G1'] == pytest.approx(np.log2(20 + 1))
