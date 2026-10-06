import pandas as pd
import pytest

from pseudospace import ortholog_pairing as op


def _tables():
    mouse = pd.DataFrame({
        'component': [1, 1, 1, 2, 2, -1, 3, 3],
        'ortholog_class': ['1:many', '1:many', '1:many', 'many:many', 'many:many', '1:1', '1:many', '1:many'],
        'pt_expressed': [False, True, True, True, False, True, False, False],
        'probe_clean': [True, True, False, True, True, True, True, True],
        'log2_cp10k': [1., 2., 5., 3., 0., 1., -1., 0.],
    }, index=['Ma', 'Mb', 'Mc', 'Md', 'Me', 'Mx', 'Mf', 'Mg'])
    human = pd.DataFrame({
        'component': [1, 2, 2, -1, 3],
        'ortholog_class': ['1:many', 'many:many', 'many:many', '1:1', '1:many'],
        'pt_expressed': [True, True, True, True, True],
        'probe_clean': [True, True, True, True, True],
        'log2_cp10k': [1., 4., 2., 1., 0.],
    }, index=['HA', 'HD', 'HE', 'HX', 'HF'])
    pairs = pd.DataFrame({'human_symbol': ['HA', 'HD', 'HE', 'HX', 'HF'],
                          'mouse_symbol': ['Ma', 'Me', 'Md', 'Mx', 'Mf'],
                          'mapping_status': ['hcop_one_to_one'] * 5})
    return pairs, mouse, human


def test_swap_uses_main_probe_clean_pt_paralog():
    pairs, mouse, human = _tables()
    swapped, log = op.swap_pairs(pairs, mouse, human)
    as_dict = dict(zip(swapped.human_symbol, swapped.mouse_symbol))
    # family 1: Mc is most expressed but probe-ambiguous, so Mb (clean, PT-expressed) replaces Ma
    assert as_dict['HA'] == 'Mb'
    # family 2: both main paralogs (HD, Md) are already kept, so its two pairs stay as they are
    assert as_dict['HD'] == 'Me' and as_dict['HE'] == 'Md'
    # 1:1 pair untouched; family 3 has no PT-expressed mouse member, so the kept member stays
    assert as_dict['HX'] == 'Mx' and as_dict['HF'] == 'Mf'
    assert len(swapped) == len(pairs) and list(log.component) == [1]
    assert log.iloc[0]['mouse changed'] and not log.iloc[0]['human changed']


def test_swap_replaces_least_expressed_kept_member():
    pairs, mouse, human = _tables()
    mouse.loc['Mh'] = [2, 'many:many', True, True, 9.]
    swapped, log = op.swap_pairs(pairs, mouse, human)
    as_dict = dict(zip(swapped.human_symbol, swapped.mouse_symbol))
    assert as_dict['HD'] == 'Mh' and as_dict['HE'] == 'Md'   # Me (0) was the least expressed kept member


def test_drop_and_probe_balance():
    pairs, mouse, human = _tables()
    kept, dropped = op.drop_pairs(pairs, mouse, human)
    assert list(kept.human_symbol) == ['HX'] and len(dropped) == 4
    table = op.probe_balance(pairs, {'E1': 2, 'E2': 3}, {'H1': 2, 'H2': 1}, {'Ma': 'E1', 'Me': 'E2'}, {'HA': 'H1', 'HD': 'H2'})
    assert table.set_index('human_symbol').loc['HA', 'probe_balance'] == 'balanced'
    assert table.set_index('human_symbol').loc['HD', 'probe_balance'] == 'human_fewer'
    assert table.set_index('human_symbol').loc['HX', 'mouse_probes'] == 0


def test_pair_across_families_raises():
    pairs, mouse, human = _tables()
    bad = pairs.copy()
    bad.loc[0, 'mouse_symbol'] = 'Md'
    with pytest.raises(ValueError):
        op.family_of_pairs(bad, mouse, human)
