import json

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

ad = pytest.importorskip('anndata')
pytest.importorskip('pseudospace_reconstruction')

from pseudospace.package_refits import (cleanup_refit, orient_like_export, perturbed_units, refit_position,
                                        runs_s1_to_s3, write_refit_input)


def _base(n=40, g=12, seed=0):
    from pseudospace_reconstruction.qc import normalize_log1p

    rng = np.random.default_rng(seed)
    counts = sparse.csr_matrix(rng.poisson(3, (n, g)).astype(np.float32))
    obs = pd.DataFrame({'sample': np.repeat(['a', 'b'], n // 2), 'condition': 'x'},
                       index=[f'u{i}' for i in range(n)])
    var = pd.DataFrame({'mt': False, 'exclude_from_harmony_hvg': [True] + [False] * (g - 1)},
                       index=[f'G{j}' for j in range(g)])
    base = ad.AnnData(X=counts.copy(), obs=obs, var=var)
    return normalize_log1p(base, 1e4)


def test_unperturbed_copy_reproduces_the_package_normalisation():
    base = _base()
    out = perturbed_units(base)
    assert out.obs_names.equals(base.obs_names) and out.var_names.equals(base.var_names)
    assert abs(out.layers['lognorm'] - base.layers['lognorm']).max() == 0
    assert abs(out.X - base.layers['lognorm']).max() == 0
    assert out.layers['counts'].dtype == base.layers['counts'].dtype


def test_floor_differs_only_at_round_off_and_counts_change_lognorm():
    base = _base()
    floor = perturbed_units(base, multiply=True)
    assert abs(floor.layers['lognorm'] - base.layers['lognorm']).max() < 1e-5
    thinned = base.layers['counts'].copy()
    thinned.data = np.floor(thinned.data / 2)
    out = perturbed_units(base, counts=thinned)
    assert abs(out.layers['counts'] - thinned).max() == 0
    totals = np.asarray(np.expm1(out.layers['lognorm']).sum(axis=1)).ravel()
    np.testing.assert_allclose(totals[totals > 0], 1e4, rtol=1e-3)


def test_fold_genes_are_flagged_and_unknown_genes_refused():
    base = _base()
    out = perturbed_units(base, exclude_genes=['G3', 'G5'])
    assert out.var.exclude_from_harmony_hvg.tolist()[:6] == [True, False, False, True, False, True]
    assert not base.var.exclude_from_harmony_hvg.iloc[3]
    with pytest.raises(ValueError):
        perturbed_units(base, exclude_genes=['nope'])
    with pytest.raises(ValueError):
        perturbed_units(base, counts=-base.layers['counts'])


def test_write_input_keeps_ids_and_records_the_perturbation(tmp_path):
    base_run = tmp_path / 'base'
    (base_run / '02_qc').mkdir(parents=True)
    (base_run / '02_qc' / 'manifest.json').write_text(json.dumps({'stage': '02', 'summary': {'n': 1}}))
    refit = write_refit_input(perturbed_units(_base()), base_run, tmp_path / 'refit', {'name': 'identity'})
    back = ad.read_h5ad(refit / '02_qc' / 'units_qc.h5ad')
    assert back.obs_names.tolist() == _base().obs_names.tolist()
    assert json.loads((refit / '02_qc' / 'manifest.json').read_text())['notebook72_perturbation'] == {'name': 'identity'}
    assert cleanup_refit(refit) > 0 and not (refit / '02_qc' / 'units_qc.h5ad').exists()


def test_extraction_payload_loads_without_pickle(tmp_path):
    from pseudospace.package_refits import extract_refit

    run = tmp_path / 'refit'
    (run / '03_integration').mkdir(parents=True)
    (run / '04_pseudospace').mkdir(parents=True)
    table = pd.DataFrame({'sample': pd.Categorical(['a', 'b']), 'segment_class': pd.Categorical(['PT-S1', 'PT-S3']),
                          'coarse_class': ['PT', 'PT'], 'pseudospace': [.1, .9], 'dpt_pseudotime': [.2, np.nan]},
                         index=pd.Index(['u0', 'u2'], name='obs_name'))
    table.to_parquet(run / '04_pseudospace' / 'pseudospace.parquet')
    base = _base(n=4, g=3)
    base.obs['segment_class'] = pd.Categorical(['PT-S1', 'PT-S2', 'PT-S3', 'PT-S1'])
    base.obs['coarse_class'] = 'PT'
    base.var['highly_variable_for_harmony'] = [False, True, True]
    base.write_h5ad(run / '03_integration' / 'pass1_labelled.h5ad')
    base[['u0', 'u2']].write_h5ad(run / '03_integration' / 'pass2.h5ad')
    for folder in ('03_integration', '04_pseudospace'):
        (run / folder / 'manifest.json').write_text(json.dumps({'summary': {'n': 2}}))
    payload = extract_refit(run)
    np.savez_compressed(tmp_path / 'p.npz', **payload)
    with np.load(tmp_path / 'p.npz', allow_pickle=False) as stored:
        back = {k: stored[k] for k in stored.files}
    assert back['unit'].tolist() == ['u0', 'u2'] and back['segment_class'].tolist() == ['PT-S1', 'PT-S3']
    assert back['pass1_features'].tolist() == ['G1', 'G2'] and len(back['pass1_unit']) == 4
    np.testing.assert_array_equal(back['dpt_pseudotime'], [.2, np.nan])


def test_matching_back_and_orientation_guards():
    payload = {'unit': np.array(['u2', 'u0']), 'pseudospace': np.array([.7, .1])}
    np.testing.assert_array_equal(refit_position(payload, ['u0', 'u1', 'u2']), [.1, np.nan, .7])
    specimen = np.array(['a'] * 4 + ['b'] * 4)
    segment = np.array(['PT-S1', 'PT-S1', 'PT-S3', 'PT-S3'] * 2)
    good = np.array([.1, .2, .8, .9, .1, .3, .7, .9])
    assert runs_s1_to_s3(good, specimen, segment)
    bad = good.copy()
    bad[4:] = bad[4:][::-1]
    assert not runs_s1_to_s3(bad, specimen, segment)
    missing = good.copy()
    missing[4:6] = np.nan
    assert not runs_s1_to_s3(missing, specimen, segment)
    flipped, was_flipped = orient_like_export(1 - good, segment)
    assert was_flipped and np.allclose(flipped, good)
