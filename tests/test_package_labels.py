"""package_labels on a synthetic package run: notebook 73's label-based inputs."""
import json

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from pseudospace.package_labels import (all_units, check_fingerprint, cortical_depth, harmony_features, label_panel_genes,
                                        label_units, native_sources)

SAMPLES = {'Ctrl1A2': 'Mouse', 'Ctrl1A4': 'Mouse', 'HUK1_COR1': 'Human', 'HUK1_MED1': 'Human'}
GENES = ['Gene0', 'Gene1', 'Gene2', 'Gene3', 'mt-Co1', 'mt-Nd1']
ORTHOLOGS = pd.DataFrame({'mouse_symbol': GENES + ['Gene9'], 'human_symbol': [g.upper() for g in GENES] + ['GENE9'],
                          'mapping_status': 'hcop_one_to_one'})
LABELS = {'0': 'PT-S1', '1': 'PT-S2', '2': 'PT-S3', '3': 'Podocyte'}
FINGERPRINT = {'n_cells': 48, 'n_clusters': 4, 'membership_sha1': 'abc123'}
PANELS = {'PT-S1': ['Gene0'], 'PT-S2': ['Gene1'], 'PT-S3': ['Gene2', 'Gene3'], 'Podocyte': ['Gene9']}


def _run(tmp_path, seed=0):
    """A package run with 12 units per sample: 9 PT (3 per segment) and 3 glomerulus units."""
    ad = pytest.importorskip('anndata')
    rng = np.random.default_rng(seed)
    results = tmp_path / 'results'
    export, run = results / 'package_pseudospace' / 'hvm', results / 'package_pseudospace' / 'runs' / 'hvm'
    for folder in (export, run / '01_units', run / '02_qc', run / '03_integration'):
        folder.mkdir(parents=True)
    obs, rows = [], []
    for sample, group in SAMPLES.items():
        native_genes = GENES + ['Extra'] if group == 'Mouse' else [g.upper() for g in GENES] + ['EXTRA']
        native = rng.poisson(4., (12, len(native_genes))).astype(float) + 1
        ad.AnnData(X=sparse.csr_matrix(native), obs=pd.DataFrame(index=[f'unit_{i}' for i in range(12)]),
                   var=pd.DataFrame({'gene_ids': [f'ID_{g}' for g in native_genes]}, index=native_genes)).write_h5ad(
            run / '01_units' / f'{sample}_units.h5ad')
        for i in range(12):
            label = ['PT-S1', 'PT-S2', 'PT-S3', 'Podocyte'][i // 3]
            obs.append({'unit': f'{sample}_unit_{i}', 'sample': sample, 'group': group, 'label': label,
                        'family': 'PT' if label.startswith('PT') else 'Glomerulus', 'feature_index': i,
                        'x_centroid': float(i), 'y_centroid': 0., 'n_spots': 10 + i, 'seg_class': 'tubule'})
            rows.append(native[i, :len(GENES)])
    obs = pd.DataFrame(obs).set_index('unit')
    counts = np.array(rows)
    obs['total_counts'] = counts.sum(axis=1)
    qc = ad.AnnData(X=sparse.csr_matrix(counts), obs=obs[['sample', 'feature_index', 'x_centroid', 'y_centroid', 'n_spots',
                                                            'total_counts', 'seg_class']], var=pd.DataFrame(index=GENES))
    qc.layers['counts'] = sparse.csr_matrix(counts)
    qc.write_h5ad(run / '02_qc' / 'units_qc.h5ad')
    obs[['sample', 'label', 'family']].to_csv(export / 'all_units_positions.csv')
    pt = obs.family.eq('PT').to_numpy()
    units = ad.AnnData(
        X=sparse.csr_matrix(counts[pt].astype(np.int32)),
        obs=pd.DataFrame({'sample': obs['sample'][pt], 'group': obs.group[pt], 'segment_class': obs.label[pt],
                          'pseudospace': np.linspace(0, 1, pt.sum()), 'feature_index': obs.feature_index[pt],
                          'x_centroid': obs.x_centroid[pt], 'y_centroid': obs.y_centroid[pt],
                          'n_spots': obs.n_spots[pt], 'total_counts': obs.total_counts[pt]}, index=obs.index[pt]),
        var=pd.DataFrame({'human_symbol': [g.upper() for g in GENES]}, index=GENES))
    units.uns['provenance'] = {'fingerprints': {'pass1': FINGERPRINT},
                               'config': {'marker_panels_yaml': {'label_panels': PANELS}}}
    units.write_h5ad(export / 'pt_units.h5ad')
    (run / '03_integration' / 'clustering_fingerprints.json').write_text(json.dumps(
        {'pass1': {**FINGERPRINT, 'n_features': 5}, 'labels_pass1': LABELS}))
    pd.DataFrame({'hvg_condition_count': [2, 2]}, index=['Gene1', 'Gene2']).to_csv(
        run / '03_integration' / 'pass1_harmony_features.csv')
    ad.AnnData(X=sparse.csr_matrix((2, len(GENES))), var=pd.DataFrame(
        {'highly_variable_for_harmony': [False, False, True, True, False, False]}, index=GENES)).write_h5ad(
        run / '03_integration' / 'pass2.h5ad')
    return results, obs, counts


def test_label_units_gene_space_library_and_lognorm(tmp_path):
    results, obs, counts = _run(tmp_path)
    pt = obs.family.eq('PT').to_numpy()
    full = label_units(results, 'hvm')
    assert list(full['gene_space']) == GENES and full['excluded'] == []
    np.testing.assert_allclose(full['counts'].toarray(), counts[pt])
    mt_free = label_units(results, 'hvm', exclude_genes=['mt-Co1', 'mt-Nd1', 'mt-Absent'])
    assert list(mt_free['gene_space']) == GENES[:4] and mt_free['excluded'] == ['mt-Co1', 'mt-Nd1']
    np.testing.assert_allclose(mt_free['library'], counts[pt][:, :4].sum(axis=1))
    np.testing.assert_allclose(mt_free['lognorm'].toarray(),
                               np.log1p(counts[pt][:, :4] / counts[pt][:, :4].sum(axis=1, keepdims=True) * 1e4))
    assert set(full['obs'].species) == {'mouse', 'human'} and set(full['obs'].segment) == {'S1', 'S2', 'S3'}
    assert full['obs'].loc['HUK1_COR1_unit_0', 'species'] == 'human'
    assert full['fingerprint'] == FINGERPRINT
    ordered = label_units(results, 'hvm', gene_order=['Gene3', 'Absent', 'Gene1'])
    assert list(ordered['gene_space']) == ['Gene3', 'Gene1', 'Gene0', 'Gene2', 'mt-Co1', 'mt-Nd1']
    np.testing.assert_allclose(ordered['counts'].toarray(), counts[pt][:, [3, 1, 0, 2, 4, 5]])


def test_all_units_joins_counts_and_library(tmp_path):
    results, obs, counts = _run(tmp_path)
    table, (marker_counts, present), library = all_units(results, 'hvm', genes=['Gene3', 'Absent', 'Gene0'],
                                                         exclude_genes=['mt-Co1', 'mt-Nd1'])
    assert list(table.index) == list(obs.index) and list(present) == ['Gene3', 'Gene0']
    np.testing.assert_allclose(marker_counts.toarray(), counts[:, [3, 0]])
    np.testing.assert_allclose(library, counts[:, :4].sum(axis=1))
    assert table.loc['Ctrl1A2_unit_10', 'family'] == 'Glomerulus' and table.feature_index.tolist()[:3] == [0, 1, 2]
    positions = pd.read_csv(results / 'package_pseudospace' / 'hvm' / 'all_units_positions.csv', index_col=0)
    pd.concat([positions, pd.DataFrame({'sample': ['Ctrl1A2'], 'label': ['x'], 'family': ['x']}, index=['ghost'])]).to_csv(
        results / 'package_pseudospace' / 'hvm' / 'all_units_positions.csv')
    with pytest.raises(ValueError, match='missing from units_qc'):
        all_units(results, 'hvm')


def test_cortical_depth_mean_of_k_nearest():
    units = pd.DataFrame({'sample': ['a', 'a', 'b'], 'x_centroid': [0., 10., 0.], 'y_centroid': [0., 0., 0.]})
    anchors = pd.DataFrame({'sample': ['a'] * 4, 'x_centroid': [1., 2., 3., 100.], 'y_centroid': [0.] * 4})
    depth = cortical_depth(units, anchors, ['a'], k=3)
    np.testing.assert_allclose(depth[:2], [2., 8.])
    assert np.isnan(depth[2])
    with pytest.raises(ValueError, match='fewer than 3'):
        cortical_depth(units, anchors.iloc[:2], ['a'], k=3)


def test_native_route_reproduces_the_export(tmp_path):
    """rebuild_pt_expression on the stage-01 matrices with the paper map gives the export's counts (guard G5)."""
    from pseudospace.pathway_inputs import rebuild_pt_expression

    results, obs, counts = _run(tmp_path)
    units = label_units(results, 'hvm')
    pt_obs, sources = native_sources(results, 'hvm', units['obs'])
    rebuilt = rebuild_pt_expression(pt_obs, sources, ORTHOLOGS)
    measured = rebuilt.var.measured_in_both_inputs.to_numpy()
    assert list(rebuilt.var_names[measured]) == GENES
    np.testing.assert_allclose(rebuilt.layers['counts'][:, measured].toarray(), units['counts'].toarray())


def test_features_panels_and_fingerprint(tmp_path):
    results, _, _ = _run(tmp_path)
    units = label_units(results, 'hvm')
    assert harmony_features(results, 'hvm') == {'Gene1', 'Gene2'}
    assert harmony_features(results, 'hvm', 'pass2') == {'Gene2', 'Gene3'}
    assert label_panel_genes(units['provenance']) == {'Gene0', 'Gene1', 'Gene2', 'Gene3'}
    confirmed = {**FINGERPRINT, 'labels_pass1': LABELS}
    assert check_fingerprint(results, 'hvm', units['provenance'], confirmed)['labels_pass1'] == LABELS
    for key, value in (('membership_sha1', 'other'), ('n_cells', 47), ('labels_pass1', {**LABELS, '3': 'PT-S1'})):
        with pytest.raises(ValueError, match=key):
            check_fingerprint(results, 'hvm', units['provenance'], {**confirmed, key: value})
    with pytest.raises(ValueError, match='not signed off'):
        check_fingerprint(results, 'hvm', units['provenance'], None)
