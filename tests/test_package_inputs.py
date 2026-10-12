"""package_inputs on synthetic exports: notebook 69's input rules, rebuilt from a pt_units.h5ad."""
import json

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from pseudospace.coordinate_inputs import EXPECTED, LIBRARIES, common_support_mask, pt_inputs
from pseudospace.package_inputs import aki_inputs, hvm_inputs, read_pt_units, support_positions

MOUSE = [f'Gene{i}' for i in range(14)] + ['mt-Co1', 'mt-Atp8']
ORTHOLOGS = pd.DataFrame({'mouse_symbol': MOUSE, 'human_symbol': [g.upper() for g in MOUSE],
                          'mapping_status': 'hcop_one_to_one'})
MT = ['mt-Co1', 'mt-Atp8']
SETS = {'transport': ['GENE0', 'GENE1', 'GENE2', 'GENE3', 'MT-CO1'], 'fatty': ['GENE4', 'GENE5', 'GENE6', 'GENE12'],
        'rare': ['GENE7', 'GENE12', 'GENE13'], 'mito': ['GENE8', 'MT-CO1', 'MT-ATP8']}


def _segment(x):
    return np.where(x < 1 / 3, 'PT-S1', np.where(x < 2 / 3, 'PT-S2', 'PT-S3'))


def _gene_sets(data):
    folder = data / 'mouse_vs_human' / 'pathway_gene_sets'
    folder.mkdir(parents=True, exist_ok=True)
    for lib in LIBRARIES:
        (folder / f'{lib}.json').write_text(json.dumps(SETS))


def _hvm(tmp_path, n=60, seed=3, scale=(1., 0.)):
    """The old inputs of pt_inputs (tubule_by_gene, the saved obs, a coordinate) and the same data as an export."""
    ad = pytest.importorskip('anndata')
    rng = np.random.default_rng(seed)
    results, data = tmp_path / 'results', tmp_path / 'data'
    (results / 'minimal_pt_scfates').mkdir(parents=True)
    (data / 'tubule_by_gene').mkdir(parents=True)
    ORTHOLOGS.to_csv(results / 'minimal_pt_scfates' / 'ortholog_map_used.csv', index=False)
    _gene_sets(data)
    rates = np.r_[rng.gamma(2., 2., 12), .4, .002, 15., 5.]      # Gene13 is rare: below the detection rule
    obs, coordinate, rows, export_ids = [], {}, [], []
    for species, names in EXPECTED.items():
        for name in names:
            # Human panels do not measure mt-Atp8 (as in the real panels) nor Gene12 (a structural zero).
            keep = [g for g in MOUSE if species == 'mouse' or g not in ('mt-Atp8', 'Gene12')]
            idx = [MOUSE.index(g) for g in keep]
            counts = rng.poisson(rates[idx] * rng.gamma(4., .25, (n, 1))).astype(float)
            ids = [f's{i}' for i in range(n)]
            symbols = keep if species == 'mouse' else [g.upper() for g in keep]
            ad.AnnData(X=sparse.csr_matrix(counts), obs=pd.DataFrame(index=ids), var=pd.DataFrame(index=symbols)).write_h5ad(
                data / 'tubule_by_gene' / f'{name}_tubule_by_gene_caleb.h5ad')
            position = np.sort(rng.uniform(0, 1, n)) * scale[0] + scale[1]
            full = np.zeros((n, len(MOUSE)))
            full[:, idx] = counts
            for i, x, row in zip(ids, position, full):
                coordinate[f'{name}_{i}'] = x
                obs.append({'id': f'{name}_{i}', 'sample': name, 'comparison_species': species,
                            'broad_tubule_marker_call': 'PT', 'segment_class': str(_segment((x - scale[1]) / scale[0]))})
                rows.append(row)
    obs = pd.DataFrame(obs).set_index('id')
    ad.AnnData(X=sparse.csr_matrix((len(obs), 1)), obs=obs).write_h5ad(
        results / 'minimal_pt_scfates' / 'cross_species_pt_scfates.h5ad')
    shared = [g for g in MOUSE if g not in ('mt-Atp8', 'Gene12')]
    export = ad.AnnData(
        X=sparse.csr_matrix(np.array(rows)[:, [MOUSE.index(g) for g in shared]].astype(np.int32)),
        obs=pd.DataFrame({'sample': obs['sample'], 'group': np.where(obs.comparison_species.eq('human'), 'Human', 'Mouse'),
                          'segment_class': obs.segment_class, 'pseudospace': pd.Series(coordinate)[obs.index].to_numpy()},
                         index=obs.index),
        var=pd.DataFrame({'human_symbol': [g.upper() for g in shared]}, index=shared))
    return results, data, pd.Series(coordinate), export


@pytest.mark.parametrize('scale', [(1., 0.), (2., -.5)])
def test_hvm_inputs_reproduce_pt_inputs(tmp_path, scale):
    """With pt_inputs' position rule, the export gives pt_inputs' dict exactly, mt genes excluded or not."""
    results, data, coordinate, export = _hvm(tmp_path, scale=scale)
    for exclude in ((), MT):
        old = pt_inputs(results, data, coordinate, set_sizes=(3, 300), exclude_genes=exclude)
        new = hvm_inputs(export, data, ORTHOLOGS, set_sizes=(3, 300), exclude_genes=exclude, rescale='unit')
        assert new['genes'] == old['genes'] and list(new['gene_sets'].items()) == list(old['gene_sets'].items())
        for key in ('position', 'human', 'specimen', 'segment', 'library', 'grid', 'strata'):
            np.testing.assert_array_equal(new[key], old[key])
        assert (new['counts'] != old['counts']).nnz == 0
        np.testing.assert_allclose(new['Y'].toarray(), old['Y'].toarray(), rtol=0, atol=1e-12)
        columns = ['mean_expression', 'detection', 'coverage', 'detection_mouse', 'detection_human']
        np.testing.assert_allclose(new['covariates'][columns], old['covariates'][columns], rtol=0, atol=1e-12)
    # Rare Gene13 fails the detection rule; Gene12 (human structural zero) and mt-Atp8 are not shared orthologs.
    assert 'Gene13' not in new['genes'] and 'Gene12' not in new['genes'] and 'mt-Co1' not in new['genes']
    assert new['excluded_library'] == ['mt-Co1'] and new['excluded_universe'] == ['mt-Co1']
    np.testing.assert_allclose(new['library_with_excluded'] - new['library'],
                               export.X[:, list(export.var_names).index('mt-Co1')].toarray().ravel()[
                                   common_support_mask(export.obs.pseudospace, export.obs['sample'])])


def test_support_rule_is_affine_invariant_and_ignores_outliers():
    rng = np.random.default_rng(0)
    specimen = np.repeat(['a', 'b', 'c'], 300)
    raw = rng.uniform(0, 1, 900)
    raw[[5, 400]] = [6., 9.]                      # far outliers: beyond every 99th percentile
    segment = _segment(np.clip(raw, 0, 1))
    mask, position = support_positions(raw, specimen, segment)
    assert position.min() == 0 and position.max() == 1 and not mask[[5, 400]].any()
    mask2, position2 = support_positions(.16 * raw + .02, specimen, segment)
    np.testing.assert_array_equal(mask, mask2)
    np.testing.assert_allclose(position, position2, rtol=0, atol=1e-12)
    # pt_inputs' rule rescales over every unit: the outliers compress the PT into a fraction of [0, 1].
    _, unit = support_positions(raw, specimen, segment, rescale='unit')
    assert unit.max() < .2
    with pytest.raises(ValueError, match='S1 to S3'):
        support_positions(-raw, specimen, segment)


def _aki_export(n=80, seed=5):
    ad = pytest.importorskip('anndata')
    rng = np.random.default_rng(seed)
    genes = MOUSE + ['Havcr1', 'NotInMap']
    rates = np.r_[rng.gamma(2., 2., 12), .4, .002, 15., 5., 2., 3.]
    rows, obs = [], []
    for sample, group in (('Ctrl1A2', 'Control'), ('Ctrl1A4', 'Control'), ('IR2A2', 'IR'), ('IR2A4', 'IR')):
        x = np.sort(rng.uniform(0, 1, n))
        counts = rng.poisson(rates * rng.gamma(4., .25, (n, 1)))
        counts[:, genes.index('Gene11')] = 0        # never measured in the PT: out of the library
        rows.append(counts)
        obs += [{'id': f'{sample}_unit_{i}', 'sample': sample, 'group': group, 'segment_class': s, 'pseudospace': .3 * p}
                for i, (p, s) in enumerate(zip(x, _segment(x)))]
    return ad.AnnData(X=sparse.csr_matrix(np.vstack(rows).astype(np.int32)), obs=pd.DataFrame(obs).set_index('id'),
                      var=pd.DataFrame(index=genes))


def test_aki_inputs_follow_notebook_69(tmp_path):
    _gene_sets(tmp_path)
    export = _aki_export()
    universe = [g for g in MOUSE if g != 'Gene5']    # the human vs mouse universe: Gene5 is not in it
    d = aki_inputs(export, tmp_path, ORTHOLOGS, universe, exclude_genes=MT, set_sizes=(3, 300))
    keep = common_support_mask(export.obs.pseudospace, export.obs['sample'])
    raw = export.X.toarray()[keep]
    names = list(export.var_names)
    panel = [g for g in MOUSE if g != 'Gene11']     # measured ortholog panel: Havcr1 and NotInMap are not in the map
    full = raw[:, [names.index(g) for g in panel]].sum(axis=1)
    np.testing.assert_allclose(d['library_with_excluded'], full)
    np.testing.assert_allclose(d['library'], full - raw[:, [names.index(g) for g in MT]].sum(axis=1))
    assert d['position'].min() == 0 and d['position'].max() == 1 and d['case'].sum() == (keep & export.obs.group.eq('IR')).sum()
    # Universe: in the human vs mouse universe, measured, detected in >= 2% (equal-specimen mean), not excluded.
    assert d['genes'] == [g for g in universe if g not in MT + ['Gene11', 'Gene13']]
    assert sorted(d['excluded_library']) == sorted(MT) and sorted(d['excluded_universe']) == sorted(MT)
    # Y keeps the full library (mt reads included).
    j = names.index('Gene0')
    np.testing.assert_allclose(d['Y'][:, 0].toarray().ravel(), np.log1p(raw[:, j] * 1e4 / full))
    np.testing.assert_array_equal(d['counts'][:, 0].toarray().ravel(), raw[:, j])
    # Pathways: sizes on the human vs mouse universe minus mt, then restricted to the AKI genes, both with 3-300.
    for lib in LIBRARIES:
        assert sorted(d['gene_sets'][f'{lib}::transport']) == ['Gene0', 'Gene1', 'Gene2', 'Gene3']
        assert d['gene_sets'][f'{lib}::fatty'] == ['Gene12', 'Gene4', 'Gene6']   # Gene5 is outside the universe
        assert f'{lib}::rare' not in d['gene_sets'] and f'{lib}::mito' not in d['gene_sets']
    assert len(d['strata']) == len(d['genes']) and d['Y'].shape == d['counts'].shape == (keep.sum(), len(d['genes']))


def test_contract_guards(tmp_path):
    export = _aki_export()
    wrong = export.copy()
    wrong.obs['group'] = wrong.obs.group.astype(str).replace({'IR': 'Control'})
    with pytest.raises(ValueError, match='cohort'):
        read_pt_units(wrong, 'aki')
    fractional = export.copy()
    fractional.X = fractional.X.astype(float) * .5
    with pytest.raises(ValueError, match='integer'):
        read_pt_units(fractional, 'aki')
    with pytest.raises(ValueError, match='cohort'):
        read_pt_units(export, 'hvm')
    _, data, _, hvm = _hvm(tmp_path)
    hvm.var['human_symbol'] = hvm.var.human_symbol.str.lower()
    with pytest.raises(ValueError, match='human symbols'):
        hvm_inputs(hvm, data, ORTHOLOGS)
