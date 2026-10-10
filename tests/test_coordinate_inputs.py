"""pt_inputs on synthetic specimens: excluded genes leave the library and the universe, nothing else moves."""
import json

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from pseudospace.coordinate_inputs import EXPECTED, LIBRARIES, pt_inputs


def _write_inputs(tmp_path, n=40, seed=7):
    ad = pytest.importorskip('anndata')
    rng = np.random.default_rng(seed)
    mouse = [f'Gene{i}' for i in range(12)] + ['mt-Co1', 'mt-Atp8']
    orthologs = pd.DataFrame({'mouse_symbol': mouse, 'human_symbol': [g.upper() for g in mouse],
                              'mapping_status': 'hcop_one_to_one'})
    results, data = tmp_path / 'results', tmp_path / 'data'
    (results / 'minimal_pt_scfates').mkdir(parents=True)
    (data / 'tubule_by_gene').mkdir(parents=True)
    (data / 'mouse_vs_human' / 'pathway_gene_sets').mkdir(parents=True)
    orthologs.to_csv(results / 'minimal_pt_scfates' / 'ortholog_map_used.csv', index=False)
    obs, coordinate = [], {}
    for species, names in EXPECTED.items():
        for name in names:
            # Human panels do not measure MT-ATP8 (as in the real panels): it never enters the library.
            symbols = mouse if species == 'mouse' else [g.upper() for g in mouse if g != 'mt-Atp8']
            counts = rng.poisson(3., (n, len(symbols))).astype(float)
            counts[:, symbols.index('mt-Co1' if species == 'mouse' else 'MT-CO1')] += 20
            ids = [f's{i}' for i in range(n)]
            ad.AnnData(X=sparse.csr_matrix(counts), obs=pd.DataFrame(index=ids),
                       var=pd.DataFrame(index=symbols)).write_h5ad(data / 'tubule_by_gene' / f'{name}_tubule_by_gene_caleb.h5ad')
            position = np.linspace(0, 1, n)
            for i, x in zip(ids, position):
                coordinate[f'{name}_{i}'] = x
                obs.append({'id': f'{name}_{i}', 'sample': name, 'comparison_species': species,
                            'broad_tubule_marker_call': 'PT',
                            'segment_class': 'PT-S1' if x < 1 / 3 else 'PT-S2' if x < 2 / 3 else 'PT-S3'})
    obs = pd.DataFrame(obs).set_index('id')
    ad.AnnData(X=sparse.csr_matrix((len(obs), 1)), obs=obs).write_h5ad(
        results / 'minimal_pt_scfates' / 'cross_species_pt_scfates.h5ad')
    sets = {'mixed': ['GENE0', 'GENE1', 'GENE2', 'GENE3', 'MT-CO1'], 'small': ['GENE4', 'GENE5', 'MT-CO1'],
            'nuclear': ['GENE6', 'GENE7', 'GENE8']}
    for lib in LIBRARIES:
        (data / 'mouse_vs_human' / 'pathway_gene_sets' / f'{lib}.json').write_text(json.dumps(sets))
    return results, data, pd.Series(coordinate)


def test_excluded_genes_leave_library_and_universe_only(tmp_path):
    results, data, coordinate = _write_inputs(tmp_path)
    full = pt_inputs(results, data, coordinate, set_sizes=(3, 300))
    reduced = pt_inputs(results, data, coordinate, set_sizes=(3, 300), exclude_genes=['mt-Co1', 'mt-Atp8'])

    # The default keeps mt-Co1 in the universe and in the library; mt-Atp8 is not measured in human, so it is in neither.
    assert 'mt-Co1' in full['genes'] and 'mt-Atp8' not in full['genes']
    mt = full['genes'].index('mt-Co1')
    np.testing.assert_allclose(full['library'], np.asarray(full['counts'].sum(axis=1)).ravel())
    # Excluding it removes exactly its counts from the library and its column from the universe.
    assert reduced['genes'] == [g for g in full['genes'] if g != 'mt-Co1']
    np.testing.assert_allclose(reduced['library'], full['library'] - full['counts'][:, mt].toarray().ravel())
    keep = [i for i, g in enumerate(full['genes']) if g != 'mt-Co1']
    np.testing.assert_array_equal(reduced['counts'].toarray(), full['counts'][:, keep].toarray())
    # Structures and the log-normalized matrix (full-library normalization) do not move.
    np.testing.assert_array_equal(reduced['position'], full['position'])
    np.testing.assert_array_equal(reduced['specimen'], full['specimen'])
    np.testing.assert_allclose(reduced['Y'].toarray(), full['Y'][:, keep].toarray())
    # Pathways are rebuilt on the reduced universe with the same size rule: 'small' falls below 3 members.
    for prefix in (f'{lib}::' for lib in LIBRARIES):
        assert sorted(full['gene_sets'][prefix + 'mixed']) == ['Gene0', 'Gene1', 'Gene2', 'Gene3', 'mt-Co1']
        assert sorted(reduced['gene_sets'][prefix + 'mixed']) == ['Gene0', 'Gene1', 'Gene2', 'Gene3']
        assert prefix + 'small' in full['gene_sets'] and prefix + 'small' not in reduced['gene_sets']
        assert reduced['gene_sets'][prefix + 'nuclear'] == full['gene_sets'][prefix + 'nuclear']
