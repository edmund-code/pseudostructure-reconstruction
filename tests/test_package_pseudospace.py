"""Notebook 70's helpers on synthetic data: shared orthologs, the ortholog-space bin matrix, polygon matching."""
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from pseudospace.package_pseudospace import (adjacent_aucs, h5_gene_names, match_polygons, pt_boundary,
                                             shared_ortholog_genes, to_shared_genes, uns_safe, validate_ortholog_map,
                                             write_shared_bin_dir)

MAP = pd.DataFrame({'human_symbol': ['GENEA', 'GENEB', 'GENEC', 'MT-CO1', 'GENED'],
                    'mouse_symbol': ['Genea', 'Geneb', 'Genec', 'mt-Co1', 'Gened'],
                    'mapping_status': 'hcop_one_to_one'})


def test_map_must_be_one_to_one():
    validate_ortholog_map(MAP)
    with pytest.raises(ValueError):
        validate_ortholog_map(pd.concat([MAP, MAP.iloc[[0]].assign(mouse_symbol='Other')]))
    with pytest.raises(ValueError):
        validate_ortholog_map(MAP.assign(human_symbol=['GENEA', 'genea', 'GENEC', 'MT-CO1', 'GENED']))


def test_shared_genes_need_every_sample_and_follow_mouse_order():
    mouse = [['mt-Co1', 'Genec', 'Genea', 'Xmouse', 'Gened'], ['Genea', 'Genec', 'mt-Co1', 'Gened']]
    human = [['GENEA', 'Genec', 'MT-CO1', 'XHUMAN', 'GENEB'], ['GENEA', 'Genec', 'MT-CO1', 'GENED']]
    table = shared_ortholog_genes(MAP, mouse, human)
    # Geneb is not a mouse gene; Gened is missing in the first human sample; human matching ignores case
    assert table.mouse_symbol.tolist() == ['mt-Co1', 'Genec', 'Genea']
    assert table.human_source.tolist() == ['MT-CO1', 'Genec', 'GENEA']
    with pytest.raises(ValueError):    # one human gene spelled two ways across samples
        shared_ortholog_genes(MAP, mouse, [['GENEC', 'GENEA'], ['Genec', 'GENEA']])


def test_to_shared_genes_renames_without_touching_counts():
    ad = pytest.importorskip('anndata')
    table = shared_ortholog_genes(MAP, [['Genea', 'Genec']], [['GENEC', 'XHUMAN', 'GENEA']])
    x = sparse.csr_matrix(np.arange(6, dtype=np.float32).reshape(2, 3))
    human = ad.AnnData(X=x, obs=pd.DataFrame(index=['unit_0', 'unit_1']), var=pd.DataFrame(index=['GENEC', 'XHUMAN', 'GENEA']))
    out = to_shared_genes(human, table, 'human')
    assert out.var_names.tolist() == ['Genea', 'Genec']
    assert out.var.source_symbol.tolist() == ['GENEA', 'GENEC']
    np.testing.assert_array_equal(out.X.toarray(), x.toarray()[:, [2, 0]])
    with pytest.raises(ValueError):
        to_shared_genes(human, table, 'mouse')


def _write_10x(path, names, counts):
    h5py = pytest.importorskip('h5py')
    csc = sparse.csc_matrix(counts.T.astype(np.int32))         # genes x barcodes
    path.mkdir(parents=True)
    (path / 'spatial').mkdir()
    (path / 'spatial' / 'scalefactors_json.json').write_text('{"spot_diameter_fullres": 4.5, "microns_per_pixel": 0.44}')
    with h5py.File(path / 'filtered_feature_bc_matrix.h5', 'w') as f:
        m = f.create_group('matrix')
        m['barcodes'] = np.array([f's_002um_{i:05d}-1'.encode() for i in range(counts.shape[0])])
        m['data'], m['indices'], m['indptr'] = csc.data, csc.indices.astype(np.int64), csc.indptr.astype(np.int64)
        m['shape'] = np.array(csc.shape, dtype=np.int32)
        g = m.create_group('features')
        g['_all_tag_keys'] = np.array([b'genome'])
        g['name'] = np.array([n.encode() for n in names])
        g['id'] = np.array([f'ID{i}'.encode() for i in range(len(names))])
        g['feature_type'] = np.array([b'Gene Expression'] * len(names))
        g['genome'] = np.array([b'GRCh38'] * len(names))


def test_shared_bin_dir_keeps_bins_and_counts(tmp_path):
    sc = pytest.importorskip('scanpy')
    rng = np.random.default_rng(0)
    counts = rng.poisson(1., (50, 4))
    names = ['GENEA', 'DUP', 'DUP', 'GENEC']               # a duplicated name becomes DUP-1, as in stage 01
    _write_10x(tmp_path / 'src', names, counts)
    assert h5_gene_names(tmp_path / 'src' / 'filtered_feature_bc_matrix.h5').tolist() == ['GENEA', 'DUP', 'DUP-1', 'GENEC']
    write_shared_bin_dir(tmp_path / 'src', tmp_path / 'dst', ['GENEC', 'DUP-1'], ['Genec', 'Dup'])
    out = sc.read_10x_h5(str(tmp_path / 'dst' / 'filtered_feature_bc_matrix.h5'))
    assert out.var_names.tolist() == ['Genec', 'Dup']
    assert out.var.gene_ids.tolist() == ['ID3', 'ID2']
    assert out.obs_names.tolist() == [f's_002um_{i:05d}-1' for i in range(50)]
    np.testing.assert_array_equal(out.X.toarray(), counts[:, [3, 2]])
    assert (tmp_path / 'dst' / 'spatial' / 'scalefactors_json.json').read_text().startswith('{"spot')
    with pytest.raises(ValueError):
        write_shared_bin_dir(tmp_path / 'src', tmp_path / 'bad', ['NOPE'], ['Nope'])


def test_match_polygons_pairs_by_iou():
    geometry = pytest.importorskip('shapely.geometry')
    box = geometry.box
    old = [box(0, 0, 10, 10), box(20, 0, 30, 10), box(40, 0, 50, 10), box(60, 0, 70, 10)]
    new = [box(21, 0, 31, 10), box(0, 0, 10, 10), box(45, 0, 55, 10), box(100, 0, 110, 10)]
    pairs = match_polygons(old, new, min_iou=.5)
    assert pairs[['old', 'new']].values.tolist() == [[0, 1], [1, 0]]     # IoU 1 and 0.82; box 2 has IoU 1/3
    np.testing.assert_allclose(pairs.iou, [1., 90 / 110])
    np.testing.assert_allclose(pairs.centroid_distance, [0., 1.])
    assert match_polygons([], new).empty


def test_adjacent_aucs_and_pt_boundary():
    position = np.array([0., .1, .2, .3, .4, .5, .6, .7])
    segment = np.array(['PT-S1', 'PT-S1', 'PT-S2', 'PT-S2', 'PT-S3', 'PT-S3', 'mTAL', 'PT-S3'])
    aucs = adjacent_aucs(position, segment)
    assert aucs == {'PT-S1->PT-S2': 1.0, 'PT-S2->PT-S3': 1.0}
    assert np.isnan(adjacent_aucs(position, segment, order=('PT-S1', 'DTL'))['PT-S1->DTL'])
    cut = pt_boundary(position, segment != 'mTAL')
    assert cut['cut'] == .5 and cut['pt_beyond'] == 1 and cut['other_before'] == 0


def test_uns_safe_writes(tmp_path):
    ad = pytest.importorskip('anndata')
    value = uns_safe({'a': None, 'b': (1, 2), 'c': [['x', 'y'], 'z'], 'd': {'e': [0.5, 1]}, 'f': True})
    assert value == {'a': 'null', 'b': [1, 2], 'c': '[["x", "y"], "z"]', 'd': {'e': [0.5, 1]}, 'f': True}
    adata = ad.AnnData(X=np.zeros((1, 1)))
    adata.uns['provenance'] = value
    adata.write_h5ad(tmp_path / 'x.h5ad')
    assert ad.read_h5ad(tmp_path / 'x.h5ad').uns['provenance']['a'] == 'null'
