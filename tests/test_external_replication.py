import gzip

import numpy as np
import pandas as pd
import pytest

from pseudospace import external_replication as er


def _write_dge(path, genes, barcodes, counts):
    with gzip.open(path, 'wt') as handle:
        handle.write('\t'.join(barcodes) + '\n')
        for g, row in zip(genes, counts):
            handle.write(g + '\t' + '\t'.join(str(v) for v in row) + '\n')


def test_dge_group_sums_streams_chunks_and_keeps_barcodes(tmp_path):
    rng = np.random.default_rng(0)
    genes = [f'G{i}' for i in range(7)]
    barcodes = [f'b{j}' for j in range(5)]
    counts = rng.integers(0, 9, (7, 5))
    path = tmp_path / 'm.dge.txt.gz'
    _write_dge(path, genes, barcodes, counts)
    groups = {'b0': 'x', 'b3': 'y', 'b1': 'x'}                 # b2 and b4 are skipped
    got_genes, sums, kept = er.dge_group_sums(path, groups, keep=['b4', 'b0'], chunksize=3)
    assert got_genes == genes and list(sums.columns) == ['x', 'y']
    assert np.array_equal(sums['x'].to_numpy(), counts[:, 0] + counts[:, 1])
    assert np.array_equal(sums['y'].to_numpy(), counts[:, 3])
    assert np.array_equal(kept.toarray(), counts[:, [4, 0]].T)
    with pytest.raises(ValueError):
        er.dge_group_sums(path, {'nope': 'x'})
    bad = tmp_path / 'bad.dge.txt.gz'
    _write_dge(bad, ['G0'], ['b0', 'b1'], [['1', '0.5']])
    with pytest.raises(ValueError):
        er.dge_group_sums(bad, {'b0': 'x'})


def test_log_cpm_matches_voom():
    y = np.array([[0, 10], [100, 30]])
    lc = er.log_cpm(y)
    assert lc[0, 0] == pytest.approx(np.log2(.5 / 101 * 1e6))
    assert lc[1, 1] == pytest.approx(np.log2(30.5 / 41 * 1e6))


def test_segment_markers_and_centroid_assignment():
    genes = [f'g{i}' for i in range(9)]
    lc = pd.DataFrame({'S1': [9, 9, 9, 1, 1, 1, 1, 1, 1], 'S2': [1, 1, 1, 9, 8, 9, 1, 1, 1],
                       'S3': [1, 1, 1, 1, 1, 1, 9, 9, 7]}, index=genes, dtype=float)
    cpm = 2. ** lc
    markers = er.segment_markers(lc, cpm, n=2, min_cpm=10)
    assert markers == {'S1': ['g0', 'g1'], 'S2': ['g3', 'g5'], 'S3': ['g6', 'g7']}
    assert er.segment_markers(lc, cpm, n=2, min_cpm=1e6) == {'S1': [], 'S2': [], 'S3': []}
    centroids = np.array([[5, 1, 1, 0], [1, 5, 1, 0], [1, 1, 5, 0]], float)
    nuclei = np.array([[3, 0, 1, 0], [0, 1, 4, 1], [2, 2, 2, 2]], float)
    choice, r = er.assign_by_centroid(nuclei, centroids)
    assert choice.tolist() == [0, 2, -1] and r.shape == (3, 3)


def test_replicated_and_matched_null():
    rep = er.replicated([1., -2., 3., 1.], [.01, .01, .2, .01], [1, 1, 1, 0])
    assert rep.tolist() == [True, False, False, False]
    hit = np.array([True] * 10 + [False] * 90)
    strata = np.array(['a'] * 20 + ['b'] * 80)                 # all hits in stratum a, half of it
    null = er.matched_null(hit, strata, ['a'] * 4, n=2000, seed=1)
    assert null.mean() == pytest.approx(2., abs=.1) and null.max() <= 4
    assert er.matched_null(hit, strata, ['b'] * 5, n=50).sum() == 0
    with pytest.raises(ValueError):
        er.matched_null(hit, strata, ['a'] * 21)
    s = er.null_summary(4, null)
    assert s['p'] == pytest.approx((1 + (null >= 4).sum()) / 2001) and s['null low'] <= s['null mean'] <= s['null high']


def test_matched_means_keep_strata():
    values = np.r_[np.zeros(50), np.ones(50) * 10.]
    strata = np.r_[np.zeros(50, int), np.ones(50, int)]
    m = er.matched_means(values, strata, [0, 0, 1, 1], n=500, seed=1)
    assert m.shape == (500, 1) and np.allclose(m, 5.)            # two from each stratum: always 5
    rng = np.random.default_rng(2)
    v = rng.normal(size=(200, 2))
    m = er.matched_means(v, np.zeros(200, int), np.zeros(20, int), n=2000, seed=3)
    assert np.allclose(m.mean(axis=0), v.mean(axis=0), atol=.02) and m.std(axis=0)[0] < v.std(axis=0)[0] / 3


def test_h5_group_sums(tmp_path):
    h5py = pytest.importorskip('h5py')
    from scipy import sparse
    X = sparse.csr_matrix(np.array([[1, 0, 2], [0, 3, 0], [4, 0, 0], [0, 0, 5]]))
    path = tmp_path / 'm.h5'
    with h5py.File(path, 'w') as f:
        g = f.create_group('raw/X')
        g.attrs['shape'] = X.shape
        for k in ('data', 'indices', 'indptr'):
            g.create_dataset(k, data=getattr(X, k))
    out = er.h5_group_sums(path, [0, 2, 3], ['a', 'a', 'b'], chunk=2)
    assert out['a'].tolist() == [5, 0, 2] and out['b'].tolist() == [0, 0, 5]
