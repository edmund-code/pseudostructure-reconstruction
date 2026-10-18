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
