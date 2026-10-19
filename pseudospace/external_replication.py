"""Scoring our call lists in an independent experiment (notebook 94; protocol in docs/results/pt-external-iri.md).

* ``dge_group_sums``: streams a dense genes x barcodes text matrix (GEO "dge" files) in row chunks and returns, per
  gene, the summed counts of each barcode group, optionally keeping some barcodes' counts as a sparse matrix. The
  dense files are too large to load whole.
* ``log_cpm``: voom-style log2((y + 0.5) / (N + 1) * 1e6) of summed counts.
* ``segment_markers`` and ``assign_by_centroid``: one-vs-rest marker genes from segment pseudobulks, and the nearest
  segment centroid of each nucleus by Spearman correlation (the assigned-nuclei sensitivity).
* ``replicated``: an item replicates when the external q passes and the external sign equals the item's sign.
* ``matched_null``: random lists drawn from a pool with the list's per-stratum counts kept; the null number replicated.
* ``matched_means``: the same draws, returning the mean of continuous per-item statistics (notebook 99).
* ``h5_group_sums``: per-label sums of rows of an AnnData h5 CSR matrix (``raw/X``), streamed in row chunks (notebook 98's
  KPMP pseudobulks; notebook 109).
"""
from __future__ import annotations

import gzip
from pathlib import Path

import h5py

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata

PRIOR_COUNT = .5


def _open(path):
    path = Path(path)
    return gzip.open(path, 'rt') if path.suffix == '.gz' else open(path)


def dge_group_sums(path, groups, *, keep=None, chunksize=1000):
    """Per gene, summed counts of each barcode group in a tab-separated genes x barcodes count file.

    The first line holds the barcodes; each later line is a gene symbol followed by one whole-number count per barcode.
    ``groups``: barcode -> group label; barcodes not in it are skipped. ``keep``: barcodes whose counts are also
    returned, as a sparse barcodes x genes matrix in ``keep`` order. Returns (genes, sums DataFrame genes x groups in
    first-seen order of ``groups``' labels, kept matrix or None). Raises on a non-integer count or a missing barcode.
    """
    with _open(path) as handle:
        barcodes = handle.readline().rstrip('\n').split('\t')
    if barcodes and barcodes[0] == '':
        barcodes = barcodes[1:]
    position = {b: i for i, b in enumerate(barcodes)}
    if len(position) != len(barcodes):
        raise ValueError(f'{path}: duplicated barcodes')
    missing = [b for b in groups if b not in position]
    if missing:
        raise ValueError(f'{path}: {len(missing)} grouped barcodes are not columns, e.g. {missing[:3]}')
    labels = list(dict.fromkeys(groups.values()))
    label_index = {g: k for k, g in enumerate(labels)}
    columns = np.array([position[b] for b in groups])
    indicator = sparse.csr_matrix((np.ones(len(columns)), (columns, [label_index[groups[b]] for b in groups])),
                                  shape=(len(barcodes), len(labels)))
    keep_columns = None if keep is None else np.array([position[b] for b in keep])
    genes, sums, kept = [], [], []
    # an integer parse fails on any non-integer count, so the dtype is the whole-number check
    dtypes = {0: str, **{i: np.int64 for i in range(1, len(barcodes) + 1)}}
    reader = pd.read_csv(path, sep='\t', header=None, skiprows=1, index_col=0, chunksize=chunksize, dtype=dtypes)
    for chunk in reader:
        values = chunk.to_numpy()
        if values.shape[1] != len(barcodes):
            raise ValueError(f'{path}: a row has {values.shape[1]} counts for {len(barcodes)} barcodes')
        if (values < 0).any():
            raise ValueError(f'{path}: counts must be nonnegative')
        genes.extend(chunk.index.astype(str))
        sums.append(np.asarray(indicator.T @ values.T).T)
        if keep_columns is not None:
            kept.append(sparse.csr_matrix(values[:, keep_columns].T))
    table = pd.DataFrame(np.vstack(sums), index=genes, columns=labels)
    matrix = None if keep_columns is None else sparse.hstack(kept).tocsr()
    return genes, table, matrix


def log_cpm(counts, library=None):
    """log2((y + 0.5) / (N + 1) * 1e6) per column; N = the column sums unless ``library`` is given."""
    y = np.asarray(counts, float)
    n = y.sum(axis=0) if library is None else np.asarray(library, float)
    return np.log2((y + PRIOR_COUNT) / (n + 1.) * 1e6)


def segment_markers(log_cpm_by_segment, cpm_by_segment, *, n=50, min_cpm=10.):
    """For each segment (column), the ``n`` genes with the largest one-vs-rest log2 CPM difference.

    One-vs-rest = the segment minus the mean of the other segments. A gene qualifies for a segment only with CPM >=
    ``min_cpm`` there. Returns {segment: [genes]} in decreasing difference; a gene can serve one segment only (its
    largest difference wins).
    """
    lc = pd.DataFrame(log_cpm_by_segment)
    cpm = pd.DataFrame(cpm_by_segment).reindex_like(lc)
    diff = pd.DataFrame({s: lc[s] - lc.drop(columns=s).mean(axis=1) for s in lc.columns})
    best = diff.idxmax(axis=1)
    out = {}
    for s in lc.columns:
        ok = (cpm[s] >= min_cpm) & (best == s)
        out[s] = diff.loc[ok, s].sort_values(ascending=False).index[:n].tolist()
    return out


def _row_ranks(x):
    r = rankdata(np.asarray(x, float), axis=1)
    r -= r.mean(axis=1, keepdims=True)
    norm = np.sqrt((r ** 2).sum(axis=1, keepdims=True))
    return r / np.where(norm > 0, norm, np.nan)


def assign_by_centroid(profiles, centroids):
    """Each row of ``profiles`` (nuclei x markers) goes to the centroid (rows of ``centroids``) of highest Spearman r.

    Returns (index of the chosen centroid per nucleus, the r matrix). A nucleus with constant markers gets -1.
    """
    a, c = _row_ranks(profiles), _row_ranks(centroids)
    r = a @ c.T
    choice = np.where(np.isfinite(r).all(axis=1), np.nanargmax(np.where(np.isfinite(r), r, -np.inf), axis=1), -1)
    return choice, r


def replicated(external_sign, external_q, internal_sign, *, alpha=.05):
    """True where the external q <= ``alpha`` and the external sign equals the internal sign (nonzero)."""
    es, q, s = (np.asarray(x, float) for x in (external_sign, external_q, internal_sign))
    return (q <= alpha) & (np.sign(es) == np.sign(s)) & (s != 0)


def matched_null(hit, pool_strata, list_strata, *, n=1000, seed=94):
    """The number of hits in ``n`` random lists drawn from the pool with the list's count in each stratum.

    ``hit``: bool per pool item; ``pool_strata``: stratum per pool item; ``list_strata``: stratum per list member.
    Members are drawn without replacement within each stratum. Raises when a stratum holds fewer pool items than list
    members.
    """
    hit = np.asarray(hit, bool)
    pool_strata, list_strata = np.asarray(pool_strata), np.asarray(list_strata)
    if hit.shape != pool_strata.shape:
        raise ValueError('Need one stratum per pool item.')
    rng = np.random.default_rng(seed)
    out = np.zeros(n, int)
    for stratum, count in zip(*np.unique(list_strata, return_counts=True)):
        members = np.flatnonzero(pool_strata == stratum)
        if len(members) < count:
            raise ValueError(f'Stratum {stratum!r} holds {len(members)} pool items for {count} list members.')
        keys = rng.random((n, len(members)))                     # a random subset per draw: the smallest keys
        chosen = np.argpartition(keys, count - 1, axis=1)[:, :count]
        out += hit[members][chosen].sum(axis=1)
    return out


def matched_means(values, pool_strata, list_strata, *, n=1000, seed=99):
    """Mean of ``values`` (pool items x statistics) over ``n`` random lists with the list's per-stratum counts."""
    values = np.asarray(values, float)
    if values.ndim == 1:
        values = values[:, None]
    pool_strata, list_strata = np.asarray(pool_strata), np.asarray(list_strata)
    if values.shape[0] != len(pool_strata):
        raise ValueError('Need one stratum per pool item.')
    rng = np.random.default_rng(seed)
    total = np.zeros((n, values.shape[1]))
    for stratum, count in zip(*np.unique(list_strata, return_counts=True)):
        members = np.flatnonzero(pool_strata == stratum)
        if len(members) < count:
            raise ValueError(f'Stratum {stratum!r} holds {len(members)} pool items for {count} list members.')
        chosen = np.argpartition(rng.random((n, len(members))), count - 1, axis=1)[:, :count]
        total += values[members][chosen].sum(axis=1)
    return total / len(list_strata)


def null_summary(k, null):
    """Precision-free summary of an observed hit count against its null: mean, 2.5/97.5% band and empirical p."""
    null = np.asarray(null, float)
    return {'null mean': float(null.mean()), 'null low': float(np.quantile(null, .025)),
            'null high': float(np.quantile(null, .975)), 'p': float((1 + (null >= k).sum()) / (1 + len(null)))}


def h5_group_sums(path, rows, labels, *, key='raw/X', chunk=20000):
    """genes x labels sums of the CSR rows ``rows`` (with one label each) of the matrix at ``key`` in an h5ad file."""
    rows = np.asarray(rows)
    lab = pd.Index(sorted(set(labels)))
    code = pd.Series(lab.get_indexer(labels), index=rows)
    with h5py.File(path, 'r') as f:
        X = f[key]
        indptr = X['indptr'][:]
        n_rows, n_genes = X.attrs['shape']
        out = np.zeros((n_genes, len(lab)))
        for start in range(0, n_rows, chunk):
            stop = min(start + chunk, n_rows)
            want = code[(code.index >= start) & (code.index < stop)]
            if want.empty:
                continue
            lo, hi = indptr[start], indptr[stop]
            block = sparse.csr_matrix((X['data'][lo:hi], X['indices'][lo:hi], indptr[start:stop + 1] - lo),
                                      shape=(stop - start, n_genes))
            sub = block[want.index.to_numpy() - start]
            ind = sparse.csr_matrix((np.ones(len(want)), (want.to_numpy(), np.arange(len(want)))), shape=(len(lab), len(want)))
            out += (ind @ sub).T.toarray()
    return pd.DataFrame(out, columns=lab)
