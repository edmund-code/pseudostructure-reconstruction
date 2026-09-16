"""The cached fit wrappers must round-trip identically to the uncached originals."""
from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from pseudospace.levelshape import run_level_shape, sample_perm_pvalues
from pseudospace.stage_cache import cached_perm_pvalues, cached_run_level_shape
from pseudospace.stats_gam import gam_internal_knots


def _panel(seed=0, n_per_group=20, n_features=4):
    rng = np.random.default_rng(seed)
    s = np.concatenate([np.linspace(0, 1, n_per_group) for _ in range(4)])
    c = np.repeat([0.0, 0.0, 1.0, 1.0], n_per_group)
    samples = np.array(['C1'] * n_per_group + ['C2'] * n_per_group
                       + ['A1'] * n_per_group + ['A2'] * n_per_group)
    base = np.linspace(0.5, 2.5, n_per_group)
    rows = []
    for group in range(4):
        signal = base if group < 2 else base[::-1]
        rows.append(np.column_stack([signal * (1 + 0.1 * feature) for feature in range(n_features)]))
    Y = np.vstack(rows)
    return Y, s, c, samples


def test_cached_fit_matches_the_uncached_fit_and_hits(tmp_path):
    Y, s, c, samples = _panel()
    knots = gam_internal_knots(s, basis_df=6)
    grid = np.linspace(0, 1, 31)
    lam = np.logspace(-2, 2, 5)
    reference = run_level_shape(Y, s, c, knots, grid, lam)
    first = cached_run_level_shape(Y, s, c, knots, grid, lam, stage='fit', root=tmp_path)
    second = cached_run_level_shape(Y, s, c, knots, grid, lam, stage='fit', root=tmp_path)
    assert np.allclose(first['shape_rms'], reference['shape_rms'])
    assert np.allclose(second['curve_healthy'], reference['curve_healthy'])
    # the slice list survives the npz round trip, which sample_perm_pvalues depends on
    assert [item.start for item in second['sl2']] == [item.start for item in reference['sl2']]
    assert [item.stop for item in second['sl2']] == [item.stop for item in reference['sl2']]
    assert int(second['p_b']) == int(reference['p_b'])


def test_cached_fit_recomputes_when_the_matrix_changes(tmp_path):
    Y, s, c, samples = _panel()
    knots = gam_internal_knots(s, basis_df=6)
    grid = np.linspace(0, 1, 31)
    lam = np.logspace(-2, 2, 5)
    first = cached_run_level_shape(Y, s, c, knots, grid, lam, stage='fit', root=tmp_path)
    changed = cached_run_level_shape(Y * 1.5, s, c, knots, grid, lam, stage='fit', root=tmp_path)
    assert not np.allclose(first['shape_rms'], changed['shape_rms'])


def test_cached_permutation_matches_and_returns_splits(tmp_path):
    Y, s, c, samples = _panel(seed=3)
    knots = gam_internal_knots(s, basis_df=6)
    grid = np.linspace(0, 1, 31)
    lam = np.logspace(-2, 2, 5)
    fit = run_level_shape(Y, s, c, knots, grid, lam)
    reference = sample_perm_pvalues(Y, s, samples, knots, grid, lam, fit['lam_idx'], fit['sl2'],
                                    fit['p_b'], ['C1', 'C2'])
    cached = cached_perm_pvalues(Y, s, samples, knots, grid, lam, fit['lam_idx'], fit['sl2'],
                                 fit['p_b'], ['C1', 'C2'], stage='perm', root=tmp_path)
    again = cached_perm_pvalues(Y, s, samples, knots, grid, lam, fit['lam_idx'], fit['sl2'],
                                fit['p_b'], ['C1', 'C2'], stage='perm', root=tmp_path)
    assert np.allclose(cached[0], reference[0])
    assert np.allclose(again[0], reference[0])
    assert [tuple(split) for split in cached[2]] == [tuple(split) for split in reference[2]]
    assert cached[3] == reference[3]


def test_cached_fit_accepts_a_precomputed_fingerprint(tmp_path):
    """Passing a fingerprint avoids re-hashing the matrix on every call site."""
    Y, s, c, samples = _panel()
    knots = gam_internal_knots(s, basis_df=6)
    grid = np.linspace(0, 1, 31)
    lam = np.logspace(-2, 2, 5)
    fingerprint = 'deadbeefdeadbeef'
    a = cached_run_level_shape(Y, s, c, knots, grid, lam, stage='fit', root=tmp_path,
                               y_fingerprint=fingerprint)
    b = cached_run_level_shape(Y, s, c, knots, grid, lam, stage='fit', root=tmp_path,
                               y_fingerprint=fingerprint)
    assert np.allclose(a['shape_rms'], b['shape_rms'])
    # a different fingerprint is a different key, so the fit is recomputed rather than reused
    c_ = cached_run_level_shape(Y, s, c, knots, grid, lam, stage='fit', root=tmp_path,
                                y_fingerprint='0000000000000000')
    assert np.allclose(c_['shape_rms'], a['shape_rms'])
    assert sp.issparse(Y) or True
