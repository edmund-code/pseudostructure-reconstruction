"""Stage cache: keying, hits/misses, code-change invalidation and purge."""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from pseudospace.stage_cache import (
    cache_status,
    cached_anndata,
    cached_frame,
    cached_payload,
    code_digest,
    digest,
    purge_stage_cache,
    stage_cache_enabled,
    stage_key,
)


def _adata(value=1.0, n=8):
    matrix = sp.csr_matrix(np.full((n, 3), value))
    obj = ad.AnnData(X=matrix)
    obj.var_names = ['g1', 'g2', 'g3']
    obj.layers['lognorm'] = matrix.copy()
    return obj


def test_digest_is_stable_and_sensitive():
    assert digest({'a': 1, 'b': [1.0, 2.0]}) == digest({'b': [1.0, 2.0], 'a': 1})
    assert digest(np.arange(4)) != digest(np.arange(4) + 1)
    assert digest(_adata(1.0).layers['lognorm']) != digest(_adata(2.0).layers['lognorm'])
    frame = pd.DataFrame({'x': [1, 2]})
    assert digest(frame) == digest(frame.copy())
    assert digest(frame) != digest(frame.assign(x=[1, 3]))


def test_code_digest_changes_with_the_implementation():
    def version_one(x):
        return x + 1

    def version_two(x):
        return x + 2

    assert code_digest(version_one) != code_digest(version_two)
    assert code_digest(version_one) == code_digest(version_one)


def test_stage_key_separates_params_inputs_and_code():
    base = dict(params={'grid': 101}, inputs={'y': 'abc'}, code='c1')
    assert stage_key('fit', **base) == stage_key('fit', **base)
    assert stage_key('fit', **base) != stage_key('fit', params={'grid': 51}, inputs={'y': 'abc'}, code='c1')
    assert stage_key('fit', **base) != stage_key('fit', params={'grid': 101}, inputs={'y': 'abd'}, code='c1')
    assert stage_key('fit', **base) != stage_key('fit', params={'grid': 101}, inputs={'y': 'abc'}, code='c2')
    assert stage_key('fit', **base) != stage_key('other', **base)


def test_anndata_cache_returns_the_stored_object_and_recomputes_on_key_change(tmp_path):
    calls = {'n': 0}

    def compute():
        calls['n'] += 1
        return _adata(0.5 + calls['n'])

    first = cached_anndata('fit', compute, root=tmp_path, params={'grid': 101},
                           inputs={'y': 'abc'}, code='c1')
    assert calls['n'] == 1
    second = cached_anndata('fit', compute, root=tmp_path, params={'grid': 101},
                            inputs={'y': 'abc'}, code='c1')
    assert calls['n'] == 1, 'an identical key must not recompute'
    assert np.allclose(first.X.toarray(), second.X.toarray())
    # a changed parameter, input or code must recompute
    cached_anndata('fit', compute, root=tmp_path, params={'grid': 51}, inputs={'y': 'abc'}, code='c1')
    cached_anndata('fit', compute, root=tmp_path, params={'grid': 101}, inputs={'y': 'abd'}, code='c1')
    cached_anndata('fit', compute, root=tmp_path, params={'grid': 101}, inputs={'y': 'abc'}, code='c2')
    assert calls['n'] == 4


def test_payload_and_frame_caches_round_trip(tmp_path):
    calls = {'n': 0}

    def compute_payload():
        calls['n'] += 1
        return {'curves': np.arange(12).reshape(3, 4), 'lam': np.array([0.1])}

    payload = cached_payload('payload', compute_payload, root=tmp_path, code='c1')
    again = cached_payload('payload', compute_payload, root=tmp_path, code='c1')
    assert calls['n'] == 1
    assert np.array_equal(payload['curves'], again['curves'])

    def compute_frame():
        calls['n'] += 1
        return pd.DataFrame({'a': [1, 2], 'b': ['x', 'y']})

    frame = cached_frame('frame', compute_frame, root=tmp_path, code='c1')
    cached_frame('frame', compute_frame, root=tmp_path, code='c1')
    assert calls['n'] == 2
    assert frame.shape == (2, 2)


def test_cache_can_be_disabled_and_env_var_is_honoured(tmp_path, monkeypatch):
    calls = {'n': 0}

    def compute():
        calls['n'] += 1
        return _adata(1.0)

    cached_anndata('off', compute, root=tmp_path, code='c1', enabled=False)
    cached_anndata('off', compute, root=tmp_path, code='c1', enabled=False)
    assert calls['n'] == 2, 'enabled=False must bypass the cache'
    monkeypatch.setenv('PSEUDOSPACE_STAGE_CACHE', '0')
    assert stage_cache_enabled() is False
    monkeypatch.setenv('PSEUDOSPACE_STAGE_CACHE', '1')
    assert stage_cache_enabled() is True
    assert stage_cache_enabled(explicit=False) is False


def test_status_and_purge(tmp_path):
    cached_anndata('fit', lambda: _adata(1.0), root=tmp_path, code='c1', verbose=False)
    cached_frame('frame', lambda: pd.DataFrame({'a': [1]}), root=tmp_path, code='c1', verbose=False)
    status = cache_status(tmp_path)
    assert set(status['stage']) == {'fit', 'frame'}
    assert status['size_mb'].notna().all()
    assert status.set_index('stage').loc['fit', 'size_mb'] > 0
    assert purge_stage_cache(tmp_path, stage='fit') == 2      # payload + sidecar
    assert set(cache_status(tmp_path)['stage']) == {'frame'}
    assert purge_stage_cache(tmp_path) == 2
    assert cache_status(tmp_path).empty
