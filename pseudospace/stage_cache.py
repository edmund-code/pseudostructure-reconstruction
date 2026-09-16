"""Content-addressed stage cache so "Run All" only recomputes what actually changed.

A stage that is expensive but pure (given its inputs and parameters) can be keyed by

* ``stage`` - a name the notebook chooses,
* ``params`` - every parameter that changes the result (grids, lambda grids, thresholds),
* ``inputs`` - digests of the data it reads (an expression matrix, a saved object, a file),
* ``code`` - the source of the function/module that implements it.

The cache stores the payload under ``<cache root>/<stage>__<key>.<ext>`` with a JSON sidecar
recording when it was written and what the key was. On the next run the same key returns the stored
payload instead of recomputing; any change to the parameters, the inputs or the implementation
produces a different key and the stage is recomputed.

Design rules, because this repository is reproducibility-first:

* A hit is always printed, never silent, so a run log still shows what was actually computed.
* Only pure compute is cached. Validation cells (Harmony version guard, input existence checks,
  fingerprint guards) must run on every execution and are deliberately left outside the cache.
* Caching is on by default and can be disabled with ``PSEUDOSPACE_STAGE_CACHE=0`` or by passing
  ``enabled=False``.
* :data:`CACHE_FORMAT_VERSION` is bumped whenever the stored payload layout changes, which
  invalidates every existing entry.
"""
from __future__ import annotations

import hashlib
import inspect
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

CACHE_FORMAT_VERSION = 1
CACHE_DIRNAME = 'stage_cache'
_SIDECAR_SUFFIX = '.key.json'

__all__ = [
    'CACHE_FORMAT_VERSION',
    'stage_is_fresh',
    'stage_mark_fresh',
    'stage_cache_enabled',
    'stage_cache_root',
    'digest',
    'code_digest',
    'stage_key',
    'cached_anndata',
    'cached_payload',
    'cached_frame',
    'cache_status',
    'purge_stage_cache',
]


def stage_cache_enabled(explicit: bool | None = None) -> bool:
    """Caching is on unless ``PSEUDOSPACE_STAGE_CACHE`` is falsy or an explicit value is given."""
    if explicit is not None:
        return bool(explicit)
    return os.environ.get('PSEUDOSPACE_STAGE_CACHE', '1').strip().lower() not in ('0', 'false', 'no', '')


def stage_cache_root(results_dir, override=None) -> Path:
    root = Path(override) if override is not None else Path(results_dir) / CACHE_DIRNAME
    root.mkdir(parents=True, exist_ok=True)
    return root


def digest(obj) -> str:
    """Stable digest of anything that can change a result.

    Arrays and sparse matrices are hashed over their data and shape, DataFrames over their values and
    axes, files over their size plus their bytes, mappings over their sorted items, everything else
    over its ``repr``. Returns a 16-hex-character string.
    """
    hasher = hashlib.sha1()
    hasher.update(f'v{CACHE_FORMAT_VERSION}|'.encode())

    def _feed(value):
        if value is None or isinstance(value, (bool, int, float, str)):
            hasher.update(repr(value).encode())
        elif isinstance(value, Path):
            _feed_file(value)
        elif sp.issparse(value):
            hasher.update(f'sparse{value.shape}{value.dtype}|'.encode())
            hasher.update(np.ascontiguousarray(value.tocsr().data).view(np.uint8))
            hasher.update(np.ascontiguousarray(value.tocsr().indices).view(np.uint8))
            hasher.update(np.ascontiguousarray(value.tocsr().indptr).view(np.uint8))
        elif isinstance(value, np.ndarray):
            hasher.update(f'array{value.shape}{value.dtype}|'.encode())
            hasher.update(np.ascontiguousarray(value).view(np.uint8))
        elif isinstance(value, pd.DataFrame):
            _feed(value.to_numpy().astype(str) if value.empty is False else value.shape)
            _feed(list(map(str, value.columns)))
            _feed(list(map(str, value.index)))
        elif isinstance(value, pd.Series):
            _feed(value.to_numpy().astype(str))
            _feed(list(map(str, value.index)))
        elif isinstance(value, dict):
            for key in sorted(map(str, value)):
                hasher.update(str(key).encode())
                _feed(value[key] if key in value else value[str(key)])
        elif isinstance(value, (list, tuple, set)):
            for item in value:
                _feed(item)
        elif hasattr(value, 'shape') and hasattr(value, 'dtype'):        # masked array and friends
            _feed(np.asarray(value))
        else:
            hasher.update(repr(value).encode())

    def _feed_file(path: Path):
        path = Path(path)
        if not path.exists():
            hasher.update(f'missing:{path.name}'.encode())
            return
        stat = path.stat()
        hasher.update(f'file:{path.name}:{stat.st_size}:{int(stat.st_mtime)}'.encode())
        # Hash the content of small files fully; for large artifacts the size+mtime pair is the
        # contract, and the expensive fingerprint of in-memory objects is taken separately.
        if stat.st_size <= 8 * 1024 * 1024:
            hasher.update(path.read_bytes())

    _feed(obj)
    return hasher.hexdigest()[:16]


def fingerprint_anndata(adata, layers=('lognorm',), use_data=True) -> str:
    """Digest of the parts of an AnnData that a downstream fit depends on.

    ``use_data=False`` hashes only the shape and the axes, which is cheap but blind to a changed
    matrix; ``use_data=True`` adds a hash of each requested layer (and of ``X`` when the layer is
    absent). Use the latter for anything that feeds a model.
    """
    parts = [f'anndata{adata.shape}', list(map(str, adata.obs_names)), list(map(str, adata.var_names))]
    if use_data:
        for name in layers:
            matrix = adata.layers[name] if name in adata.layers else adata.X
            parts.append((name, matrix))
    return digest(parts)


def code_digest(*objects) -> str:
    """Digest of the implementation behind a stage.

    Accepts functions, classes or modules; the source is read with :mod:`inspect`. Anything whose
    source cannot be retrieved (an object defined interactively without a linecache entry) falls back
    to its qualified name, which still changes when the name changes.
    """
    chunks = []
    for obj in objects:
        if obj is None:
            continue
        try:
            chunks.append(inspect.getsource(obj))
        except (OSError, TypeError):
            chunks.append(getattr(obj, '__qualname__', repr(obj)))
    if not chunks:
        chunks.append('no-code-fingerprint')
    return digest(chunks)


def stage_key(stage: str, *, params=None, inputs=None, code=None) -> str:
    return digest({'stage': stage, 'params': params, 'inputs': inputs, 'code': code})


def _ensure_dir(root) -> Path:
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    return root


def _sidecar(path: Path) -> Path:
    return path.with_suffix(path.suffix + _SIDECAR_SUFFIX)


def _record(path: Path, *, stage, key, params, inputs, code):
    sidecar = _sidecar(path)
    sidecar.write_text(json.dumps({
        'stage': stage,
        'key': key,
        'cache_format_version': CACHE_FORMAT_VERSION,
        'written_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'params': repr(params),
        'inputs': repr(inputs),
        'code_digest': code[:16] if isinstance(code, str) else repr(code),
    }, indent=2) + '\n')


def _hit(path: Path, key: str) -> bool:
    sidecar = _sidecar(path)
    if not path.exists() or not sidecar.exists():
        return False
    try:
        record = json.loads(sidecar.read_text())
    except json.JSONDecodeError:
        return False
    return record.get('key') == key and record.get('cache_format_version') == CACHE_FORMAT_VERSION


def _load_anndata(path: Path):
    import anndata as ad

    return ad.read_h5ad(path)


def cached_anndata(stage, compute, *, root, params=None, inputs=None, code=None, enabled=None,
                   verbose=True):
    """Return ``compute()`` or the cached AnnData for the same stage/params/inputs/code.

    ``compute`` is a zero-argument callable returning the AnnData to cache.
    """
    if not stage_cache_enabled(enabled):
        return compute()
    root = _ensure_dir(root)
    key = stage_key(stage, params=params, inputs=inputs, code=code_digest(code))
    path = root / f'{stage}__{key}.h5ad'
    if _hit(path, key):
        if verbose:
            written = json.loads(_sidecar(path).read_text()).get('written_utc', '?')
            print(f'[stage cache] hit {stage} (key {key}, written {written})')
        return _load_anndata(path)
    if verbose:
        print(f'[stage cache] computing {stage} (key {key})')
    result = compute()
    result.write(path)
    _record(path, stage=stage, key=key, params=params, inputs=inputs, code=code_digest(code))
    return result


def cached_payload(stage, compute, *, root, params=None, inputs=None, code=None, enabled=None,
                   verbose=True, filename_suffix='.npz'):
    """Cache a dict of arrays / scalars, stored with :func:`numpy.savez_compressed`."""
    if not stage_cache_enabled(enabled):
        return compute()
    root = _ensure_dir(root)
    key = stage_key(stage, params=params, inputs=inputs, code=code_digest(code))
    path = root / f'{stage}__{key}{filename_suffix}'
    if _hit(path, key):
        if verbose:
            print(f'[stage cache] hit {stage} (key {key})')
        with np.load(path, allow_pickle=False) as stored:
            return {name: stored[name] for name in stored.files}
    if verbose:
        print(f'[stage cache] computing {stage} (key {key})')
    payload = compute()
    np.savez_compressed(path, **{str(name): np.asarray(value) for name, value in payload.items()})
    _record(path, stage=stage, key=key, params=params, inputs=inputs, code=code_digest(code))
    return payload


def cached_frame(stage, compute, *, root, params=None, inputs=None, code=None, enabled=None,
                 verbose=True):
    """Cache a DataFrame as CSV (read back with ``index_col=False``)."""
    if not stage_cache_enabled(enabled):
        return compute()
    root = _ensure_dir(root)
    key = stage_key(stage, params=params, inputs=inputs, code=code_digest(code))
    path = root / f'{stage}__{key}.csv'
    if _hit(path, key):
        if verbose:
            print(f'[stage cache] hit {stage} (key {key})')
        return pd.read_csv(path)
    if verbose:
        print(f'[stage cache] computing {stage} (key {key})')
    frame = compute()
    frame.to_csv(path, index=False)
    _record(path, stage=stage, key=key, params=params, inputs=inputs, code=code_digest(code))
    return frame


def stage_is_fresh(path, key) -> bool:
    """True when the artifact at ``path`` was produced with this stage key.

    Artifact-freshness caching: the notebook already writes e.g. the pass-1 object, so a matching
    sidecar lets the next run load it instead of recomputing. Any change to the key (parameters,
    inputs, code) makes the artifact stale and the stage recomputes and re-marks it.
    """
    return _hit(Path(path), key)


def stage_mark_fresh(path, key, *, stage=None, params=None, inputs=None, code=None) -> None:
    """Record the key that produced ``path`` (writes ``<path>.key.json``)."""
    _record(Path(path), stage=stage or Path(path).name, key=key, params=params, inputs=inputs,
            code=code or '')


def cached_run_level_shape(Y, s, c, knots, grid, lambda_grid, *, stage, root, y_fingerprint=None,
                           enabled=None, verbose=True):
    """``levelshape.run_level_shape`` with the nested fit stored in the cache.

    The payload is the fit dictionary; ``sl2`` (a list of slices) is stored as a (start, stop) array
    and rebuilt, so callers such as ``sample_perm_pvalues`` keep working unchanged.
    """
    from . import stats_gam
    from .levelshape import run_level_shape

    params = {'grid': np.asarray(grid), 'lambda_grid': np.asarray(lambda_grid),
              'knots': np.asarray(knots)}
    inputs = {'s': np.asarray(s), 'c': np.asarray(c),
              'y': y_fingerprint if y_fingerprint is not None else digest(Y)}
    key = stage_key(stage, params=params, inputs=inputs, code=code_digest(run_level_shape, stats_gam))
    root = _ensure_dir(root) if stage_cache_enabled(enabled) else None
    path = root / f'{stage}__{key}.npz' if root is not None else None
    if path is not None and _hit(path, key):
        if verbose:
            print(f'[stage cache] hit {stage} (key {key})')
        with np.load(path, allow_pickle=False) as stored:
            payload = {name: stored[name] for name in stored.files}
        payload['sl2'] = [slice(int(start), int(stop)) for start, stop in payload['sl2_bounds']]
        return payload
    if verbose:
        print(f'[stage cache] computing {stage} (key {key})')
    fit = run_level_shape(Y, s, c, knots, grid, lambda_grid)
    if path is not None:
        payload = {name: np.asarray(value) for name, value in fit.items() if name != 'sl2'}
        payload['sl2_bounds'] = np.asarray([[item.start, item.stop] for item in fit['sl2']])
        np.savez_compressed(path, **payload)
        _record(path, stage=stage, key=key, params=params, inputs={k: repr(v)[:80] for k, v in inputs.items()},
                code=code_digest(run_level_shape, stats_gam))
    return fit


def cached_perm_pvalues(Y, s, samples, knots, grid, lambda_grid, lam_idx, sl2, p_b,
                        control_samples, *, stage, root, y_fingerprint=None, enabled=None,
                        verbose=True):
    """``levelshape.sample_perm_pvalues`` with the exact permutation block cached."""
    from . import stats_gam as _stats_gam
    from .levelshape import sample_perm_pvalues

    params = {'grid': np.asarray(grid), 'lambda_grid': np.asarray(lambda_grid),
              'knots': np.asarray(knots), 'lam_idx': np.asarray(lam_idx), 'p_b': int(p_b)}
    inputs = {'s': np.asarray(s), 'samples': np.asarray(samples).astype(str),
              'y': y_fingerprint if y_fingerprint is not None else digest(Y),
              'sl2': [[item.start, item.stop] for item in sl2],
              'control': list(map(str, control_samples))}
    key = stage_key(stage, params=params, inputs=inputs, code=code_digest(sample_perm_pvalues, _stats_gam))
    root = _ensure_dir(root) if stage_cache_enabled(enabled) else None
    path = root / f'{stage}__{key}.npz' if root is not None else None
    if path is not None and _hit(path, key):
        if verbose:
            print(f'[stage cache] hit {stage} (key {key})')
        with np.load(path, allow_pickle=False) as stored:
            return (stored['shape_p'], stored['level_p'],
                    [tuple(group) for group in stored['splits'].tolist()], int(stored['true_idx']))
    if verbose:
        print(f'[stage cache] computing {stage} (key {key})')
    shape_p, level_p, splits, true_idx = sample_perm_pvalues(
        Y, s, samples, knots, grid, lambda_grid, lam_idx, sl2, p_b, control_samples)
    if path is not None:
        np.savez_compressed(path, shape_p=np.asarray(shape_p), level_p=np.asarray(level_p),
                            splits=np.asarray([list(map(str, split)) for split in splits]),
                            true_idx=np.asarray(true_idx))
        _record(path, stage=stage, key=key, params=params, inputs={k: repr(v)[:80] for k, v in inputs.items()},
                code=code_digest(sample_perm_pvalues, _stats_gam))
    return shape_p, level_p, splits, true_idx


def cached_neighbor_graph(adata, key, compute, *, stage, root, params=None, inputs=None,
                          code=None, enabled=None, verbose=True):
    """``compute()`` builds ``sc.pp.neighbors(adata, key_added=key, ...)`` once, then reuses it.

    A neighbour graph is the most expensive object in the pipeline to build and the cheapest to
    store, so the connectivity and distance matrices are cached together with the ``uns`` parameters
    Scanpy needs to find them again. On a hit the two ``obsp`` entries and the ``uns`` record are
    restored, so downstream PAGA, diffusion map and DPT calls behave exactly as after a fresh build.
    """
    payload_dir = _ensure_dir(root) if stage_cache_enabled(enabled) else None
    digest_of = code_digest(code)
    cache_key = stage_key(stage, params=params,
                          inputs={**(inputs or {}), 'key': key, 'n_obs': int(adata.n_obs)},
                          code=digest_of)
    path = payload_dir / f'{stage}__{cache_key}.npz' if payload_dir is not None else None
    if path is not None and _hit(path, cache_key):
        if verbose:
            print(f'[stage cache] hit {stage} (key {cache_key})')
        with np.load(path, allow_pickle=False) as stored:
            connectivities = sp.csr_matrix(
                (stored['conn_data'], stored['conn_indices'], stored['conn_indptr']),
                shape=tuple(stored['shape']))
            distances = sp.csr_matrix(
                (stored['dist_data'], stored['dist_indices'], stored['dist_indptr']),
                shape=tuple(stored['shape']))
            uns_record = json.loads(str(stored['uns_json']))
        # Scanpy stores the connectivities and distances under the keys named in `uns[key]`
        # (`<key>_connectivities` / `<key>_distances` for a vanilla `sc.pp.neighbors` call), so the
        # restored entries have to use those names rather than the neighbours key itself.
        adata.obsp[uns_record.get('connectivities_key', key)] = connectivities
        adata.obsp[uns_record.get('distances_key', f'{key}_distances')] = distances
        adata.uns[key] = uns_record
        return adata
    if verbose:
        print(f'[stage cache] computing {stage} (key {cache_key})')
    compute()
    if path is not None:
        record = dict(adata.uns.get(key, {}))
        connectivities_key = record.get('connectivities_key', key)
        distances_key = record.get('distances_key', f'{key}_distances')
        connectivities = adata.obsp[connectivities_key].tocsr()
        distances = (adata.obsp[distances_key].tocsr() if distances_key in adata.obsp
                     else connectivities.copy())
        np.savez_compressed(
            path,
            conn_data=connectivities.data, conn_indices=connectivities.indices,
            conn_indptr=connectivities.indptr,
            dist_data=distances.data, dist_indices=distances.indices,
            dist_indptr=distances.indptr,
            shape=np.asarray(connectivities.shape),
            uns_json=np.asarray(json.dumps({
                'connectivities_key': record.get('connectivities_key', key),
                'distances_key': distances_key,
                'params': {k: (v.item() if hasattr(v, 'item') else v)
                           for k, v in record.get('params', {}).items()},
            })),
        )
        _record(path, stage=stage, key=cache_key, params=params, inputs=inputs, code=digest_of)
    return adata


def cache_status(root) -> pd.DataFrame:
    """One row per cached entry: what is stored, when, and how large."""
    root = Path(root)
    if not root.exists():
        return pd.DataFrame(columns=['stage', 'key', 'written_utc', 'size_mb'])
    rows = []
    for sidecar in sorted(root.glob(f'*{_SIDECAR_SUFFIX}')):
        payload = Path(str(sidecar)[: -len(_SIDECAR_SUFFIX)])
        try:
            record = json.loads(sidecar.read_text())
        except json.JSONDecodeError:
            continue
        rows.append({
            'stage': record.get('stage', payload.stem),
            'key': record.get('key', ''),
            'written_utc': record.get('written_utc', ''),
            'size_mb': round(payload.stat().st_size / 1e6, 4) if payload.exists() else np.nan,
        })
    frame = pd.DataFrame(rows, columns=['stage', 'key', 'written_utc', 'size_mb'])
    return frame.sort_values('stage').reset_index(drop=True)


def purge_stage_cache(root, stage=None) -> int:
    """Delete cached entries (all of them, or only those whose name contains ``stage``)."""
    root = Path(root)
    if not root.exists():
        return 0
    removed = 0
    for sidecar in sorted(root.glob(f'*{_SIDECAR_SUFFIX}')):
        payload = Path(str(sidecar)[: -len(_SIDECAR_SUFFIX)])
        name = payload.stem.split('__')[0]
        if stage is not None and stage not in name:
            continue
        for path in (payload, sidecar):
            if path.exists():
                path.unlink()
                removed += 1
    return removed
