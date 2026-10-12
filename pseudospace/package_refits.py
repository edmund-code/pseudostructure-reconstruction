"""Refits of the package pseudospace on perturbed stage-02 inputs (notebook 72).

The ``pseudospace_reconstruction`` package (stages 01-04) is the frozen pseudospace. Notebook 72 checks its coordinate
the way notebooks 35-36 checked SCF13: perturb the counts (gene folds, count splitting, equal depth, an injected group
effect) and rebuild the coordinate. The package is never modified. A refit is a new run folder holding

- ``02_qc/units_qc.h5ad``: the base run's QC-retained units (same obs, var and order) with perturbed counts, and the
  package's own ``normalize_log1p`` recomputing X / ``layers['lognorm']``; gene folds set
  ``var['exclude_from_harmony_hvg']`` instead (notebook 35's withholding);
- ``02_qc/manifest.json``: the base manifest plus the perturbation record (stage 03 only checks that it exists);

after which the package's command line runs stages 03 and 04 on it. Stages 01-02 (units, QC) are held fixed. Unit
ids survive stages 03-04, so refit positions are matched back by id; a unit the refit drops from pass 2 gets NaN.
``extract_refit`` keeps what notebook 72 needs and ``cleanup_refit`` then deletes the large h5ad files.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

LARGE_FILES = ('02_qc/units_qc.h5ad', '03_integration/pass1_labelled.h5ad', '03_integration/pass2.h5ad')
EXCLUDE_COL = 'exclude_from_harmony_hvg'
TABLE_COLUMNS = ('sample', 'segment_class', 'coarse_class', 'pseudospace', 'curve_coordinate', 'marker_axis',
                 'dpt_pseudotime')


def _integer_counts(counts):
    counts = sparse.csr_matrix(counts)
    data = counts.data
    if not np.isfinite(data).all() or (data < 0).any() or not np.allclose(data, np.round(data)):
        raise ValueError('Expected nonnegative integer counts.')
    return counts


def perturbed_units(base, *, counts=None, exclude_genes=(), multiply=False, target_sum=1e4):
    """A stage-02 object like ``base`` with ``counts`` (default: the base counts) and lognorm rebuilt.

    ``normalize_log1p`` is the package's own stage-02 normalisation. ``multiply=True`` is the numerical floor:
    log1p(x * (target / total)) instead of scanpy's division, a change of order 1e-16. ``exclude_genes`` are flagged in
    ``var['exclude_from_harmony_hvg']`` so stage 03 keeps them out of its integration features.
    """
    import anndata as ad
    from pseudospace_reconstruction.qc import normalize_log1p

    template = base.layers['counts']
    raw = _integer_counts(template if counts is None else counts).astype(template.dtype)
    if raw.shape != base.shape:
        raise ValueError(f'Counts have shape {raw.shape}, the base units {base.shape}.')
    out = ad.AnnData(X=raw.copy(), obs=base.obs.copy(), var=base.var.copy(),
                     uns={k: v for k, v in base.uns.items() if k != 'log1p'},
                     obsm={k: np.asarray(v).copy() for k, v in base.obsm.items()})
    if multiply:
        out.layers['counts'] = out.X.copy()
        total = np.asarray(out.X.sum(axis=1), dtype=np.float64).ravel()
        if (total <= 0).any():
            raise ValueError('A unit has no counts.')
        scaled = sparse.csr_matrix(out.X.multiply((target_sum / total)[:, None]), dtype=template.dtype)
        scaled.data = np.log1p(scaled.data)
        out.X = scaled
        out.layers['lognorm'] = scaled.copy()
        out.uns['log1p'] = {'base': None}
    else:
        normalize_log1p(out, target_sum)
    # The package's stage-02 file stores lognorm with sorted column indices (X and counts keep the counts' order).
    # Stage 03 is chaotic: a different summation order alone changes its clusters, so match that layout exactly.
    lognorm = sparse.csr_matrix(out.layers['lognorm'], copy=True)
    lognorm.sort_indices()
    out.layers['lognorm'] = lognorm
    exclude = pd.Index(list(map(str, exclude_genes)))
    if len(exclude):
        missing = exclude.difference(out.var_names)
        if len(missing):
            raise ValueError(f'{len(missing)} genes to withhold are not in the units, e.g. {list(missing[:3])}.')
        out.var[EXCLUDE_COL] = out.var[EXCLUDE_COL].to_numpy(bool) | out.var_names.isin(exclude)
    return out


def write_refit_input(units, base_run, refit_run, record):
    """Write ``units`` as ``refit_run/02_qc/units_qc.h5ad`` with the base stage-02 manifest plus ``record``."""
    refit_run = Path(refit_run)
    qc = refit_run / '02_qc'
    qc.mkdir(parents=True, exist_ok=True)
    units.write_h5ad(qc / 'units_qc.h5ad')
    manifest = json.loads((Path(base_run) / '02_qc' / 'manifest.json').read_text())
    manifest['notebook72_perturbation'] = record
    (qc / 'manifest.json').write_text(json.dumps(manifest, indent=1, default=str))
    return refit_run


THREAD_VARIABLES = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS')


def run_stages_03_04(config_dir, data_root, refit_run, *, threads=None, env=None):
    """Run the package's command line (``run 03 04``) on a refit folder; return its output, raise on failure.

    ``threads=None`` removes the thread-count variables, so the package runs with its libraries' defaults, as the
    base run did. Stage 03 is chaotic: a different BLAS thread count alone changes its clusters and positions.
    """
    command = [sys.executable, '-m', 'pseudospace_reconstruction.cli', 'run', '03', '04', '--config', str(config_dir),
               '--data-root', str(data_root), '--out', str(refit_run)]
    environment = {**os.environ, **(env or {})}
    for name in THREAD_VARIABLES:
        if threads:
            environment[name] = str(int(threads))
        else:
            environment.pop(name, None)
    result = subprocess.run(command, capture_output=True, text=True, env=environment)
    output = (result.stdout + result.stderr).strip()
    (Path(refit_run) / 'notebook72_cli.log').write_text(output + '\n')
    if result.returncode != 0:
        raise RuntimeError(f'pseudospace-recon run 03 04 failed for {refit_run} (exit {result.returncode}):\n'
                           + '\n'.join(output.splitlines()[-20:]))
    return output


def _text(values):
    """Fixed-width unicode array: the stage cache stores payloads without pickle, so no object arrays."""
    return np.asarray(list(map(str, values)), dtype=str)


def extract_refit(refit_run):
    """What notebook 72 keeps from a finished refit, as a dict of arrays (stage-cache payload, no object arrays)."""
    import anndata as ad

    refit_run = Path(refit_run)
    table = pd.read_parquet(refit_run / '04_pseudospace' / 'pseudospace.parquet')
    out = {'unit': _text(table.index)}
    for column in TABLE_COLUMNS:
        if column in table:
            values = table[column]
            out[column] = (_text(values) if values.dtype == object or str(values.dtype) == 'category'
                           else values.to_numpy(float))
    pass1 = ad.read_h5ad(refit_run / '03_integration' / 'pass1_labelled.h5ad', backed='r')
    out['pass1_unit'] = _text(pass1.obs_names)
    out['pass1_segment_class'] = _text(pass1.obs['segment_class'])
    out['pass1_coarse_class'] = _text(pass1.obs['coarse_class'])
    out['pass1_features'] = _text(pass1.var_names[pass1.var['highly_variable_for_harmony'].to_numpy(bool)])
    pass1.file.close()
    pass2 = ad.read_h5ad(refit_run / '03_integration' / 'pass2.h5ad', backed='r')
    out['pass2_features'] = _text(pass2.var_names[pass2.var['highly_variable_for_harmony'].to_numpy(bool)])
    pass2.file.close()
    for stage, folder in (('03', '03_integration'), ('04', '04_pseudospace')):
        out[f'manifest_{stage}'] = np.array(json.dumps(json.loads((refit_run / folder / 'manifest.json').read_text())
                                                       .get('summary', {}), default=str))
    return out


def cleanup_refit(refit_run):
    """Delete the large h5ad files of a refit (input, pass 1, pass 2); return the bytes freed."""
    freed = 0
    for name in LARGE_FILES:
        path = Path(refit_run) / name
        if path.is_file():
            freed += path.stat().st_size
            path.unlink()
    return freed


def free_gb(path):
    """Free disk space at ``path`` in GB."""
    return shutil.disk_usage(Path(path)).free / 1e9


def refit_position(payload, units, column='pseudospace'):
    """A refit's ``column`` for ``units`` (unit ids), NaN for units the refit dropped from pass 2."""
    values = pd.Series(np.asarray(payload[column], float), index=pd.Index(np.asarray(payload['unit']).astype(str)))
    return values.reindex(pd.Index(units).astype(str)).to_numpy(float)


def orient_like_export(dpt, segment):
    """Notebook 70's DPT rule: flip (1 - x) when the pooled PT-S1 median exceeds the PT-S3 median."""
    dpt = np.asarray(dpt, float)
    segment = np.asarray(segment).astype(str)
    s1, s3 = np.nanmedian(dpt[segment == 'PT-S1']), np.nanmedian(dpt[segment == 'PT-S3'])
    return (1 - dpt, True) if s1 > s3 else (dpt, False)


def runs_s1_to_s3(position, specimen, segment):
    """True when, in every specimen, the PT-S1 median lies below the PT-S3 median (NaN positions ignored)."""
    frame = pd.DataFrame({'z': np.asarray(position, float), 'specimen': np.asarray(specimen).astype(str),
                          'segment': np.asarray(segment).astype(str)}).dropna()
    medians = frame.groupby(['specimen', 'segment']).z.median().unstack()
    if not {'PT-S1', 'PT-S3'} <= set(medians.columns) or medians[['PT-S1', 'PT-S3']].isna().any().any():
        return False
    return bool((medians['PT-S1'] < medians['PT-S3']).all())
