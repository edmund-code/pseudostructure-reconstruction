"""Uncertainty of reconstructed positions: spatially blocked delete-half jackknife.

The pipeline (integration + elastic curve) is rerun on half-samples of each section. Resampling is by spatial
block because neighbouring tubules share position; half-sampling (not the bootstrap) because duplicated units
would add artificial neighbours to the kNN graph. With d = n/2 deleted, the delete-d jackknife variance of a
full-data statistic is the mean squared deviation of the half-sample statistics (Shao & Wu 1989).

Protocol: docs/results/pt-position-uncertainty.md.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd


def spatial_blocks(xy, section, per_block: float = 25.0) -> np.ndarray:
    """Square-grid block id per unit ("section|i|j"); the cell edge gives `per_block` units per occupied cell
    on average, set from each section's own density so the pixel scale does not matter."""
    xy = np.asarray(xy, float)
    section = np.asarray(section).astype(str)
    out = np.empty(len(xy), dtype=object)
    for s in np.unique(section):
        k = section == s
        p = xy[k]
        lo = p.min(0)
        area = np.prod(np.ptp(p, axis=0)) or 1.0
        edge = np.sqrt(area * per_block / k.sum())
        for _ in range(20):  # occupied cells hold more units than the bounding-box average; shrink to target
            ij = np.floor((p - lo) / edge).astype(int)
            mean_occ = k.sum() / len({(a, b) for a, b in ij})
            if mean_occ <= per_block * 1.05:
                break
            edge *= np.sqrt(per_block / mean_occ)
        out[k] = [f'{s}|{a}|{b}' for a, b in ij]
    return out


def half_sample(blocks, section, rng) -> np.ndarray:
    """Boolean mask keeping floor(n_blocks / 2) uniformly chosen blocks of every section, all units in them."""
    blocks = np.asarray(blocks).astype(str)
    section = np.asarray(section).astype(str)
    keep = np.zeros(len(blocks), bool)
    for s in np.unique(section):
        names = np.unique(blocks[section == s])
        chosen = rng.choice(names, size=len(names) // 2, replace=False)
        keep |= (section == s) & np.isin(blocks, chosen)
    return keep


def jackknife_se(draws, full=None) -> float:
    """Delete-half jackknife SE of a statistic from its half-sample values (NaN ignored)."""
    v = np.asarray(draws, float)
    v = v[np.isfinite(v)]
    if len(v) < 2:
        return np.nan
    return float(np.sqrt(np.mean((v - v.mean()) ** 2)))


def interval(full: float, draws, z: float = 1.96) -> dict:
    """Full-data estimate, jackknife SE, normal interval and the number of usable draws."""
    se = jackknife_se(draws)
    return {'estimate': full, 'se': se, 'lo': full - z * se, 'hi': full + z * se,
            'n_draws': int(np.isfinite(np.asarray(draws, float)).sum())}


def unit_summary(positions: pd.DataFrame) -> pd.DataFrame:
    """Per-unit summary of a units x draws position table (NaN where a unit was not drawn)."""
    p = positions.to_numpy(float)
    n = np.isfinite(p).sum(1)
    with warnings.catch_warnings():   # units never drawn give all-NaN rows; they stay NaN
        warnings.simplefilter('ignore', RuntimeWarning)
        mean = np.nanmean(p, 1)
        se = np.sqrt(np.nanmean((p - mean[:, None]) ** 2, 1))
        q = np.nanquantile(p, [.025, .975], axis=1)
    return pd.DataFrame({'n_draws': n, 'mean': mean, 'se': np.where(n >= 2, se, np.nan), 'q025': q[0], 'q975': q[1]},
                        index=positions.index)


def call_share(draws, predicate) -> float:
    """Share of usable draws in which `predicate(value)` holds."""
    v = [x for x in draws if x is not None and not (isinstance(x, float) and np.isnan(x))]
    return float(np.mean([bool(predicate(x)) for x in v])) if v else np.nan
