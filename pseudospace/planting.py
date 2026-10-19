"""Planted pathway-level differences for a ground-truth benchmark (notebook 100; docs/results/pt-planted-benchmark.md).

* ``choose_sets``: k sets of a size range, drawn in random order and accepted when their Jaccard with every set already
  chosen is below a limit and enough unplanted members remain; a gene is planted in one set at most.
* ``effect_shape``: the per-unit shape f(x) of a scenario: 'shape' and 'bump' are centred over the planted specimens'
  units (pure positional), 'located' is an S3-only logistic, 'level' is constant.
* ``plant_counts``: multiply expected counts by m: binomial thinning when m < 1, Poisson addition of c (m - 1) when
  m > 1. Ponytail: the added counts are Poisson, so increases gain slightly less overdispersion than an NB would;
  regenerating planted genes from the fitted NB is the upgrade.
"""
from __future__ import annotations

import numpy as np

SCENARIOS = ('shape', 'located', 'level', 'bump')


def choose_sets(sets, rng, *, k=15, sizes=(20, 150), max_jaccard=.1, frac=.6, min_members=10):
    """{set id: planted genes} for k sets; raises when the pool cannot supply k."""
    pool = [s for s, m in sets.items() if sizes[0] <= len(m) <= sizes[1]]
    order = rng.permutation(len(pool))
    chosen, used = {}, set()
    for i in order:
        s = pool[i]
        members = set(sets[s])
        if any(len(members & set(sets[c])) / len(members | set(sets[c])) >= max_jaccard for c in chosen):
            continue
        free = sorted(members - used)
        n = max(min_members, int(round(frac * len(members))))
        if len(free) < n:
            continue
        planted = sorted(rng.choice(free, n, replace=False).tolist())
        chosen[s] = planted
        used.update(planted)
        if len(chosen) == k:
            return chosen
    raise ValueError(f'Only {len(chosen)} sets could be chosen of {k}.')


def effect_shape(x, scenario, planted_units, *, centre=.65, width=.06, bump_centre=.5, bump_sd=.1):
    """f(x) per unit for a scenario; centred shapes have mean 0 over ``planted_units``."""
    x = np.asarray(x, float)
    planted_units = np.asarray(planted_units, bool)
    if scenario == 'level':
        return np.ones_like(x)
    if scenario in ('shape', 'located'):
        f = 1. / (1. + np.exp(-(x - centre) / width))
    elif scenario == 'bump':
        f = np.exp(-(x - bump_centre) ** 2 / (2 * bump_sd ** 2))
    else:
        raise ValueError(f'Unknown scenario {scenario!r}; one of {SCENARIOS}.')
    if scenario in ('shape', 'bump'):
        f = f - f[planted_units].mean()
    return f


def plant_counts(counts, multiplier, rng):
    """Counts (units x genes, whole numbers) with expectation scaled by ``multiplier`` (same shape, > 0)."""
    c = np.asarray(counts, float)
    m = np.asarray(multiplier, float)
    if c.shape != m.shape or (m <= 0).any() or (c < 0).any():
        raise ValueError('Need nonnegative counts and positive multipliers of the same shape.')
    out = c.copy()
    down = m < 1
    out[down] = rng.binomial(c[down].astype(np.int64), m[down])
    up = m > 1
    out[up] = c[up] + rng.poisson(c[up] * (m[up] - 1.))
    return out
