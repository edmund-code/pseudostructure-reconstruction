"""Genes outside the one-to-one ortholog space along the package's elastic curve (notebook 86).

Protocol: ``docs/results/pt-non-one-to-one.md``. Notebook 62 asked these questions on the retired scFates coordinate;
its coordinate-free parts (HCOP orthogroups and the probe audit, ``paralog_zonation``) are reused unchanged. The helpers
here are the rules the protocol fixes:

* ``family_clean`` and ``family_sums``: a family unit is the summed counts of its members whose probes hit nothing
  outside the family, so cross-hybridisation inside a family cannot matter;
* ``species_zonation``: one species' per-specimen NB spline curves, the typical gene's shape, and a likelihood-ratio
  test of each gene against that shape (not against a flat line, which most genes reject);
* ``zonation_calls``: the zonated and strongly zonated rule on two specimens;
* ``family_pairs``: unordered paralog pairs inside families, the units of the co-positioning tests.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import chi2, false_discovery_control

from .nb_gam import nb_dispersion, nb_fit
from .stats_gam import make_gam_design

FLOOR_CP10K = .05          # notebook 62's pseudocount: log2(CP10k + 0.05) bounds the log curve of near-absent genes
ZONATED_Q, ZONATED_AMPLITUDE, STRONG_AMPLITUDE, REPLICATION_R = .01, 1., 2., .5   # notebook 62's rule
USABLE_DETECTION = .05     # genes that define the typical shape (notebook 81's usable rule)


def family_clean(inventory):
    """Per gene: probe-clean, or ambiguous only through members of its own orthogroup (same species).

    ``inventory`` is one species' notebook-62 gene inventory: gene, component (-1 = none), probe_status
    (clean / near / ambiguous / unverified) and offtarget_genes (';'-joined symbols). Inside a family unit the counts are
    summed, so a probe that also hits a family member moves counts within the unit, never into it from outside.
    """
    inv = inventory.set_index('gene') if 'gene' in inventory.columns else inventory
    members = {k: set(v) for k, v in pd.Series(inv.index, index=inv.index).groupby(inv.component.to_numpy())}
    out = {}
    for gene, row in inv.iterrows():
        status = str(row.probe_status)
        if status in ('clean', 'near'):
            out[gene] = True
            continue
        if status != 'ambiguous' or int(row.component) < 0:
            out[gene] = False
            continue
        targets = {g for g in str(row.offtarget_genes if pd.notna(row.offtarget_genes) else '').split(';') if g}
        out[gene] = bool(targets) and targets <= members.get(int(row.component), set())
    return pd.Series(out, name='family_clean').reindex(inv.index)


def family_sums(counts, genes, families):
    """Units x families sparse matrix of summed counts; ``families`` is {name: member genes} (all in ``genes``)."""
    index = pd.Index(genes)
    rows, cols = [], []
    for k, members in enumerate(families.values()):
        where = index.get_indexer(list(members))
        if (where < 0).any() or len(set(where)) != len(where) or not len(where):
            raise ValueError(f'Family {list(families)[k]}: members missing, repeated or empty.')
        rows.extend(where)
        cols.extend([k] * len(where))
    membership = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(index), len(families)))
    return sparse.csr_matrix(sparse.csr_matrix(counts, dtype=float) @ membership)


def _in_range(grid, position):
    lo, hi = np.quantile(position, [.01, .99])
    return (grid >= lo) & (grid <= hi)


def species_zonation(counts, library, position, specimen, knots, grid, *, usable_detection=USABLE_DETECTION,
                     floor_cp10k=FLOOR_CP10K, shapes=None, n_jobs=1, block_size=200):
    """One species' zonation along the curve, specimen by specimen.

    counts: units x genes (raw); library: per unit; position: in [0, 1]; specimen: two (or more) names; knots and grid:
    the cubic B-spline knots and evaluation grid. Steps, all at one ML dispersion per gene (design: specimen intercepts
    plus the shared basis):
    1. per specimen, an NB spline of the rate (offset log library); phi = log2(CP10k curve + ``floor_cp10k``);
    2. the typical shape: the median over usable genes (detection >= ``usable_detection`` in that specimen, converged,
       not separated) of their mean-centred phi (``published_findings.reference_curve`` on phi);
    3. per specimen, a likelihood-ratio test of the spline against that shape: both models carry the offset
       log library + ln 2 x shape(position), the null has an intercept only (df = basis columns - 1, chi-square:
       unit weights and a fixed dispersion);
    4. centred curves c = phi - mean(phi) - shape, in log2, NaN outside the specimen's [1st, 99th] percentile range.
    ``shapes`` ({specimen: grid}) replaces step 2 with given shapes, so that extra columns (family sums) are read
    against the genes' typical shape.
    Returns {'alpha', 'alpha_converged', 'shape': {name: grid}, 'phi', 'centred', 'lr', 'p', 'q' (BH within specimen),
    'converged', 'separated', 'detection': {name: genes}} with per-specimen arrays (genes x grid or genes).
    """
    x = sparse.csr_matrix(counts, dtype=float)
    lib = np.asarray(library, float)
    s = np.asarray(position, float)
    names = list(dict.fromkeys(np.asarray(specimen).astype(str)))
    spec = np.asarray(specimen).astype(str)
    if x.shape[0] != len(s) or lib.shape != s.shape or spec.shape != s.shape:
        raise ValueError('Counts, library, position and specimen must align.')
    if (lib <= 0).any() or not np.isfinite(s).all():
        raise ValueError('Need positive libraries and finite positions.')
    basis_units = make_gam_design(s, knots)
    base = make_gam_design(grid, knots)
    dummies = np.column_stack([(spec == n).astype(float) for n in names])
    run = {'block_size': block_size, 'n_jobs': n_jobs}
    disp = nb_dispersion(x, np.column_stack([dummies, basis_units[:, 1:]]), np.ones(len(s)), np.log(lib), **run)
    out = {'alpha': disp['alpha'], 'alpha_converged': disp['converged']}
    for key in ('shape', 'phi', 'centred', 'lr', 'p', 'q', 'converged', 'separated', 'detection'):
        out[key] = {}
    df = basis_units.shape[1] - 1
    for name in names:
        rows = spec == name
        xs, offset = x[rows], np.log(lib[rows])
        ones = np.ones(rows.sum())
        fit = nb_fit(xs, basis_units[rows], ones, offset, disp, **run)
        phi = np.log2(np.exp(fit['beta'] @ base.T) * 1e4 + floor_cp10k)
        inside = _in_range(grid, s[rows])
        detection = np.asarray((xs > 0).mean(axis=0)).ravel()
        if shapes is None:
            usable = fit['converged'] & ~fit['separated'] & (detection >= usable_detection)
            if usable.sum() < 10:
                raise ValueError(f'{name}: fewer than 10 usable genes define the typical shape.')
            y = phi[usable][:, inside]
            shape = np.full(len(grid), np.nan)
            shape[inside] = np.median(y - y.mean(axis=1, keepdims=True), axis=0)
        else:
            shape = np.asarray(shapes[name], float)
            if shape.shape != grid.shape or not np.isfinite(shape[inside]).all():
                raise ValueError(f'{name}: a given shape must be finite on the specimen range of the grid.')
        shifted = offset + np.log(2) * np.interp(s[rows], grid[inside], shape[inside])
        full = nb_fit(xs, basis_units[rows], ones, shifted, disp, **run)
        null = nb_fit(xs, np.ones((rows.sum(), 1)), ones, shifted, disp, **run)
        lr = np.maximum(null['deviance'] - full['deviance'], 0.)
        ok = fit['converged'] & full['converged'] & null['converged'] & np.isfinite(lr)
        p = np.where(ok, chi2.sf(lr, df), np.nan)
        q = np.full(len(p), np.nan)
        q[ok] = false_discovery_control(p[ok])
        centred = np.full(phi.shape, np.nan)
        c = phi[:, inside] - phi[:, inside].mean(axis=1, keepdims=True) - shape[inside]
        centred[:, inside] = c
        out['shape'][name], out['phi'][name], out['centred'][name] = shape, phi, centred
        out['lr'][name], out['p'][name], out['q'][name] = lr, p, q
        out['converged'][name] = ok
        out['separated'][name] = fit['separated'] | full['separated']
        out['detection'][name] = detection
    return out


def zonation_calls(zonation, genes, grid_units, *, q=ZONATED_Q, amplitude=ZONATED_AMPLITUDE, strong=STRONG_AMPLITUDE,
                   replication=REPLICATION_R):
    """Notebook 62's two-specimen rule on ``species_zonation`` output.

    Zonated: q < ``q`` in every specimen, amplitude (max - min of the centred curve) >= ``amplitude`` in every specimen,
    Pearson r between the two specimens' centred curves >= ``replication`` (on the grid points both cover), converged
    and not separated everywhere. Strong: zonated with amplitude >= ``strong`` in every specimen. Z = the smallest
    amplitude. Peak and trough: argmax / argmin of the specimens' mean centred curve, in ``grid_units``.
    """
    names = list(zonation['centred'])
    if len(names) != 2:
        raise ValueError('The rule needs exactly two specimens.')
    a, b = (zonation['centred'][n] for n in names)
    both = np.isfinite(a).all(axis=0) & np.isfinite(b).all(axis=0)
    if both.sum() < 5:
        raise ValueError('The specimens share fewer than 5 grid points.')
    amp = {n: np.nanmax(zonation['centred'][n], axis=1) - np.nanmin(zonation['centred'][n], axis=1) for n in names}
    ac, bc = a[:, both] - a[:, both].mean(1, keepdims=True), b[:, both] - b[:, both].mean(1, keepdims=True)
    den = np.sqrt((ac ** 2).sum(1) * (bc ** 2).sum(1))
    r = np.divide((ac * bc).sum(1), den, out=np.full(len(den), np.nan), where=den > 0)
    mean = (a[:, both] + b[:, both]) / 2
    units = np.asarray(grid_units, float)[both]
    ok = np.logical_and.reduce([zonation['converged'][n] & ~zonation['separated'][n] for n in names])
    sig = np.logical_and.reduce([np.nan_to_num(zonation['q'][n], nan=1.) < q for n in names])
    big = np.logical_and.reduce([amp[n] >= amplitude for n in names])
    table = pd.DataFrame({f'q_{n}': zonation['q'][n] for n in names}, index=pd.Index(genes, name='gene'))
    for n in names:
        table[f'amplitude_{n}'] = amp[n]
        table[f'detection_{n}'] = zonation['detection'][n]
    table['replication_r'] = r
    table['Z'] = np.minimum(*(amp[n] for n in names))
    table['peak'] = units[np.argmax(mean, axis=1)]
    table['trough'] = units[np.argmin(mean, axis=1)]
    table['usable'] = ok
    table['zonated'] = ok & sig & big & (np.nan_to_num(r, nan=-1.) >= replication)
    table['strong'] = table.zonated & (table.Z >= strong)
    return table


def family_pairs(members):
    """All unordered pairs inside each family: {family: [genes]} -> DataFrame family, a, b (a < b by name)."""
    rows = []
    for family, genes in members.items():
        g = sorted(set(genes))
        rows.extend((family, g[i], g[j]) for i in range(len(g)) for j in range(i + 1, len(g)))
    return pd.DataFrame(rows, columns=['family', 'a', 'b'])
