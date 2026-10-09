"""One standardized gene model under two pathway pipelines (rank-first and clustering-first).

The gene model is the frozen NB specification (docs/results/pt-gene-model-nb.md), built from ``nb_gam``: an NB GLM
on raw counts with log link and offset log(library), a cubic B-spline basis (6 df plus intercept, internal knots at
the position quartiles), contrast_designs' prior weights w_i = n / (2 J_s n_j), and one ML dispersion per gene under
the species-split M_full, reused for every fit. ``standard_gene_model`` returns everything both pipelines consume:
- the NB likelihood ratios T_level, T_spatial and T_total of one two-group comparison (the rank-first input);
- group curves and the group difference Delta(s) with its HC3 SE, natural-log rates on a common grid;
- one curve per specimen on the same grid (the clustering-first input: per-sample curves and their gaps).

A relabeled grouping (``pathway_calibration.balanced_partitions``) keeps the species' own smooth shape as a nuisance
term in every model and the species-split dispersion, exactly as notebooks 37 and 65 do; ``partition_groups``
builds those inputs. ``rank_first_gsea`` is a thin wrapper around gseapy's preranked GSEA.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata

from .nb_gam import nb_dispersion, nb_fit, nested_nb_lr
from .pathway_calibration import COMPARISONS, balanced_partitions, contrast_designs
from .stats_gam import gam_internal_knots, make_gam_design

STATISTICS = ('T_level', 'T_spatial', 'T_total')
SUPPORT = (.01, .99)          # each specimen's own position range: its 1st-99th percentiles


def specimen_ranges(position, specimen, quantiles=SUPPORT):
    """Each specimen's [1st, 99th] percentile position interval: DataFrame (lo, hi) indexed by specimen."""
    s = pd.Series(np.asarray(position, float))
    groups = s.groupby(np.asarray(specimen).astype(str))
    return pd.DataFrame({'lo': groups.quantile(quantiles[0]), 'hi': groups.quantile(quantiles[1])})


def common_grid(position, specimen, n_points=101, quantiles=SUPPORT):
    """``n_points`` equally spaced positions over the intersection of the specimens' [1st, 99th] percentile ranges."""
    ranges = specimen_ranges(position, specimen, quantiles)
    lo, hi = ranges.lo.max(), ranges.hi.min()
    if not lo < hi:
        raise ValueError("The specimens' position ranges do not overlap.")
    return np.linspace(lo, hi, int(n_points))


def partition_groups(specimen, human):
    """{label: (group, nuisance_shape)} per structure for the species split and its two balanced relabelings.

    Notebooks 37/65: the species split is group = human with no nuisance; a relabeling puts the specimens of
    ``balanced_partitions`` in group 1 and keeps species as the nuisance shape.
    """
    specimen, human = np.asarray(specimen).astype(str), np.asarray(human, bool)
    species = human.astype(float)
    return {label: (species, None) if label == 'species' else (np.isin(specimen, members).astype(float), species)
            for label, members in balanced_partitions(specimen, np.where(human, 'human', 'mouse')).items()}


@dataclass(frozen=True)
class GeneModel:
    """Output of ``standard_gene_model``. Curves are genes x grid natural-log rates per library unit (add log 1e4
    for log CP10k). Group curves are the equal-specimen mean of each group's specimens; Delta = case - reference.
    Specimen curves are NaN outside that specimen's own [1st, 99th] percentile range. Values of fits that did not
    converge or are separated are kept: mask them with the flags."""
    statistics: pd.DataFrame            # genes x T_level, T_spatial, T_total (NB LR), alpha, converged, separated
    grid: np.ndarray
    knots: np.ndarray
    reference: np.ndarray               # group 0
    case: np.ndarray                    # group 1
    delta: np.ndarray
    se: np.ndarray                      # HC3 SE of delta
    specimen_curves: dict               # specimen -> genes x grid
    specimen_converged: pd.DataFrame    # genes x specimens
    specimen_separated: pd.DataFrame
    groups: pd.Series                   # specimen -> 0/1
    dispersion: dict                    # pass to the relabeled runs: alpha and converged (nb_dispersion's output)


def standard_gene_model(counts, library, position, group, specimen, grid=None, *, genes=None, nuisance_shape=None,
                        dispersion=None, basis_df=6, block_size=200, n_jobs=1, fill='median'):
    """The frozen NB gene model of one two-group comparison (group 1 = case, group 0 = reference).

    counts: structures x genes raw counts (dense or sparse); library: per-structure library size (offset log
    library); position: coordinate in [0, 1]; group: 0/1 per structure (whole specimens); grid: defaults to
    ``common_grid(position, specimen)``. ``nuisance_shape``: the species 0/1 when ``group`` is a relabeling (see
    ``partition_groups``). ``dispersion``: nb_dispersion's dict or one alpha per gene; when None it is estimated under
    the M_full of the species split, which is ``nuisance_shape`` when given and ``group`` otherwise. Statistics of
    genes that did not converge are filled by ``fill`` (nb_gam's rule); separated genes are flagged, never filled.
    """
    s, g = np.asarray(position, float), np.asarray(group, float)
    samples = np.asarray(specimen).astype(str)
    lib = np.asarray(library, float)
    n = len(s)
    if counts.ndim != 2 or counts.shape[0] != n or g.shape != (n,) or samples.shape != (n,) or lib.shape != (n,):
        raise ValueError('Count rows, library, position, group and specimen must align.')
    if not np.isfinite(lib).all() or (lib <= 0).any() or not np.isfinite(s).all() or s.min() < 0 or s.max() > 1:
        raise ValueError('Need positive finite library sizes and positions in [0, 1].')
    grid = common_grid(s, samples) if grid is None else np.asarray(grid, float)
    if grid.ndim != 1 or grid.size < 2 or np.any(np.diff(grid) <= 0) or grid.min() < s.min() or grid.max() > s.max():
        raise ValueError('Need an increasing grid within the observed coordinate range.')
    genes = pd.RangeIndex(counts.shape[1]) if genes is None else pd.Index(genes)
    if len(genes) != counts.shape[1]:
        raise ValueError('One gene name per count column.')
    knots = gam_internal_knots(s, basis_df=basis_df)
    designs, weights = contrast_designs(s, g, samples, knots, nuisance_shape=nuisance_shape)
    offset = np.log(lib)
    run = {'block_size': block_size, 'n_jobs': n_jobs}
    if dispersion is None:
        split = g if nuisance_shape is None else np.asarray(nuisance_shape, float)
        species_designs, species_weights = contrast_designs(s, split, samples, knots)
        dispersion = nb_dispersion(counts, species_designs['full'], species_weights, offset, **run)
    if not isinstance(dispersion, dict):
        dispersion = {'alpha': np.asarray(dispersion, float), 'converged': np.ones(counts.shape[1], bool)}
    lr = nested_nb_lr(counts, designs, weights, offset, dispersion, {k: COMPARISONS[k] for k in STATISTICS},
                      fill=fill, **run)

    # Group curves: the full design's rows at each group's equal-specimen mean. Column layout (contrast_designs):
    # [intercept, basis, group, group x basis, within-group specimen contrasts (mean 0 in each group),
    #  nuisance x basis (last p_b - 1 columns, at the group's equal-specimen nuisance mean)].
    base = make_gam_design(grid, knots)
    p_b = base.shape[1]
    groups = pd.Series(g, index=samples).groupby(level=0).first()
    rows = {}
    for level in (0., 1.):
        r = np.zeros((len(grid), designs['full'].shape[1]))
        r[:, :p_b], r[:, p_b], r[:, p_b + 1:2 * p_b] = base, level, level * base[:, 1:]
        if nuisance_shape is not None:
            share = pd.Series(np.asarray(nuisance_shape, float), index=samples).groupby(level=0).mean()
            r[:, -(p_b - 1):] = share[groups.index[groups.eq(level)]].mean() * base[:, 1:]
        rows[level] = r
    full = nb_fit(counts, designs['full'], weights, offset, dispersion,
                  contrasts={'reference': rows[0.], 'case': rows[1.], 'delta': rows[1.] - rows[0.]}, **run)
    curves = full['contrasts']

    # Specimen curves: intercept + the same basis on that specimen's structures, unit weights, the gene's alpha.
    ranges = specimen_ranges(s, samples)
    own, converged, separated = {}, {}, {}
    for name, (lo, hi) in ranges.iterrows():
        rows_j = samples == name
        fit = nb_fit(counts[rows_j], make_gam_design(s[rows_j], knots), np.ones(rows_j.sum()), offset[rows_j],
                     dispersion, **run)
        own[name] = fit['beta'] @ base.T
        own[name][:, (grid < lo) | (grid > hi)] = np.nan
        converged[name], separated[name] = fit['converged'], fit['separated']

    statistics = pd.DataFrame({k: lr[k] for k in STATISTICS}, index=genes)
    statistics['alpha'] = np.asarray(dispersion['alpha'], float)
    statistics['converged'] = lr['converged'] & full['converged']
    statistics['separated'] = lr['separated'] | full['separated']
    return GeneModel(statistics=statistics, grid=grid, knots=knots, reference=curves['reference']['value'],
                     case=curves['case']['value'], delta=curves['delta']['value'], se=curves['delta']['se_hc3'],
                     specimen_curves=own, specimen_converged=pd.DataFrame(converged, index=genes),
                     specimen_separated=pd.DataFrame(separated, index=genes), groups=groups.astype(int),
                     dispersion=dispersion)


def rank_normal_scores(statistic):
    """Rank-normal score Phi^-1((rank - 0.5) / G) of a gene statistic (rank 1 = smallest; ties share their mean rank).

    The pre-specified sensitivity metric for preranked GSEA: the NB LR is heavy-tailed, and GSEA's weighted running
    sum would otherwise be driven by a few extreme genes.
    """
    values = pd.Series(statistic, dtype=float)
    if not np.isfinite(values.to_numpy()).all():
        raise ValueError('The statistic must be finite.')
    return pd.Series(norm.ppf((rankdata(values) - .5) / len(values)), index=values.index, name=values.name)


def rank_first_gsea(statistic, gene_sets, *, metric='raw', min_size=10, max_size=300, permutation_num=1000, seed=12,
                    threads=4, fdr=.05):
    """Rank all genes by ``statistic`` (descending) and run gseapy's preranked GSEA on ``gene_sets``.

    gseapy defaults (weight 1, gene-set permutation, its own FDR) except sizes, permutations, seed and threads; no
    plots or files. ``metric='rank_normal'`` ranks by ``rank_normal_scores(statistic)`` instead of the raw values.
    Genes are passed sorted by (score descending, name) so ties are ordered reproducibly. Returns one row per tested
    pathway, in ``gene_sets`` order: pathway_id, n (genes in the ranking), ES, NES, p (GSEA's nominal p, computed
    within the sign of ES), q (GSEA's FDR), fwer_p, leading_edge (list), p_one_sided (upper tail, notebook 61's
    convention: p/2 when NES > 0, else 1 - p/2; exact only for a symmetric null ES) and call (NES > 0, q <= fdr).
    """
    import gseapy as gp

    values = pd.Series(statistic, dtype=float)
    if values.index.has_duplicates or not np.isfinite(values.to_numpy()).all():
        raise ValueError('Need one finite score per uniquely named gene.')
    if metric == 'rank_normal':
        values = rank_normal_scores(values)
    elif metric != 'raw':
        raise ValueError("metric must be 'raw' or 'rank_normal'.")
    ranking = (values.rename_axis('gene').reset_index(name='score').astype({'gene': str})
               .sort_values(['score', 'gene'], ascending=[False, True], kind='stable'))
    result = gp.prerank(rnk=ranking, gene_sets=gene_sets, min_size=min_size, max_size=max_size,
                        permutation_num=permutation_num, seed=seed, threads=threads, outdir=None, no_plot=True)
    table = result.res2d.set_index(result.res2d.Term.astype(str))
    universe = set(ranking.gene)
    tested = [p for p in gene_sets if p in table.index]
    table = table.loc[tested]
    out = pd.DataFrame({'pathway_id': tested, 'n': [len(universe.intersection(map(str, gene_sets[p]))) for p in tested],
                        'ES': pd.to_numeric(table.ES).to_numpy(), 'NES': pd.to_numeric(table.NES).to_numpy(),
                        'p': pd.to_numeric(table['NOM p-val']).to_numpy(),
                        'q': pd.to_numeric(table['FDR q-val']).to_numpy(),
                        'fwer_p': pd.to_numeric(table['FWER p-val']).to_numpy(),
                        'leading_edge': [[x for x in str(v).split(';') if x] for v in table.Lead_genes]})
    out['p_one_sided'] = np.where(out.NES > 0, out.p / 2, 1 - out.p / 2)
    out['call'] = out.NES.gt(0) & out.q.le(fdr)
    return out
