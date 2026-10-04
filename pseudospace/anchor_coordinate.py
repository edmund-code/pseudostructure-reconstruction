"""External segment zonation and conserved-anchor gene selection (notebooks 35 and later).

External data resolve only coarse PT segments: S1/S2/S3 in mouse and rat, convoluted PT versus S3
in human. A conserved anchor is a gene whose S3-versus-early direction agrees across species in
those data. Selection uses external tables plus position-free filters from our data only; it
never sees our segment labels, coordinate or test statistics.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

SEGMENTS = ('S1', 'S2', 'S3')


def bulk_segment_profiles(table, pattern):
    """log2(value + 1) averaged over replicate columns whose name contains ``pattern.format(seg=...)``.

    Adds ``early`` as the mean of S1 and S2 (both microdissected separately).
    """
    out = {}
    for seg in SEGMENTS:
        columns = [c for c in table.columns if pattern.format(seg=seg) in c]
        if not columns:
            raise ValueError(f'No replicate columns for {seg}.')
        values = table[columns].apply(pd.to_numeric, errors='coerce').clip(lower=0)
        out[seg] = np.log2(values + 1).mean(axis=1)
    frame = pd.DataFrame(out)
    frame['early'] = frame[['S1', 'S2']].mean(axis=1)
    return frame


def census_segment_profiles(counts, meta, segment_map, *, early=None, min_cells=50,
                            exclude_stage=r'embryo|E16|Theiler'):
    """Per donor and sex, log2 CPM of counts pooled within each segment.

    ``meta`` rows align with ``counts`` rows and carry donor_id, sex, cell_type, development_stage,
    n_cells and total_counts_all_genes (the cells' totals over all genes, so normalization does not
    depend on the gene panel). Groups below ``min_cells`` are dropped, and a donor is kept only with
    every mapped segment present. ``early`` lists segments pooled by counts into an extra column.
    """
    meta = meta.loc[~meta.development_stage.astype(str).str.contains(exclude_stage, case=False)].copy()
    meta['segment'] = meta.cell_type.map(segment_map)
    meta = meta[meta.segment.notna() & meta.n_cells.ge(min_cells)]
    needed = set(segment_map.values())
    complete = meta.groupby('donor_id').segment.apply(lambda s: needed.issubset(set(s)))
    meta = meta[meta.donor_id.isin(complete[complete].index)]
    columns = {}
    for (donor, sex), block in meta.groupby(['donor_id', 'sex']):
        groups = {seg: block.index[block.segment.eq(seg)] for seg in needed}
        if early:
            groups['early'] = block.index[block.segment.isin(early)]
        for seg, idx in groups.items():
            pooled = counts.loc[idx].sum(axis=0)
            columns[(donor, sex, seg)] = np.log2(pooled / meta.loc[idx, 'total_counts_all_genes'].sum() * 1e6 + 1)
    frame = pd.DataFrame(columns)
    if frame.empty:
        raise ValueError('No donor passed the segment and cell-count rules.')
    frame.columns.names = ['donor', 'sex', 'segment']
    return frame


def segment_contrast(frame, *, late='S3', early='early', sex=None, donors=None):
    """Return (late - early, max(late, early)) averaged over the selected donor columns."""
    keep = np.ones(frame.shape[1], bool)
    if sex is not None:
        keep &= frame.columns.get_level_values('sex') == sex
    if donors is not None:
        keep &= frame.columns.get_level_values('donor').isin(list(donors))
    sub = frame.loc[:, keep]
    if sub.shape[1] == 0:
        raise ValueError('No external columns match the requested donors or sex.')
    late_mean = sub.xs(late, axis=1, level='segment').mean(axis=1)
    early_mean = sub.xs(early, axis=1, level='segment').mean(axis=1)
    return late_mean - early_mean, pd.concat([late_mean, early_mean], axis=1).max(axis=1)


def conserved_anchor_mask(contrasts, levels, eligible, *, threshold, min_level=1.0):
    """Genes with |contrast| >= threshold and the same sign in every dataset column.

    ``contrasts`` and ``levels`` are genes x datasets (same columns); ``eligible`` is a boolean
    gene mask of position-free filters (measured, detected, PT-specific).
    """
    if list(contrasts.columns) != list(levels.columns) or not contrasts.index.equals(levels.index):
        raise ValueError('Contrasts and levels must share genes and dataset columns.')
    eligible = pd.Series(eligible, index=contrasts.index).astype(bool)
    finite = contrasts.notna().all(axis=1) & levels.notna().all(axis=1)
    signs = np.sign(contrasts.fillna(0))
    same = signs.nunique(axis=1).eq(1) & signs.iloc[:, 0].ne(0)
    big = contrasts.abs().ge(threshold).all(axis=1)
    expressed = levels.ge(min_level).all(axis=1)
    return eligible & finite & same & big & expressed


def pt_specific_mask(class_means, *, pt='PT', max_ratio=2.0):
    """False for genes whose mean in any non-PT class exceeds ``max_ratio`` x the PT mean.

    ``class_means`` maps species to a genes x classes frame of linear-scale means. A gene must pass
    in every species. This blocks placement driven by neighbouring-segment contamination.
    """
    masks = []
    for frame in class_means.values():
        other = frame.drop(columns=pt).max(axis=1)
        masks.append(other <= max_ratio * frame[pt].clip(lower=1e-12))
    return pd.concat(masks, axis=1).all(axis=1)


def positive_control_masks(human, mouse, mouse_check, human_level, eligible, *,
                           discordant=.5, flat_mouse=1.0, flat_human=.25, min_level=1.0):
    """Externally species-dependent zonation: opposite-sign and mouse-zonated/human-flat genes.

    ``mouse_check`` is a second mouse contrast whose sign must agree with ``mouse``.
    """
    eligible = pd.Series(eligible, index=human.index).astype(bool)
    finite = human.notna() & mouse.notna() & mouse_check.notna() & human_level.notna()
    mouse_agree = np.sign(mouse) == np.sign(mouse_check)
    opposite = (eligible & finite & mouse_agree & human.abs().ge(discordant) & mouse.abs().ge(discordant)
                & (np.sign(human) != np.sign(mouse)) & human_level.ge(min_level))
    flat = (eligible & finite & mouse_agree & mouse.abs().ge(flat_mouse) & mouse_check.abs().ge(flat_mouse)
            & human.abs().lt(flat_human) & human_level.ge(min_level))
    return opposite, flat
