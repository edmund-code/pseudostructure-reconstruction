"""Pseudospace-binned marker expression heatmaps and their scaling/densify helpers.

Extracted verbatim from ``6_mouse_only_pseudospace.ipynb`` (Section 3).
``plot_marker_heatmap`` takes explicit ``output_dir``/``project_dir`` arguments
instead of reading notebook globals.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy import sparse
from scipy.ndimage import gaussian_filter1d

from .markers import resolve_available_marker_groups


def matrix_to_dense(matrix):
    if sparse.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix)


def robust_scale_rows(values):
    values = np.asarray(values, dtype=float)
    mean = np.nanmean(values, axis=1, keepdims=True)
    std = np.nanstd(values, axis=1, keepdims=True)
    std[~np.isfinite(std) | (std == 0)] = 1.0
    return (values - mean) / std


def binned_interpolated_expression(adata_obj, gene_idx, pseudotime_col, n_bins=120, smooth_sigma=2.5, expression_layer='lognorm'):
    pos = adata_obj.obs[pseudotime_col].to_numpy(dtype=float)
    finite = np.isfinite(pos)
    if finite.sum() < 10:
        raise ValueError(f'Too few finite {pseudotime_col} values to build heatmap')
    adata_obj = adata_obj[finite].copy()
    pos = pos[finite]

    x_min = float(np.nanmin(pos))
    x_max = float(np.nanmax(pos))
    if np.isclose(x_min, x_max):
        x_min, x_max = x_min - 0.005, x_max + 0.005

    sort_idx = np.argsort(pos)
    pos = pos[sort_idx]
    X_obj = adata_obj.layers[expression_layer] if expression_layer in adata_obj.layers else adata_obj.X
    expr = matrix_to_dense(X_obj[:, gene_idx])[sort_idx].T

    bin_edges = np.linspace(x_min, x_max, n_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_ids = np.digitize(pos, bin_edges) - 1
    bin_ids = np.clip(bin_ids, 0, n_bins - 1)

    heatmap = np.full((len(gene_idx), n_bins), np.nan)
    occupied = np.zeros(n_bins, dtype=bool)
    for bin_i in range(n_bins):
        mask = bin_ids == bin_i
        count = int(mask.sum())
        if count > 0:
            occupied[bin_i] = True
            heatmap[:, bin_i] = expr[:, mask].mean(axis=1)

    occupied_idx = np.flatnonzero(occupied)
    if len(occupied_idx) == 0:
        raise ValueError(f'No occupied bins for {pseudotime_col}')
    first_bin = int(occupied_idx[0])
    last_bin = int(occupied_idx[-1])

    for row_i in range(heatmap.shape[0]):
        valid = np.isfinite(heatmap[row_i])
        valid[:first_bin] = False
        valid[last_bin + 1:] = False
        if valid.sum() == 0:
            continue
        if valid.sum() == 1:
            heatmap[row_i, first_bin:last_bin + 1] = heatmap[row_i, valid][0]
        else:
            heatmap[row_i, first_bin:last_bin + 1] = np.interp(
                bin_centers[first_bin:last_bin + 1],
                bin_centers[valid],
                heatmap[row_i, valid],
            )

    heatmap = heatmap[:, first_bin:last_bin + 1]
    bin_centers = bin_centers[first_bin:last_bin + 1]
    bin_edges = bin_edges[first_bin:last_bin + 2]

    if smooth_sigma and smooth_sigma > 0:
        heatmap = gaussian_filter1d(heatmap, sigma=smooth_sigma, axis=1, mode='nearest')

    gene_mean = np.nanmean(heatmap, axis=1, keepdims=True)
    gene_std = np.nanstd(heatmap, axis=1, keepdims=True)
    gene_std[~np.isfinite(gene_std) | (gene_std == 0)] = 1
    heatmap_z = (heatmap - gene_mean) / gene_std
    return heatmap, heatmap_z, bin_centers, bin_edges, occupied[first_bin:last_bin + 1]


def plot_marker_heatmap(
    adata_obj,
    marker_groups,
    group_order,
    group_colors,
    pseudotime_col,
    cluster_col,
    cluster_value=None,
    title=None,
    output_name='marker_heatmap.png',
    mask=None,
    target_label=None,
    strip_col=None,
    strip_order=None,
    strip_colors=None,
    n_bins=120,
    smooth_sigma=2.5,
    trim_fraction=0.05,
    z_clip=2,
    expression_layer='lognorm',
    show_module_legend=True,
    show_group_labels=True,
    figsize=(16, 10),
    output_dir=None,
    project_dir=None,
):
    if pseudotime_col not in adata_obj.obs.columns:
        raise KeyError(f'Missing pseudotime column: {pseudotime_col}')
    if cluster_col is not None and cluster_col not in adata_obj.obs.columns:
        raise KeyError(f'Missing cluster label column: {cluster_col}')

    plot_mask = np.ones(adata_obj.n_obs, dtype=bool)
    if cluster_col is not None and cluster_value is not None:
        plot_mask &= adata_obj.obs[cluster_col].astype(str).eq(str(cluster_value)).to_numpy()
    if mask is not None:
        plot_mask &= np.asarray(mask, dtype=bool)
    plot_mask &= adata_obj.obs[pseudotime_col].notna().to_numpy()
    if plot_mask.sum() == 0:
        raise ValueError(f'No tubules available for {title or output_name}')

    plot_adata = adata_obj[plot_mask].copy()
    n_cells_before_trim = int(plot_adata.n_obs)
    if not 0 <= trim_fraction < 0.5:
        raise ValueError('trim_fraction must be in [0, 0.5)')
    n_trim_per_end = int(np.floor(n_cells_before_trim * trim_fraction))
    if n_trim_per_end > 0:
        pseudotime_order = np.argsort(plot_adata.obs[pseudotime_col].to_numpy(dtype=float), kind='stable')
        keep_idx = pseudotime_order[n_trim_per_end:n_cells_before_trim - n_trim_per_end]
        plot_adata = plot_adata[keep_idx].copy()
    if plot_adata.n_obs < 10:
        raise ValueError(f'Too few tubules remain after pseudospace trimming for {title or output_name}')

    available_groups, gene_order, gene_idx, missing = resolve_available_marker_groups(plot_adata, marker_groups, group_order)
    heatmap, heatmap_z, bin_centers, bin_edges, occupied_bins = binned_interpolated_expression(
        plot_adata, gene_idx, pseudotime_col, n_bins=n_bins, smooth_sigma=smooth_sigma, expression_layer=expression_layer,
    )
    heatmap_z = np.clip(heatmap_z, -z_clip, z_clip)

    module_names = [group for group in group_order if group in available_groups]
    gene_pos = {gene: i for i, gene in enumerate(gene_order)}
    module_scores = []
    for module_name in module_names:
        idx = [gene_pos[gene] for gene in available_groups[module_name]]
        module_scores.append(np.nanmean(heatmap_z[idx], axis=0))
    module_scores = np.vstack(module_scores)
    dominant_module = np.nanargmax(module_scores, axis=0)

    strip_values = dominant_module
    strip_cmap = ListedColormap([group_colors[name] for name in module_names])
    if strip_col is not None:
        if strip_col not in plot_adata.obs.columns:
            raise KeyError(f'Missing cluster strip column: {strip_col}')
        labels = plot_adata.obs[strip_col].astype(str).to_numpy()
        pos = plot_adata.obs[pseudotime_col].to_numpy(dtype=float)
        if strip_order is None:
            if hasattr(plot_adata.obs[strip_col], 'cat'):
                strip_order = [str(x) for x in plot_adata.obs[strip_col].cat.categories]
            else:
                strip_order = sorted(pd.unique(labels).astype(str))
        strip_order = [str(x) for x in strip_order if str(x) not in {'nan', 'None'}]
        strip_lookup = {label: i for i, label in enumerate(strip_order)}
        strip_values = np.full(len(bin_centers), np.nan, dtype=float)
        for bin_i in range(len(bin_centers)):
            if bin_i == len(bin_centers) - 1:
                in_bin = (pos >= bin_edges[bin_i]) & (pos <= bin_edges[bin_i + 1])
            else:
                in_bin = (pos >= bin_edges[bin_i]) & (pos < bin_edges[bin_i + 1])
            bin_labels = pd.Series(labels[in_bin]).replace({'nan': np.nan, 'None': np.nan}).dropna()
            if len(bin_labels):
                label = str(bin_labels.value_counts().idxmax())
                if label in strip_lookup:
                    strip_values[bin_i] = strip_lookup[label]
        if strip_colors is None:
            palette = plt.get_cmap('tab20').colors
            strip_colors = {label: palette[i % len(palette)] for i, label in enumerate(strip_order)}
        strip_cmap = ListedColormap([strip_colors.get(label, '#4D4D4D') for label in strip_order])

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 1, height_ratios=[1.2, 0.25, 8], hspace=0.08)
    ax_top = fig.add_subplot(gs[0])
    ax_strip = fig.add_subplot(gs[1], sharex=ax_top)
    ax_heat = fig.add_subplot(gs[2], sharex=ax_top)

    for module_idx, module_name in enumerate(module_names):
        ax_top.plot(bin_centers, module_scores[module_idx], color=group_colors[module_name], linewidth=2.5, label=module_name)
    ax_top.set_xlim(float(bin_edges[0]), float(bin_edges[-1]))
    ax_top.set_ylabel('Module score')
    ax_top.set_xticks([])
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)
    if show_module_legend:
        ax_top.legend(loc='upper right', frameon=False, ncol=min(3, max(1, len(module_names))), fontsize=9)

    ax_strip.imshow(
        strip_values[np.newaxis, :], aspect='auto', cmap=strip_cmap, interpolation='nearest',
        extent=[float(bin_edges[0]), float(bin_edges[-1]), 0, 1],
    )
    ax_strip.set_yticks([])
    ax_strip.set_xticks([])
    for spine in ax_strip.spines.values():
        spine.set_visible(False)

    im = ax_heat.imshow(
        heatmap_z, aspect='auto', cmap='bwr', vmin=-z_clip, vmax=z_clip, interpolation='nearest',
        extent=[float(bin_edges[0]), float(bin_edges[-1]), len(gene_order) - 0.5, -0.5],
    )
    ax_heat.set_xlim(float(bin_edges[0]), float(bin_edges[-1]))
    ax_heat.set_xlabel(pseudotime_col)
    xticks = np.linspace(float(bin_edges[0]), float(bin_edges[-1]), 6)
    ax_heat.set_xticks(xticks)
    ax_heat.set_xticklabels([f'{x:.2f}' for x in xticks])
    ax_heat.set_yticks(np.arange(len(gene_order)))
    ax_heat.set_yticklabels(gene_order, fontsize=10)
    ax_heat.yaxis.tick_right()
    ax_heat.tick_params(axis='y', length=0)

    if target_label:
        ax_heat.text(
            1.01, 1.01, target_label, transform=ax_heat.transAxes, ha='left', va='bottom', fontsize=10,
            bbox={'boxstyle': 'round,pad=0.3', 'facecolor': 'white', 'edgecolor': 'black', 'linewidth': 0.8},
        )

    offset = 0
    for module_name in module_names:
        genes = available_groups[module_name]
        center = offset + (len(genes) - 1) / 2
        if show_group_labels:
            ax_heat.text(
                -0.02, center, module_name, transform=ax_heat.get_yaxis_transform(), ha='right', va='center',
                fontsize=11, fontweight='bold', color=group_colors[module_name],
            )
        offset += len(genes)
        if offset < len(gene_order):
            ax_heat.axhline(offset - 0.5, color='black', linewidth=0.8)

    fig.subplots_adjust(left=0.10, right=0.76, top=0.92, bottom=0.08, hspace=0.08)
    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.03, pad=0.12)
    cbar.set_label('Per-gene z-score')
    fig.suptitle(title or output_name, fontsize=16, fontweight='bold', y=0.98)
    fig.savefig(output_dir / output_name, dpi=200, bbox_inches='tight')
    plt.show()

    return {
        'n_cells': int(plot_adata.n_obs),
        'n_cells_before_trim': n_cells_before_trim,
        'n_cells_trimmed_per_end': n_trim_per_end,
        'trim_fraction_per_end': float(trim_fraction),
        'cluster_col': cluster_col,
        'cluster_value': cluster_value,
        'pseudotime_col': pseudotime_col,
        'dpt_min': float(bin_edges[0]),
        'dpt_max': float(bin_edges[-1]),
        'genes_used': gene_order,
        'missing_markers': missing,
        'module_names': module_names,
        'output_path': (
            str((output_dir / output_name).relative_to(project_dir))
            if project_dir is not None
            else str(output_dir / output_name)
        ),
    }
