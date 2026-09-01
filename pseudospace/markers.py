"""Gene-symbol lookup / alias resolution, the total marker axis, and DE-based
cluster annotation.

Extracted verbatim from ``6_mouse_only_pseudospace.ipynb`` (Sections 2-3). Marker
panels themselves stay in the notebook and are passed in as arguments.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import spearmanr


def build_gene_lookup(adata):
    var_names = pd.Index(adata.var_names.astype(str))
    lookup = {name.upper(): name for name in var_names}
    for alias_col in ['gene_symbol', 'gene_symbols', 'symbol', 'gene_name', 'feature_name', 'features']:
        if alias_col in adata.var.columns:
            for idx_name, alias in zip(var_names, adata.var[alias_col].astype(str)):
                alias_upper = alias.upper()
                if alias_upper and alias_upper != 'NAN' and alias_upper not in lookup:
                    lookup[alias_upper] = idx_name
    return lookup


def resolve_gene_aliases(adata_obj, marker_groups):
    var_names = pd.Index(adata_obj.var_names.astype(str))
    lookup = {name.upper(): name for name in var_names}
    alias_columns = [
        col for col in ['gene_symbol', 'gene_symbols', 'symbol', 'gene_name', 'feature_name', 'features']
        if col in adata_obj.var.columns
    ]
    for col in alias_columns:
        for idx_name, alias in zip(var_names, adata_obj.var[col].astype(str)):
            alias_upper = alias.upper()
            if alias_upper and alias_upper != 'NAN' and alias_upper not in lookup:
                lookup[alias_upper] = idx_name
    resolved = {}
    missing = []
    for group_name, markers in marker_groups.items():
        for marker in markers:
            match = None
            for candidate in marker['candidates']:
                match = lookup.get(candidate.upper())
                if match is not None:
                    break
            if match is None:
                missing.append(marker['label'])
            else:
                resolved[marker['label']] = match
    return resolved, missing


def resolve_marker_list(adata_obj, markers):
    resolved, missing = resolve_gene_aliases(adata_obj, {'all': markers})
    ordered = [marker['label'] for marker in markers if marker['label'] in resolved]
    return ordered, resolved, missing


def marker_labels_for_group(marker_groups, group_order):
    labels = []
    for group in group_order:
        for marker in marker_groups.get(group, []):
            labels.append(marker['label'])
    return labels


def resolve_available_marker_groups(adata_obj, marker_groups, group_order):
    resolved_map, missing = resolve_gene_aliases(adata_obj, marker_groups)
    available_groups = {}
    for group in group_order:
        keep = [marker['label'] for marker in marker_groups.get(group, []) if marker['label'] in resolved_map]
        if keep:
            available_groups[group] = keep
    if not available_groups:
        preview_cols = [
            col for col in ['gene_symbol', 'gene_symbols', 'symbol', 'gene_name', 'feature_name', 'features']
            if col in adata_obj.var.columns
        ]
        raise ValueError(
            'None of the requested marker genes were found. '
            f'Checked adata.var_names and alias columns: {preview_cols if preview_cols else "none"}.'
        )
    gene_order = [gene for group in group_order if group in available_groups for gene in available_groups[group]]
    gene_idx = [adata_obj.var_names.get_loc(resolved_map[gene]) for gene in gene_order]
    return available_groups, gene_order, gene_idx, missing


def compute_total_marker_axis(
    adata,
    marker_groups,
    output_col="total_marker_axis",
    layer="lognorm",
    min_markers_per_group=2,
):
    """
    Compute a PT-like to CD-like positional axis.

    Each marker is z-scored across tubules before averaging, preventing
    highly expressed genes from dominating the module score.

    axis = mean_z(late markers) - mean_z(early markers)
    """

    gene_lookup = {str(gene).upper(): gene for gene in adata.var_names}

    resolved = {}
    missing = {}

    for group, requested_genes in marker_groups.items():
        resolved[group] = [
            gene_lookup[gene.upper()]
            for gene in requested_genes
            if gene.upper() in gene_lookup
        ]
        missing[group] = [
            gene
            for gene in requested_genes
            if gene.upper() not in gene_lookup
        ]

        print(
            f"{group}: resolved {len(resolved[group])}/{len(requested_genes)} "
            f"markers: {resolved[group]}"
        )
        if missing[group]:
            print(f"  Missing: {missing[group]}")

        if len(resolved[group]) < min_markers_per_group:
            raise ValueError(
                f"{group!r} resolved only {len(resolved[group])} markers; "
                f"at least {min_markers_per_group} are required."
            )

    all_genes = resolved["early"] + resolved["late"]

    if layer in adata.layers:
        expression = adata[:, all_genes].layers[layer]
        expression_source = f"layers[{layer!r}]"
    else:
        expression = adata[:, all_genes].X
        expression_source = "X"
        print(
            f"WARNING: layer {layer!r} not found; using adata.X to compute "
            "the marker axis."
        )

    if sparse.issparse(expression):
        expression = expression.toarray()

    # Shape: tubules × marker genes.
    expression = np.asarray(expression, dtype=np.float64)

    # Z-score each gene across tubules.
    gene_mean = np.nanmean(expression, axis=0, keepdims=True)
    gene_std = np.nanstd(expression, axis=0, keepdims=True)
    gene_std[~np.isfinite(gene_std) | (gene_std < 1e-8)] = 1.0

    marker_z = (expression - gene_mean) / gene_std
    marker_z = np.nan_to_num(marker_z, nan=0.0, posinf=0.0, neginf=0.0)

    n_early = len(resolved["early"])

    early_score = np.nanmean(marker_z[:, :n_early], axis=1)
    late_score = np.nanmean(marker_z[:, n_early:], axis=1)
    marker_axis = late_score - early_score

    if not np.isfinite(marker_axis).all():
        raise ValueError("The computed total marker axis contains nonfinite values.")

    adata.obs["total_early_marker_score"] = early_score
    adata.obs["total_late_marker_score"] = late_score
    adata.obs[output_col] = marker_axis

    print()
    print("Total marker axis computed:")
    print(f"  Expression source: {expression_source}")
    print(f"  Output column: adata.obs[{output_col!r}]")
    print(f"  Mean: {marker_axis.mean():.4f}")
    print(f"  SD: {marker_axis.std():.4f}")
    print(f"  Range: {marker_axis.min():.4f} to {marker_axis.max():.4f}")
    print("  Direction: low = PT-like; high = CD-like")

    # Diagnostic only: this does not define the axis or determine its direction.
    if "total_scanpy_dpt" in adata.obs:
        dpt = adata.obs["total_scanpy_dpt"].to_numpy(dtype=float)
        finite = np.isfinite(dpt) & np.isfinite(marker_axis)
        rho = spearmanr(dpt[finite], marker_axis[finite]).correlation
        print(f"  Spearman(marker axis, existing DPT): {rho:+.4f}")

    return {
        "resolved_markers": resolved,
        "missing_markers": missing,
        "early_score": early_score,
        "late_score": late_score,
        "marker_axis": marker_axis,
    }


def available_markers_by_group(adata, marker_dict, gene_lookup=None):
    """Resolve a ``{group: [symbols]}`` panel to the genes actually present in
    ``adata`` (case-insensitive), preserving group and gene order and dropping
    empty groups. The returned ordered dict feeds ``sc.pl.dotplot(var_names=...)``
    directly (scanpy draws group brackets from a dict) for classic marker-based
    cluster annotation. Missing genes are returned separately for reporting.
    """
    if gene_lookup is None:
        gene_lookup = build_gene_lookup(adata)
    resolved, missing = {}, {}
    for group, genes in marker_dict.items():
        present = [gene_lookup[g.upper()] for g in genes if g.upper() in gene_lookup]
        absent = [g for g in genes if g.upper() not in gene_lookup]
        if present:
            resolved[group] = present
        if absent:
            missing[group] = absent
    return resolved, missing


def assign_cluster_labels(adata, cluster_col, mapping, out_col, unknown='Unknown'):
    """Apply an explicit ``{cluster_id: label}`` dict to per-cell labels (classic
    manual annotation). Cluster ids are matched as strings; any cluster absent
    from ``mapping`` is labeled ``unknown`` and reported. Writes ``out_col`` and
    returns the per-cell label Series.
    """
    cluster_str = adata.obs[cluster_col].astype(str)
    str_mapping = {str(k): v for k, v in mapping.items()}
    present_clusters = set(cluster_str.unique())
    unmapped = sorted(present_clusters - set(str_mapping))
    if unmapped:
        print(f'assign_cluster_labels: {len(unmapped)} cluster(s) not in mapping -> '
              f"'{unknown}': {unmapped}")
    labels = cluster_str.map(str_mapping).fillna(unknown)
    adata.obs[out_col] = labels.to_numpy()
    print(f"{out_col} value counts:\n{labels.value_counts()}")
    return labels


def annotate_clusters_by_de(adata, cluster_col, marker_genes_by_family, margin_threshold, gene_lookup):
    """
    DE-primary cluster annotation.

    Primary: Wilcoxon rank_genes_groups score for each family's marker panel, averaged
    per cluster. Cross-check: mean broad marker score per cluster (already computed as
    obs[f'{family}_broad_score']).

    Family DE scores are min-max normalized per cluster (across families) before computing
    the assignment margin, so `margin_threshold` is comparable across resolutions and
    clusters regardless of the raw Wilcoxon z-score scale.

    Returns (cell_labels, cluster_table) where cluster_table has one row per cluster with
    columns: leiden_cluster, n_cells, top_markers, mean_score_<family>..., de_score_<family>...,
    assigned_label, assignment_margin.
    """
    rank_key = f'rank_genes_{cluster_col}'
    sc.tl.rank_genes_groups(adata, groupby=cluster_col, method='wilcoxon', key_added=rank_key)

    families = list(marker_genes_by_family)
    family_gene_map = {
        family: [gene_lookup[g.upper()] for g in genes if g.upper() in gene_lookup]
        for family, genes in marker_genes_by_family.items()
    }

    if hasattr(adata.obs[cluster_col], 'cat'):
        clusters = adata.obs[cluster_col].cat.categories
    else:
        clusters = sorted(adata.obs[cluster_col].astype(str).unique())

    rows = []
    for cluster in clusters:
        de_df = sc.get.rank_genes_groups_df(adata, group=str(cluster), key=rank_key)
        score_lookup = dict(zip(de_df['names'], de_df['scores']))

        de_scores = {}
        for family, genes in family_gene_map.items():
            vals = [score_lookup.get(g, 0.0) for g in genes]
            de_scores[family] = float(np.mean(vals)) if vals else float('-inf')

        cell_mask = adata.obs[cluster_col].astype(str).eq(str(cluster)).to_numpy()
        mean_scores = {
            family: float(adata.obs.loc[cell_mask, f'{family}_broad_score'].mean())
            for family in families
            if f'{family}_broad_score' in adata.obs.columns
        }

        top_markers = de_df.sort_values('scores', ascending=False)['names'].head(5).tolist()

        row = {
            'leiden_cluster': str(cluster),
            'n_cells': int(cell_mask.sum()),
            'top_markers': '; '.join(top_markers),
        }
        for family in families:
            row[f'de_score_{family}'] = de_scores.get(family, float('-inf'))
            row[f'mean_score_{family}'] = mean_scores.get(family, np.nan)
        rows.append(row)

    cluster_table = pd.DataFrame(rows).set_index('leiden_cluster')

    de_cols = [f'de_score_{family}' for family in families]
    de_matrix = cluster_table[de_cols]
    row_min = de_matrix.min(axis=1)
    row_max = de_matrix.max(axis=1)
    row_range = (row_max - row_min).replace(0, 1.0)
    norm = de_matrix.sub(row_min, axis=0).div(row_range, axis=0)
    norm.columns = families

    ranked = norm.to_numpy().argsort(axis=1)[:, ::-1]
    top_family_idx = ranked[:, 0]
    second_family_idx = ranked[:, 1] if norm.shape[1] > 1 else ranked[:, 0]
    norm_values = norm.to_numpy()
    top_val = norm_values[np.arange(len(norm)), top_family_idx]
    second_val = norm_values[np.arange(len(norm)), second_family_idx]
    margin = top_val - second_val

    assigned = np.array(families)[top_family_idx]
    assigned = np.where(margin < margin_threshold, 'Transitional', assigned)

    cluster_table['assigned_label'] = assigned
    cluster_table['assignment_margin'] = margin
    cluster_table = cluster_table.reset_index()

    label_map = dict(zip(cluster_table['leiden_cluster'], cluster_table['assigned_label']))
    cell_labels = adata.obs[cluster_col].astype(str).map(label_map)

    return cell_labels, cluster_table
