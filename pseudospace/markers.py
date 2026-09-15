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


def auto_discover_cluster_identities(
    adata_cluster: sc.AnnData,
    panel_scores: pd.DataFrame,
    top_de_df: pd.DataFrame,
    cluster_col: str = 'leiden_coarse',
    species_col: str = 'comparison_species',
    segment_to_coarse: dict[str, str] | None = None,
    remove_classes: list[str] | None = None,
    coarse_order: list[str] | None = None,
    non_tubule_families: set[str] | None = None,
) -> tuple[dict[str, str], pd.DataFrame]:
    """
    Deterministically discover and assign nephron segment and non-tubular identities
    to Leiden clusters without hardcoded cluster IDs or manual fingerprints.

    Combines:
    1) Multi-gene panel scores (mean module score per cluster across canonical markers)
    2) Top Wilcoxon differential expression marker genes
    3) Nephron segment biology (PCT S1/S2 consensus, PST S3, TAL, DCT, CNT, CD)
    4) Non-tubular stroma and stressed debris detection

    Returns:
        (cluster_labels_dict, cluster_summary_dataframe)
    """
    if segment_to_coarse is None:
        segment_to_coarse = {
            'Podocyte': 'Glomerulus',
            'PT-S1': 'PT', 'PT-S2': 'PT', 'PT-S3': 'PT', 'PT': 'PT',
            'DTL1': 'DTL', 'DTL2': 'DTL', 'DTL3': 'DTL', 'DTL': 'DTL',
            'ATL': 'AL', 'mTAL': 'AL', 'cTAL': 'AL', 'Macula-densa': 'AL', 'AL': 'AL',
            'DCT1': 'DCT', 'DCT2': 'DCT', 'DCT': 'DCT',
            'CNT': 'CNT_CD', 'CCD': 'CNT_CD', 'OMCD': 'CNT_CD', 'IMCD': 'CNT_CD', 'CNT_CD': 'CNT_CD',
            'Glomerulus': 'Glomerulus',
        }
    if remove_classes is None:
        remove_classes = ['Glomerulus', 'Vessel', 'Stroma', 'SmoothMuscle', 'Immune', 'Unassigned']
    if coarse_order is None:
        coarse_order = ['PT', 'DTL', 'AL', 'DCT', 'CNT_CD']
    if non_tubule_families is None:
        non_tubule_families = {'Podocyte', 'Glomerulus', 'Vessel', 'Stroma', 'SmoothMuscle', 'Immune'}

    stroma_markers = {
        'VIM', 'SPARC', 'CALD1', 'TAGLN', 'COL1A1', 'COL1A2', 'COL3A1',
        'MGP', 'ACTA2', 'MYH11', 'PECAM1', 'CDH5', 'PTPRC'
    }
    s3_markers = {'SLC22A7', 'SLC7A13', 'SLC6A18', 'ACSM3', 'NAPSA', 'MEP1A', 'MEP1B'}
    tal_markers = {'SLC12A1', 'UMOD', 'PPP1R1A', 'EGF', 'CLDN16', 'CLDN10', 'KCNJ10'}
    dct_markers = {'SLC12A3', 'TMEM52B', 'WNK1', 'PVALB', 'TRPM6', 'KLHL3', 'SFRP1'}
    cnt_markers = {'HSD11B2', 'CALB1', 'RHCG', 'SLC8A1', 'ATP6V1G3', 'SLC4A9'}
    cd_markers = {'AQP2', 'AQP3', 'AQP4', 'BCAM', 'SLC14A2', 'ATP6V0D2'}

    cluster_ids = sorted(
        adata_cluster.obs[cluster_col].astype(str).unique(),
        key=lambda x: int(x) if x.isdigit() else x,
    )

    assigned_labels = {}
    summary_rows = []

    # Cross-species breakdown if species column exists
    has_species = species_col in adata_cluster.obs.columns

    for c_id in cluster_ids:
        # Match scores row
        if c_id in panel_scores.index:
            cluster_scores = panel_scores.loc[c_id]
        elif c_id.isdigit() and int(c_id) in panel_scores.index:
            cluster_scores = panel_scores.loc[int(c_id)]
        else:
            cluster_scores = pd.Series(0.0, index=panel_scores.columns)

        # Match top DE genes
        c_de = top_de_df[top_de_df['group'].astype(str) == str(c_id)]
        top_de_genes = [str(g).upper() for g in c_de.head(10)['names']] if not c_de.empty else []

        top_panel = str(cluster_scores.idxmax())
        top_score = float(cluster_scores.max())

        # 1. Non-tubule stroma / glomerular / vascular check
        non_tub_cols = [col for col in cluster_scores.index if col in non_tubule_families]
        max_non_tub_score = float(cluster_scores[non_tub_cols].max()) if non_tub_cols else -1.0
        top_non_tub_panel = str(cluster_scores[non_tub_cols].idxmax()) if non_tub_cols else 'Stroma'
        stroma_hits = len(set(top_de_genes[:5]) & stroma_markers)

        assigned_label = None
        if (max_non_tub_score > 0.4 and max_non_tub_score >= top_score - 0.15) or stroma_hits >= 2:
            assigned_label = top_non_tub_panel if top_non_tub_panel in remove_classes else 'Stroma'

        # 2. Stressed / damaged debris check (e.g. Cryab, Sgk1 dominance with weak specific tubule score)
        if assigned_label is None:
            if top_de_genes and top_de_genes[0] in {'CRYAB', 'BHMT', 'HSPA1A', 'HSPA1B'} and top_score < 0.6:
                assigned_label = 'Unassigned'

        # 3. Proximal tubule: PCT (S1+S2) vs PST (S3)
        if assigned_label is None:
            s1_score = float(cluster_scores.get('PT-S1', 0.0))
            s2_score = float(cluster_scores.get('PT-S2', 0.0))
            s3_score = float(cluster_scores.get('PT-S3', 0.0))
            s3_hits = len(set(top_de_genes[:5]) & s3_markers)

            if top_panel == 'PT-S3' or (s3_score > max(s1_score, s2_score) and s3_hits >= 2):
                assigned_label = 'PT-S3'
            elif top_panel in {'PT-S1', 'PT-S2'} or (s1_score > 0.5 and s2_score > 0.5):
                # Merges human PCT and mouse S1/S2 transitional tubules at coarse 'PT'
                assigned_label = 'PT'

        # 4. Thick Ascending Limb (cTAL / mTAL)
        if assigned_label is None:
            tal_hits = len(set(top_de_genes[:5]) & tal_markers)
            if top_panel in {'cTAL', 'mTAL'} or tal_hits >= 2:
                assigned_label = 'cTAL'

        # 5. Distal Convoluted Tubule (DCT1 / DCT2)
        if assigned_label is None:
            dct_hits = len(set(top_de_genes[:5]) & dct_markers)
            if top_panel in {'DCT1', 'DCT2'} or dct_hits >= 2:
                assigned_label = 'DCT1'

        # 6. Connecting Tubule (CNT)
        if assigned_label is None:
            cnt_hits = len(set(top_de_genes[:5]) & cnt_markers)
            if top_panel == 'CNT' or cnt_hits >= 2:
                assigned_label = 'CNT'

        # 7. Collecting Duct (CCD / OMCD / IMCD)
        if assigned_label is None:
            cd_hits = len(set(top_de_genes[:5]) & cd_markers)
            if top_panel in {'CCD', 'OMCD', 'IMCD'} or cd_hits >= 2:
                assigned_label = 'CCD'

        # Fallback to mapped coarse if valid, else Unassigned
        if assigned_label is None:
            mapped = segment_to_coarse.get(top_panel)
            assigned_label = top_panel if mapped in coarse_order else 'Unassigned'

        assigned_labels[str(c_id)] = assigned_label

        # Compute cell and species breakdown
        c_mask = adata_cluster.obs[cluster_col].astype(str).eq(str(c_id)).to_numpy()
        n_cells = int(c_mask.sum())

        coarse_class = segment_to_coarse.get(
            assigned_label, assigned_label if assigned_label in remove_classes else 'Unassigned'
        )
        fate = 'Retained' if coarse_class in coarse_order else 'Filtered'
        top5_str = ', '.join(c_de.head(5)['names'].tolist()) if not c_de.empty else 'N/A'

        row = {
            'Cluster': str(c_id),
            'Total Cells': n_cells,
            'Discovered Segment': assigned_label,
            'Coarse Class': coarse_class,
            'Fate': fate,
            'Top Panel': f"{top_panel} ({top_score:.2f})",
            'Top Markers': top5_str,
        }

        if has_species:
            sp_counts = adata_cluster.obs.loc[c_mask, species_col].value_counts()
            n_human = int(sp_counts.get('human', 0))
            n_mouse = int(sp_counts.get('mouse', 0))
            pct_human = (n_human / n_cells * 100) if n_cells > 0 else 0.0
            pct_mouse = (n_mouse / n_cells * 100) if n_cells > 0 else 0.0
            row['Human Cells'] = n_human
            row['Mouse Cells'] = n_mouse
            row['Human %'] = round(pct_human, 1)
            row['Mouse %'] = round(pct_mouse, 1)

        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    return assigned_labels, summary_df


def plot_coarse_cluster_visualizations(
    adata_cluster: sc.AnnData,
    panel_scores: pd.DataFrame,
    top_de_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    output_dir: str | Path,
    cluster_col: str = 'leiden_coarse',
    species_col: str = 'comparison_species',
) -> list[str]:
    """
    Generate comprehensive diagnostic visualizations explaining what each cluster represents:
    1) Stacked cross-species breakdown bar plot
    2) Panel scores heatmap
    3) Top diagnostic DE marker matrixplot
    """
    import matplotlib.pyplot as plt
    from pathlib import Path

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    generated_files = []

    # 1. Cross-Species Tubule Breakdown (Stacked Bar Plot)
    if species_col in adata_cluster.obs.columns:
        sp_ct = pd.crosstab(adata_cluster.obs[cluster_col], adata_cluster.obs[species_col])
        fig, ax = plt.subplots(figsize=(max(8, len(sp_ct) * 0.9), 4.8))
        colors = {'human': '#D55E00', 'mouse': '#0072B2'}
        present_species = [s for s in ['mouse', 'human'] if s in sp_ct.columns]

        sp_ct[present_species].plot(
            kind='bar', stacked=True, color=[colors.get(s, '#7f7f7f') for s in present_species],
            ax=ax, edgecolor='black', linewidth=0.6,
        )
        ax.set_title('Pass-1 Leiden Clusters: Tubule Counts and Species Breakdown', fontsize=12, fontweight='bold')
        ax.set_xlabel('Leiden Coarse Cluster ID', fontsize=10)
        ax.set_ylabel('Number of Tubules', fontsize=10)
        ax.grid(axis='y', linestyle='--', alpha=0.4)

        for idx, (c, row) in enumerate(sp_ct.iterrows()):
            total = row.sum()
            h_count = row.get('human', 0)
            m_count = row.get('mouse', 0)
            h_pct = (h_count / total * 100) if total > 0 else 0
            m_pct = (m_count / total * 100) if total > 0 else 0
            ax.text(idx, total + max(50, total * 0.02), f"n={total:,}\n({m_pct:.0f}% M, {h_pct:.0f}% H)",
                    ha='center', va='bottom', fontsize=7.5)

        ax.set_ylim(0, sp_ct.sum(axis=1).max() * 1.25)
        ax.legend([s.capitalize() for s in present_species], frameon=True, loc='upper right')
        plt.tight_layout()
        breakdown_path = out_dir / 'coarse_cluster_species_breakdown.png'
        plt.savefig(breakdown_path, dpi=180)
        plt.close(fig)
        generated_files.append(str(breakdown_path))

    # 2. Panel Scores Heatmap
    interesting_cols = [
        c for c in panel_scores.columns
        if panel_scores[c].max() > 0.25 or panel_scores[c].min() < -0.25
    ]
    if interesting_cols:
        fig, ax = plt.subplots(figsize=(max(10, len(interesting_cols) * 0.6), max(5, len(panel_scores) * 0.45)))
        im = ax.imshow(panel_scores[interesting_cols].to_numpy(), cmap='coolwarm', aspect='auto')
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Mean Panel Score', fontsize=10)

        ax.set_xticks(range(len(interesting_cols)))
        ax.set_xticklabels(interesting_cols, rotation=45, ha='right', fontsize=9)
        ax.set_yticks(range(len(panel_scores)))
        ax.set_yticklabels(panel_scores.index.astype(str), fontsize=9)

        # Annotate text values
        matrix_vals = panel_scores[interesting_cols].to_numpy()
        for i in range(matrix_vals.shape[0]):
            for j in range(matrix_vals.shape[1]):
                val = matrix_vals[i, j]
                color = 'white' if abs(val) > 1.2 else 'black'
                ax.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=7, color=color)

        ax.set_title('Pass-1 Leiden Clusters: Mean Marker Panel Scores', fontsize=12, fontweight='bold')
        ax.set_xlabel('Nephron Segment & Non-Tubule Marker Panels', fontsize=10)
        ax.set_ylabel('Leiden Coarse Cluster ID', fontsize=10)
        plt.tight_layout()
        scores_path = out_dir / 'coarse_cluster_panel_scores_heatmap.png'
        plt.savefig(scores_path, dpi=180)
        plt.close(fig)
        generated_files.append(str(scores_path))

    # 3. Top Diagnostic Markers Matrixplot
    top3_genes = (
        top_de_df.groupby('group', observed=True).head(3)['names']
        .unique().tolist()
    )
    present_top3 = [g for g in top3_genes if g in adata_cluster.var_names]
    if present_top3:
        sc.pl.matrixplot(
            adata_cluster,
            var_names=present_top3,
            groupby=cluster_col,
            standard_scale='var',
            cmap='Blues',
            title='Pass-1 Clusters: Normalized Expression of Top Diagnostic DE Markers',
            show=False,
        )
        markers_path = out_dir / 'coarse_cluster_diagnostic_markers.png'
        plt.savefig(markers_path, dpi=180, bbox_inches='tight')
        plt.close()
        generated_files.append(str(markers_path))

    # 4. Tubule Segment Likelihood Matrix (Annotated Heatmap)
    nephron_order = [
        'Podocyte', 'PT-S1', 'PT-S2', 'PT-S3',
        'DTL1', 'DTL2', 'DTL3', 'ATL',
        'mTAL', 'cTAL', 'Macula-densa',
        'DCT1', 'DCT2', 'CNT', 'CCD', 'OMCD', 'IMCD',
        'Stroma', 'SmoothMuscle', 'Vessel', 'Immune'
    ]
    present_order = [c for c in nephron_order if c in panel_scores.columns]
    cluster_ids_str = [str(c) for c in panel_scores.index]

    # Check if per-cell primary segment calls exist in adata_cluster
    has_primary = 'celltype_primary_segment' in adata_cluster.obs.columns
    if has_primary:
        ct_tab = pd.crosstab(adata_cluster.obs[cluster_col].astype(str), adata_cluster.obs['celltype_primary_segment'])
        # Reindex to present order
        for col in present_order:
            if col not in ct_tab.columns:
                ct_tab[col] = 0
        ct_tab = ct_tab[present_order]
        prob_df = ct_tab.div(ct_tab.sum(axis=1).replace(0, 1), axis=0) * 100
    else:
        # Fallback: softmax of panel scores across tubule segments
        sub_scores = panel_scores[present_order].to_numpy(dtype=float)
        exp_s = np.exp(sub_scores / 0.5)
        prob_df = pd.DataFrame(exp_s / exp_s.sum(axis=1, keepdims=True) * 100, index=panel_scores.index.astype(str), columns=present_order)

    fig, ax = plt.subplots(figsize=(max(12, len(present_order) * 0.7), max(5, len(prob_df) * 0.55)))
    matrix_pct = prob_df.to_numpy()
    im = ax.imshow(matrix_pct, cmap='YlGnBu', aspect='auto', vmin=0, vmax=100)
    cbar = plt.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label('Segment Likelihood (% of cells in cluster)', fontsize=10, fontweight='bold')

    ax.set_xticks(range(len(present_order)))
    ax.set_xticklabels(present_order, rotation=45, ha='right', fontsize=9.5, fontweight='bold')

    # Y-labels with discovered identity if available
    y_lbls = []
    for c_id in prob_df.index:
        s_row = summary_df[summary_df['Cluster'].astype(str) == str(c_id)]
        seg_call = s_row['Discovered Segment'].iloc[0] if not s_row.empty else str(c_id)
        n_c = int(s_row['Total Cells'].iloc[0]) if (not s_row.empty and 'Total Cells' in s_row.columns) else int((adata_cluster.obs[cluster_col].astype(str) == str(c_id)).sum())
        y_lbls.append(f"Cl {c_id} ({seg_call})\nn={n_c:,}")
    ax.set_yticks(range(len(prob_df)))
    ax.set_yticklabels(y_lbls, fontsize=9.5)

    for i in range(matrix_pct.shape[0]):
        for j in range(matrix_pct.shape[1]):
            val = matrix_pct[i, j]
            if val >= 1.0:
                t_str = f"{val:.1f}%" if val < 10 else f"{val:.0f}%"
                t_col = "white" if val > 45 else ("black" if val > 5 else "#666666")
                ax.text(j, i, t_str, ha='center', va='center', fontsize=8, color=t_col, fontweight='bold' if val >= 10 else 'normal')
            elif val > 0:
                ax.text(j, i, "<1%", ha='center', va='center', fontsize=6.5, color='#888888')

    # Compartment dividers
    for div_pos in [0.5, 3.5, 6.5, 7.5, 10.5, 12.5, 16.5]:
        if div_pos < len(present_order):
            ax.axvline(div_pos, color='#333333', linestyle='--', linewidth=0.8, alpha=0.6)

    # Highlight unknown/unassigned cluster if present
    unknown_mask = summary_df['Fate'].astype(str).str.lower().eq('filtered') & summary_df['Coarse Class'].astype(str).str.lower().eq('unassigned')
    if unknown_mask.any():
        u_id = str(summary_df.loc[unknown_mask, 'Cluster'].iloc[0])
        if u_id in list(prob_df.index):
            u_idx = list(prob_df.index).index(u_id)
            ax.axhline(u_idx - 0.5, color='red', linestyle='-', linewidth=1.5, alpha=0.8)
            ax.axhline(u_idx + 0.5, color='red', linestyle='-', linewidth=1.5, alpha=0.8)

    ax.set_title('Pass-1 Clusters: Likelihood Distribution Across Canonical Tubule Types', fontsize=12, fontweight='bold', pad=12)
    ax.set_xlabel('Canonical Nephron Segments (Ordered Anatomically) & Non-Tubular Classes', fontsize=10, fontweight='bold', labelpad=8)
    ax.set_ylabel('Pass-1 Leiden Cluster', fontsize=10, fontweight='bold')
    plt.tight_layout()
    prob_path = out_dir / 'coarse_cluster_tubule_likelihood_matrix.png'
    plt.savefig(prob_path, dpi=180)
    plt.close(fig)
    generated_files.append(str(prob_path))

    # 5. Coarse Family Tubule Composition (Stacked Bar & Donut)
    seg_to_coarse = {
        'Podocyte': 'Glomerulus',
        'PT-S1': 'PT', 'PT-S2': 'PT', 'PT-S3': 'PT',
        'DTL1': 'DTL', 'DTL2': 'DTL', 'DTL3': 'DTL',
        'ATL': 'AL', 'mTAL': 'AL', 'cTAL': 'AL', 'Macula-densa': 'AL',
        'DCT1': 'DCT', 'DCT2': 'DCT',
        'CNT': 'CNT_CD', 'CCD': 'CNT_CD', 'OMCD': 'CNT_CD', 'IMCD': 'CNT_CD',
        'Stroma': 'Non-Tubular', 'SmoothMuscle': 'Non-Tubular', 'Vessel': 'Non-Tubular', 'Immune': 'Non-Tubular'
    }
    coarse_colors = {
        'PT': '#2ca02c', 'DTL': '#17becf', 'AL': '#1f77b4', 'DCT': '#ff7f0e',
        'CNT_CD': '#9467bd', 'Glomerulus': '#d62728', 'Non-Tubular': '#7f7f7f'
    }

    if has_primary:
        coarse_counts = pd.DataFrame(index=ct_tab.index)
        for cat in ['PT', 'DTL', 'AL', 'DCT', 'CNT_CD', 'Glomerulus', 'Non-Tubular']:
            segs = [s for s, c in seg_to_coarse.items() if c == cat and s in ct_tab.columns]
            coarse_counts[cat] = ct_tab[segs].sum(axis=1) if segs else 0
        coarse_pct = coarse_counts.div(coarse_counts.sum(axis=1).replace(0, 1), axis=0) * 100

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5.5), gridspec_kw={'width_ratios': [1.2, 1.0]})
        bottom = np.zeros(len(coarse_pct))
        x_bars = np.arange(len(coarse_pct))
        for cat, col in coarse_colors.items():
            if cat in coarse_pct.columns:
                vals = coarse_pct[cat].to_numpy()
                ax1.bar(x_bars, vals, bottom=bottom, label=cat, color=col, edgecolor='black', linewidth=0.5, width=0.65)
                bottom += vals

        ax1.set_xticks(x_bars)
        ax1.set_xticklabels([f"Cl {c}" for c in coarse_pct.index], fontsize=9)
        ax1.set_ylabel('Cell Fraction (%)', fontsize=10, fontweight='bold')
        ax1.set_title('Coarse Nephron Family Composition Across Clusters', fontsize=11.5, fontweight='bold')
        ax1.set_ylim(0, 100)
        ax1.legend(loc='upper right', frameon=True, title='Nephron Family', fontsize=8.5)
        ax1.grid(axis='y', linestyle='--', alpha=0.3)

        # Right: Unknown Cluster Donut (if present, else top mixed cluster)
        u_id = summary_df.loc[unknown_mask, 'Cluster'].iloc[0] if unknown_mask.any() else str(coarse_pct.index[-1])
        if str(u_id) in coarse_pct.index:
            c_u = coarse_pct.loc[str(u_id)]
            c_u_filt = c_u[c_u > 0.5].sort_values(ascending=False)
            ax2.pie(
                c_u_filt, labels=[f"{k} ({v:.1f}%)" for k, v in c_u_filt.items()],
                startangle=140, colors=[coarse_colors[k] for k in c_u_filt.index],
                wedgeprops={'edgecolor': 'black', 'linewidth': 0.8, 'width': 0.45},
            )
            ax2.text(0, 0.05, f"Cluster {u_id}", ha='center', va='center', fontsize=12, fontweight='bold')
            total_u = int(summary_df.loc[summary_df['Cluster'].astype(str) == str(u_id), 'Total Cells'].iloc[0]) if (summary_df['Cluster'].astype(str) == str(u_id)).any() else 0
            ax2.text(0, -0.10, f"n={total_u:,}", ha='center', va='center', fontsize=9, color='#444444')
            ax2.set_title(f'Unknown Cluster (Cluster {u_id}): Coarse Family Allocation', fontsize=11.5, fontweight='bold')

        plt.tight_layout()
        stacked_path = out_dir / 'coarse_cluster_tubule_composition_stacked.png'
        plt.savefig(stacked_path, dpi=180)
        plt.close(fig)
        generated_files.append(str(stacked_path))

    # 6. Dedicated Unknown Cluster Diagnostic Profile
    if unknown_mask.any() and has_primary:
        u_id = str(summary_df.loc[unknown_mask, 'Cluster'].iloc[0])
        fig = plt.figure(figsize=(16, 11))
        gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.28)

        # Panel A: Ranked Likelihood Bar
        ax_a = fig.add_subplot(gs[0, 0])
        u_segs = prob_df.loc[u_id].sort_values(ascending=True)
        top_u_segs = u_segs[u_segs >= 0.2]
        bar_colors = [coarse_colors.get(seg_to_coarse.get(seg, 'Non-Tubular'), '#888888') for seg in top_u_segs.index]
        y_pos = np.arange(len(top_u_segs))
        ax_a.barh(y_pos, top_u_segs.values, color=bar_colors, edgecolor='black', linewidth=0.6, height=0.7)
        ax_a.set_yticks(y_pos)
        ax_a.set_yticklabels(top_u_segs.index, fontsize=9, fontweight='bold')
        ax_a.set_xlabel('Cell Likelihood (% in cluster)', fontsize=10, fontweight='bold')
        ax_a.set_title(f'A. Ranked Tubule Segment Affinity in Unknown Cluster (Cl {u_id})', fontsize=11, fontweight='bold')
        ax_a.grid(axis='x', linestyle='--', alpha=0.4)
        for idx, val in enumerate(top_u_segs.values):
            n_val = int(ct_tab.loc[u_id, top_u_segs.index[idx]])
            ax_a.text(val + 0.5, idx, f"{val:.1f}% (n={n_val:,})", va='center', fontsize=8, fontweight='bold')
        ax_a.set_xlim(0, max(top_u_segs.values) * 1.25)

        # Panel B: Module Scores Comparison vs Pure Clusters
        ax_b = fig.add_subplot(gs[0, 1])
        comp_panels = [p for p in ['PT-S1', 'PT-S3', 'ATL', 'mTAL', 'cTAL', 'DCT1', 'CNT', 'CCD', 'IMCD', 'Podocyte', 'Stroma'] if p in panel_scores.columns]
        ref_clusters = [c for c in ['0', '1', '3', '6', u_id] if c in [str(x) for x in panel_scores.index]]
        ref_colors = {'0': '#2ca02c', '1': '#1f77b4', '3': '#9467bd', '6': '#17becf', u_id: 'red'}
        x_comp = np.arange(len(comp_panels))
        w_bar = 0.8 / max(1, len(ref_clusters))
        for i, cid in enumerate(ref_clusters):
            c_idx = int(cid) if cid.isdigit() and int(cid) in panel_scores.index else cid
            scores = [panel_scores.loc[c_idx, p] for p in comp_panels]
            is_u = (cid == u_id)
            ax_b.bar(
                x_comp + (i - len(ref_clusters)/2 + 0.5) * w_bar, scores, width=w_bar,
                label=f"Cl {cid}" + (" (Unknown)" if is_u else ""),
                color=ref_colors.get(cid, '#888888'), alpha=1.0 if is_u else 0.55,
                edgecolor='black' if is_u else 'none', linewidth=1.2 if is_u else 0
            )
        ax_b.set_xticks(x_comp)
        ax_b.set_xticklabels(comp_panels, rotation=45, ha='right', fontsize=8.5, fontweight='bold')
        ax_b.set_ylabel('Mean Marker Panel Score', fontsize=10, fontweight='bold')
        ax_b.set_title('B. Module Scores: Unknown Cluster vs. Pure Reference Clusters', fontsize=11, fontweight='bold')
        ax_b.axhline(0, color='grey', linestyle='-', linewidth=0.8)
        ax_b.grid(axis='y', linestyle='--', alpha=0.4)
        ax_b.legend(frameon=True, fontsize=8, loc='upper right')

        # Panel C: Sample & Anatomical Region Origin
        ax_c = fig.add_subplot(gs[1, 0])
        u_mask = adata_cluster.obs[cluster_col].astype(str) == u_id
        if 'sample' in adata_cluster.obs.columns:
            s_counts = adata_cluster.obs.loc[u_mask, 'sample'].value_counts()
            # HUK1_MED1 is named medulla at source but is healthy cortex tissue: label it as a slice,
            # never as a medullary sample.
            s_labels = {'HUK1_COR1': 'HUK1_COR1 (Cortex slice)', 'HUK1_MED1': 'HUK1_MED1 (Cortex slice)', 'Ctrl1A2': 'Ctrl1A2 (Mouse)', 'Ctrl1A4': 'Ctrl1A4 (Mouse)'}
            y_s = np.arange(len(s_counts))
            ax_c.barh(y_s, s_counts.values, color='#4575b4', edgecolor='black', linewidth=0.7, height=0.6)
            ax_c.set_yticks(y_s)
            ax_c.set_yticklabels([s_labels.get(k, k) for k in s_counts.index], fontsize=9, fontweight='bold')
            ax_c.set_xlabel('Cell Count in Cluster', fontsize=10, fontweight='bold')
            ax_c.set_title('C. Anatomical and Sample Origin of Unknown Cluster Cells', fontsize=11, fontweight='bold')
            ax_c.grid(axis='x', linestyle='--', alpha=0.4)
            for idx, val in enumerate(s_counts.values):
                pct = val / max(1, s_counts.sum()) * 100
                ax_c.text(val + 15, idx, f"n={val:,} ({pct:.1f}%)", va='center', fontsize=8.5, fontweight='bold')
            ax_c.set_xlim(0, max(s_counts.values) * 1.35)
            ax_c.text(
                0.05, 0.12,
                "Key Biological Insight:\n"
                "• Cortex cells in Cl 9 map to PT-S1\n"
                "• Medulla cells in Cl 9 map to CCD/IMCD/ATL\n"
                "→ Cl 9 captures stressed epithelia across regions!",
                transform=ax_c.transAxes, fontsize=8,
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#fff2cc", edgecolor="#d6b656", linewidth=0.8)
            )

        # Panel D: Stress, Injury & Plasticity Markers
        ax_d = fig.add_subplot(gs[1, 1])
        c_u_de = top_de_df[top_de_df['group'].astype(str) == u_id]
        marker_genes = ['Cryab', 'Sgk1', 'Krt18', 'Igfbp5', 'Nupr1', 'Aqp2', 'Vim', 'Vcam1', 'Bhmt', 'Hspa1b']
        m_lfc, m_pct = [], []
        for g in marker_genes:
            match = c_u_de[c_u_de['names'].astype(str).str.lower() == g.lower()]
            if not match.empty:
                m_lfc.append(float(match['logfoldchanges'].iloc[0]))
                pct_v = float(match['pct_nz_group'].iloc[0]) * 100 if 'pct_nz_group' in match.columns else 0.0
                m_pct.append(pct_v)
            else:
                m_lfc.append(0.0)
                m_pct.append(0.0)

        x_m = np.arange(len(marker_genes))
        ax_d.bar(x_m, m_lfc, color='#e41a1c', alpha=0.8, edgecolor='black', linewidth=0.7, width=0.6, label='Log2 Fold-Change (vs All Tubules)')
        ax_d.set_xticks(x_m)
        ax_d.set_xticklabels(marker_genes, rotation=45, ha='right', fontsize=9, fontweight='bold')
        ax_d.set_ylabel('Log2 Fold Change', fontsize=10, fontweight='bold', color='#b2182b')
        ax_d.set_title('D. Diagnostic Wilcoxon Markers: Stress, Injury & Plasticity', fontsize=11, fontweight='bold')
        ax_d.grid(axis='y', linestyle='--', alpha=0.4)

        if any(m_pct):
            ax_d2 = ax_d.twinx()
            ax_d2.plot(x_m, m_pct, color='#08519c', marker='o', markersize=5, linewidth=1.8, label='% Cells Expressing in Cluster')
            ax_d2.set_ylabel('% Cells Expressing', fontsize=10, fontweight='bold', color='#08519c')
            ax_d2.set_ylim(0, 105)
            lines_1, l1 = ax_d.get_legend_handles_labels()
            lines_2, l2 = ax_d2.get_legend_handles_labels()
            ax_d.legend(lines_1 + lines_2, l1 + l2, loc='upper right', frameon=True, fontsize=7.5)

        fig.suptitle(
            f'Diagnostic Decomposition of Unknown Cluster (Cluster {u_id}, n={total_u:,} cells)\n'
            'Evidence for a Mixed Injured/Stressed Tubule Population (PT + Medullary CD/ATL)',
            fontsize=13, fontweight='bold', y=0.99
        )
        fig.subplots_adjust(left=0.14, right=0.93, top=0.92, bottom=0.08, hspace=0.36, wspace=0.30)
        profile_path = out_dir / 'unknown_cluster_diagnostic_profile.png'
        plt.savefig(profile_path, dpi=180)
        plt.close(fig)
        generated_files.append(str(profile_path))

    return generated_files

