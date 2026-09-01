"""Diffusion-pseudotime helpers: representation/root/component selection, subset
DPT, and saving the canonical pseudotime AnnData.

Extracted verbatim from ``6_mouse_only_pseudospace.ipynb`` (Section 2-3).
Functions that previously read notebook globals now take them as explicit
arguments (defaults preserve the notebook's values).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

from .markers import resolve_available_marker_groups
from .heatmaps import matrix_to_dense, robust_scale_rows


def orient_and_normalize(values, marker_axis, min_valid=100):
    """Flip sign so values increase with marker axis, then clip to p5-p95 and rescale to [0, 1]."""
    values = np.asarray(values, dtype=float)
    marker_axis = np.asarray(marker_axis, dtype=float)
    valid = np.isfinite(values) & np.isfinite(marker_axis)
    oriented = values.copy()
    if valid.sum() >= min_valid:
        corr = spearmanr(oriented[valid], marker_axis[valid]).correlation
        if np.isfinite(corr) and corr < 0:
            oriented = -oriented
    finite = np.isfinite(oriented)
    out = np.full(oriented.shape, np.nan, dtype=float)
    if finite.any():
        p5, p95 = np.nanpercentile(oriented[finite], [5, 95])
        out[finite] = np.clip((oriented[finite] - p5) / max(p95 - p5, 1e-8), 0, 1)
    return out


def choose_representation(adata, prefix, n_pcs_fallback=30, random_state=0):
    preferred_reps = [f'X_pca_harmony_{prefix}', 'X_pca_harmony', 'X_harmony', 'X_pca']
    for rep in preferred_reps:
        if rep in adata.obsm:
            return rep
    n_comps = min(n_pcs_fallback, adata.n_obs - 1, adata.n_vars - 1)
    sc.pp.pca(adata, n_comps=n_comps, random_state=random_state)
    return 'X_pca'


def choose_root_global(marker_axis, adata, rep_key, bottom_quantile=0.01):
    """Fallback root: medoid of the bottom-quantile marker-axis cells (notebook 2's original choice)."""
    finite_mask = np.isfinite(marker_axis)
    if not finite_mask.any():
        raise ValueError('marker_axis has no finite values; cannot choose a root.')
    threshold = np.nanpercentile(marker_axis[finite_mask], 100 * bottom_quantile)
    candidates = np.where(finite_mask & (marker_axis <= threshold))[0]
    if candidates.size == 1:
        return int(candidates[0])
    coords = adata.obsm[rep_key][candidates]
    centroid = coords.mean(axis=0, keepdims=True)
    dists = np.linalg.norm(coords - centroid, axis=1)
    return int(candidates[int(np.argmin(dists))])


def choose_root_pt_cluster(adata, marker_axis, leiden_col, rep_key, family_col='broad_tubule_marker_call', pt_label='PT'):
    """
    Root at the medoid of the lowest-marker-axis PT-labelled Leiden cluster.
    Returns (None, None) if no clean PT cluster exists.
    """
    pt_mask = adata.obs[family_col].astype(str).eq(pt_label).to_numpy()
    if not pt_mask.any():
        return None, None

    clusters_for_pt = adata.obs.loc[pt_mask, leiden_col].astype(str).to_numpy()
    axis_for_pt = marker_axis[pt_mask]
    cluster_means = pd.Series(axis_for_pt, index=clusters_for_pt).groupby(level=0).mean()
    if cluster_means.empty or cluster_means.isna().all():
        return None, None

    best_cluster = cluster_means.idxmin()
    best_cluster_mask = pt_mask & (adata.obs[leiden_col].astype(str) == best_cluster).to_numpy()
    candidate_positions = np.flatnonzero(best_cluster_mask)
    if candidate_positions.size == 0:
        return None, None

    coords = adata.obsm[rep_key][candidate_positions]
    centroid = coords.mean(axis=0, keepdims=True)
    dists = np.linalg.norm(coords - centroid, axis=1)
    root_idx = int(candidate_positions[int(np.argmin(dists))])
    return root_idx, best_cluster


def choose_diffusion_component(adata, marker_axis, n_components_to_test, eigenvalue_floor=0.5, min_valid=100):
    diffmap = adata.obsm['X_diffmap']
    eigvals = np.asarray(adata.uns['diffmap_evals'], dtype=float)
    n_components = min(n_components_to_test, diffmap.shape[1], eigvals.size)

    rows = []
    for i in range(n_components):
        values = diffmap[:, i].astype(float)
        valid = np.isfinite(values) & np.isfinite(marker_axis)
        corr = (
            spearmanr(values[valid], marker_axis[valid]).correlation
            if valid.sum() >= min_valid
            else np.nan
        )
        rows.append({
            'component': i,
            'eigenvalue': float(eigvals[i]),
            'spearman_with_marker_axis': corr,
            'score': (abs(corr) * eigvals[i]) if (np.isfinite(corr) and eigvals[i] >= eigenvalue_floor) else np.nan,
        })
    component_df = pd.DataFrame(rows)

    if component_df['score'].notna().any():
        best_component = int(component_df.loc[component_df['score'].idxmax(), 'component'])
    else:
        fallback = component_df['spearman_with_marker_axis'].abs().fillna(-1)
        best_component = int(component_df.loc[fallback.idxmax(), 'component'])

    return best_component, component_df


def recompute_subset_dpt(adata_obj, segment, marker_groups, group_order, output_col, n_neighbors=30, random_state=0):
    subset_mask = adata_obj.obs['broad_tubule_marker_call'].astype(str).eq(segment).to_numpy()
    if subset_mask.sum() < 30:
        raise ValueError(f'{segment}: too few tubules for subset-specific DPT')

    subset = adata_obj[subset_mask].copy()
    representation = next((key for key in ['X_harmony', 'X_pca'] if key in subset.obsm), None)
    if representation is None:
        raise KeyError(f'{segment}: neither X_harmony nor X_pca is available')

    available_groups, gene_order, gene_idx, missing = resolve_available_marker_groups(subset, marker_groups, group_order)
    available_order = [group for group in group_order if group in available_groups]
    if len(available_order) < 2:
        raise ValueError(f'{segment}: at least two ordered marker modules are required')

    expression = subset.layers['lognorm'] if 'lognorm' in subset.layers else subset.X
    marker_expression = matrix_to_dense(expression[:, gene_idx]).T
    marker_z = robust_scale_rows(marker_expression)
    gene_position = {label: i for i, label in enumerate(gene_order)}

    early_group = available_order[0]
    late_group = available_order[-1]
    early_idx = [gene_position[label] for label in available_groups[early_group]]
    late_idx = [gene_position[label] for label in available_groups[late_group]]
    early_score = np.nanmean(marker_z[early_idx], axis=0)
    late_score = np.nanmean(marker_z[late_idx], axis=0)
    marker_axis = late_score - early_score

    neighbors_key = f'{segment.lower()}_subset_dpt_neighbors'
    n_neighbors_actual = min(int(n_neighbors), subset.n_obs - 1)
    sc.pp.neighbors(subset, n_neighbors=n_neighbors_actual, use_rep=representation, key_added=neighbors_key, random_state=random_state)
    sc.tl.diffmap(subset, neighbors_key=neighbors_key)

    finite_axis = np.isfinite(marker_axis)
    if finite_axis.sum() < 10:
        raise ValueError(f'{segment}: too few finite marker-axis values')
    early_cutoff = np.nanquantile(marker_axis[finite_axis], 0.05)
    root_candidates = np.flatnonzero(finite_axis & (marker_axis <= early_cutoff))
    embedding = np.asarray(subset.obsm[representation], dtype=float)
    root_centroid = np.nanmedian(embedding[root_candidates], axis=0)
    root_distance = np.linalg.norm(embedding[root_candidates] - root_centroid, axis=1)
    root_index = int(root_candidates[np.nanargmin(root_distance)])

    subset.uns['iroot'] = root_index
    sc.tl.dpt(subset, neighbors_key=neighbors_key)
    raw_dpt = subset.obs['dpt_pseudotime'].to_numpy(dtype=float)
    finite = np.isfinite(raw_dpt) & np.isfinite(marker_axis)
    if finite.sum() < 10:
        raise ValueError(f'{segment}: DPT produced too few finite values')

    orientation_corr = spearmanr(raw_dpt[finite], marker_axis[finite]).correlation
    oriented = raw_dpt.copy()
    if np.isfinite(orientation_corr) and orientation_corr < 0:
        oriented = -oriented
    finite_oriented = np.isfinite(oriented)
    p5, p95 = np.nanpercentile(oriented[finite_oriented], [5, 95])
    scaled = np.full(oriented.shape, np.nan, dtype=float)
    scaled[finite_oriented] = np.clip((oriented[finite_oriented] - p5) / max(p95 - p5, 1e-8), 0, 1)

    full_values = np.full(adata_obj.n_obs, np.nan, dtype=float)
    full_values[subset_mask] = scaled
    adata_obj.obs[output_col] = full_values

    return {
        'segment': segment,
        'n_tubules': int(subset.n_obs),
        'n_finite_dpt': int(np.isfinite(scaled).sum()),
        'representation': representation,
        'n_neighbors': int(n_neighbors_actual),
        'early_module': early_group,
        'late_module': late_group,
        'root_obs_name': str(subset.obs_names[root_index]),
        'root_marker_axis': float(marker_axis[root_index]),
        'raw_dpt_marker_spearman': float(orientation_corr),
        'scaled_dpt_marker_spearman': float(spearmanr(scaled[np.isfinite(scaled)], marker_axis[np.isfinite(scaled)]).correlation),
        'missing_markers': '; '.join(missing),
        'pseudotime_col': output_col,
    }


def save_total_pseudotime_anndata(
    source_adata,
    pseudotime_col,
    output_path,
    expression_filtered_adata,
    annotation_adata,
    project_dir,
    forbidden_labels=None,
    assert_no_nan=True,
    extra_obs_cols=None,
):
    if pseudotime_col not in source_adata.obs.columns:
        raise KeyError(f'Missing required pseudotime column: {pseudotime_col}')

    # source_adata (adata_total) already excludes Glomerulus/Vessel (see the continuum filter
    # above). Intersecting on obs_names and SUBSETTING out (not just NaN-filling the excluded
    # rows) is what makes the saved DPT h5ad contain only the retained (tubule-family +
    # Transitional) cells, per the requirement that this file never carries Glomerulus/Vessel.
    shared = expression_filtered_adata.obs_names.intersection(source_adata.obs_names)
    out = expression_filtered_adata[shared].copy()

    out.obs[pseudotime_col] = source_adata.obs.loc[shared, pseudotime_col].to_numpy(dtype=float)
    out.obs['total_segment_marker_call'] = source_adata.obs.loc[
        shared, 'total_segment_marker_call'
    ].astype(str).values

    if 'broad_tubule_marker_call' not in annotation_adata.obs.columns:
        raise KeyError('annotation_adata must carry broad_tubule_marker_call')
    for col in ['broad_tubule_marker_call', 'broad_tubule_marker_call_percell', 'celltype_primary', 'celltype_secondary']:
        if col in annotation_adata.obs.columns:
            out.obs[col] = annotation_adata.obs.loc[shared, col].astype(str).values

    if 'celltype_margin' in annotation_adata.obs.columns:
        out.obs['celltype_margin'] = annotation_adata.obs.loc[shared, 'celltype_margin'].to_numpy(dtype=float)

    keep_obs = [
        'sample', 'species', 'condition',
        'broad_tubule_marker_call', 'broad_tubule_marker_call_percell',
        'celltype_primary', 'celltype_secondary', 'celltype_margin',
        'total_segment_marker_call', pseudotime_col,
    ]
    if extra_obs_cols:
        keep_obs += [c for c in extra_obs_cols if c not in keep_obs]
    keep_obs = [col for col in keep_obs if col in out.obs.columns]
    out.obs = out.obs[keep_obs].copy()

    # --- Guards: the canonical pseudospace file must contain only nephron-tubule
    # families with a valid label (see notebook6_bioskills_review.md, the
    # Glomerulus/Glomeruli leak). These asserts fail loudly if excluded classes or
    # unlabeled cells slip through the continuum filter. ---
    labels = out.obs['broad_tubule_marker_call'].astype(str)
    if assert_no_nan:
        n_nan = int(labels.isin(['nan', 'NaN', 'None']).sum() + out.obs['broad_tubule_marker_call'].isna().sum())
        assert n_nan == 0, (
            f'{n_nan} tubules have a missing broad_tubule_marker_call in the DPT file; '
            'unlabeled cells must not leak into the pseudospace continuum.'
        )
    if forbidden_labels:
        present_forbidden = sorted(set(labels) & {str(x) for x in forbidden_labels})
        assert not present_forbidden, (
            f'DPT file contains forbidden classes {present_forbidden}; these must be '
            'excluded from the continuum before saving.'
        )

    out.write(output_path)
    print(
        f'Saved {pseudotime_col} AnnData to: {output_path.relative_to(project_dir)} '
        f'({out.n_obs:,} tubules; Glomerulus/Vessel excluded from this file)'
    )
    return out
