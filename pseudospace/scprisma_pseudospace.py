"""scPrisma linear-spectral pseudospace reconstruction (alternative to Scanpy DPT)
plus equal-count binning helpers.

Extracted verbatim from ``6_mouse_only_pseudospace.ipynb`` (Section 2). The
``scPrisma.algorithms`` dependency is imported lazily inside the functions that
use it, so importing this module never requires scPrisma to be on the path.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata, spearmanr
from sklearn.neighbors import NearestNeighbors


def _to_dense_float32(matrix):
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=np.float32)


def _sample_stratified_indices(
    adata,
    target_size,
    group_key="sample",
    random_state=0,
):
    """Draw a sample-stratified subsample with approximately preserved proportions."""
    if adata.n_obs <= target_size:
        return np.arange(adata.n_obs, dtype=int)

    if group_key not in adata.obs:
        raise KeyError(f"{group_key!r} is not present in adata.obs.")

    rng = np.random.default_rng(random_state)
    groups = adata.obs[group_key].astype(str)
    counts = groups.value_counts()
    proportions = counts / counts.sum()

    selected = []
    for group, proportion in proportions.items():
        group_idx = np.flatnonzero(groups.to_numpy() == group)
        n_group = max(1, int(round(target_size * proportion)))
        n_group = min(n_group, len(group_idx))
        selected.append(rng.choice(group_idx, size=n_group, replace=False))

    selected = np.concatenate(selected)

    if len(selected) > target_size:
        selected = rng.choice(selected, size=target_size, replace=False)
    elif len(selected) < target_size:
        remaining = np.setdiff1d(
            np.arange(adata.n_obs),
            selected,
            assume_unique=False,
        )
        n_extra = min(target_size - len(selected), len(remaining))
        selected = np.concatenate(
            [selected, rng.choice(remaining, size=n_extra, replace=False)]
        )

    return np.sort(selected.astype(int))


def _select_expression_features(adata, subsample_idx, n_genes):
    """
    Prefer the Harmony intersection-HVG feature set when available.

    Otherwise choose the most variable eligible genes within the scPrisma
    subsample. Mitochondrial/ribosomal genes are excluded when corresponding
    annotations are available.
    """
    eligible = np.ones(adata.n_vars, dtype=bool)

    if "exclude_from_harmony_hvg" in adata.var:
        eligible &= ~adata.var["exclude_from_harmony_hvg"].fillna(False).to_numpy(
            dtype=bool
        )

    if "highly_variable_for_harmony" in adata.var:
        harmony_hvg = adata.var["highly_variable_for_harmony"].fillna(False).to_numpy(
            dtype=bool
        )
        selected = eligible & harmony_hvg

        if selected.sum() >= 500:
            idx = np.flatnonzero(selected)

            # If more features are present than requested, rank by variance.
            if len(idx) > n_genes:
                layer = "lognorm" if "lognorm" in adata.layers else None
                matrix = (
                    adata.layers[layer][subsample_idx][:, idx]
                    if layer is not None
                    else adata.X[subsample_idx][:, idx]
                )
                values = _to_dense_float32(matrix)
                variance = np.nanvar(values, axis=0)
                idx = idx[np.argsort(variance)[-n_genes:]]

            return np.sort(idx)

    eligible_idx = np.flatnonzero(eligible)
    layer = "lognorm" if "lognorm" in adata.layers else None
    matrix = (
        adata.layers[layer][subsample_idx][:, eligible_idx]
        if layer is not None
        else adata.X[subsample_idx][:, eligible_idx]
    )
    values = _to_dense_float32(matrix)
    variance = np.nanvar(values, axis=0)

    n_keep = min(n_genes, len(eligible_idx))
    keep_local = np.argsort(variance)[-n_keep:]
    return np.sort(eligible_idx[keep_local])


def _prepare_expression_matrix(adata, subsample_idx, gene_idx):
    """
    Prepare cells × genes expression for scPrisma.

    Steps:
      1. Read log-normalized expression.
      2. Remove zero-variance genes.
      3. Z-score each gene across the subsample.
      4. L2-normalize each cell.
    """
    layer = "lognorm" if "lognorm" in adata.layers else None
    matrix = (
        adata.layers[layer][subsample_idx][:, gene_idx]
        if layer is not None
        else adata.X[subsample_idx][:, gene_idx]
    )
    values = _to_dense_float32(matrix)

    gene_mean = np.nanmean(values, axis=0, keepdims=True)
    gene_std = np.nanstd(values, axis=0, keepdims=True)

    valid = np.isfinite(gene_std[0]) & (gene_std[0] > 1e-8)
    values = values[:, valid]
    gene_idx = gene_idx[valid]
    gene_mean = gene_mean[:, valid]
    gene_std = gene_std[:, valid]

    values = (values - gene_mean) / gene_std
    values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)

    cell_norm = np.linalg.norm(values, axis=1, keepdims=True)
    cell_norm[cell_norm < 1e-8] = 1.0
    values /= cell_norm

    return np.ascontiguousarray(values, dtype=np.float64), gene_idx


def _random_bistochastic_matrix(n, rng, n_sinkhorn_iterations=100):
    """
    Generate a random positive doubly-stochastic initialization for E.

    A uniform E can be too symmetric, so use random positive entries followed
    by alternating row/column normalization (Sinkhorn normalization).
    """
    E = rng.random((n, n), dtype=np.float64) + 1e-8

    for _ in range(n_sinkhorn_iterations):
        E /= np.maximum(E.sum(axis=1, keepdims=True), 1e-12)
        E /= np.maximum(E.sum(axis=0, keepdims=True), 1e-12)

    return np.ascontiguousarray(E, dtype=np.float64)


def _extract_linear_eigenvectors(eig_data, n_cells):
    """
    Extract V from get_linear_eig_data output across API variants.

    V is expected to have n_cells rows and approximately n_cells - 1 columns.
    If a full n × n eigenvector matrix is returned, remove its most nearly
    constant column.
    """
    if isinstance(eig_data, (tuple, list)):
        candidates = [
            np.asarray(item, dtype=np.float64)
            for item in eig_data
            if np.asarray(item).ndim == 2
            and np.asarray(item).shape[0] == n_cells
        ]
    else:
        item = np.asarray(eig_data, dtype=np.float64)
        candidates = [item] if item.ndim == 2 and item.shape[0] == n_cells else []

    if not candidates:
        shapes = (
            [np.asarray(item).shape for item in eig_data]
            if isinstance(eig_data, (tuple, list))
            else [np.asarray(eig_data).shape]
        )
        raise ValueError(
            "Could not identify the linear eigenvector matrix V from "
            f"get_linear_eig_data output. Returned shapes: {shapes}"
        )

    V = max(candidates, key=lambda x: x.shape[1])

    if V.shape[1] == n_cells:
        # The constant eigenvector carries no positional information.
        constant_col = int(np.argmin(np.nanstd(V, axis=0)))
        V = np.delete(V, constant_col, axis=1)

    if V.shape != (n_cells, n_cells - 1):
        raise ValueError(
            f"Expected V shape {(n_cells, n_cells - 1)}, got {V.shape}."
        )

    return np.ascontiguousarray(V, dtype=np.float64)


def _hard_permutation(E):
    """Convert soft E to a hard permutation matrix when reconstruct_e is available."""
    from scPrisma.algorithms import reconstruct_e

    reconstructed = reconstruct_e(E)

    if isinstance(reconstructed, (tuple, list)):
        candidates = [
            np.asarray(item)
            for item in reconstructed
            if np.asarray(item).shape == E.shape
        ]
        if not candidates:
            raise ValueError("reconstruct_e returned no matrix matching E.shape.")
        permutation = candidates[0]
    else:
        permutation = np.asarray(reconstructed)

    if permutation.shape != E.shape:
        raise ValueError(
            f"reconstruct_e returned {permutation.shape}; expected {E.shape}."
        )

    return permutation


def _choose_permutation_orientation(E, permutation, marker_axis_sub):
    """
    Recover cell → linear-position mapping robustly.

    scPrisma versions differ in whether permutation rows encode cells or
    positions. Evaluate row-based, column-based, and soft-expectation mappings,
    then retain the candidate with the strongest absolute marker-axis
    correlation. This chooses matrix semantics, not biological direction.
    """
    n = E.shape[0]
    positions = np.arange(n, dtype=float)

    candidates = {
        "hard_row_argmax": np.argmax(permutation, axis=1).astype(float),
        "hard_column_argmax": np.argmax(permutation, axis=0).astype(float),
        "soft_row_expectation": E @ positions,
        "soft_column_expectation": E.T @ positions,
    }

    diagnostics = []
    for name, order in candidates.items():
        finite = np.isfinite(order) & np.isfinite(marker_axis_sub)
        rho = spearmanr(order[finite], marker_axis_sub[finite]).correlation
        diagnostics.append((name, float(rho), order))

    diagnostics.sort(
        key=lambda item: abs(item[1]) if np.isfinite(item[1]) else -np.inf,
        reverse=True,
    )

    best_name, best_rho, best_order = diagnostics[0]

    print("Permutation-orientation candidates:")
    for name, rho, _ in diagnostics:
        print(f"  {name:24s}: marker-axis rho={rho:+.4f}")
    print(f"Selected mapping: {best_name}")

    # Convert arbitrary positions to unique rank-based coordinates.
    ranked = rankdata(best_order, method="average") - 1.0
    return ranked, best_name, best_rho


def _propagate_ordering(X_full, X_sub, pseudo_sub, k):
    """Propagate subsample ordering to all cells through X_harmony kNN."""
    if k >= len(X_sub):
        k = max(1, len(X_sub) - 1)

    neighbors = NearestNeighbors(n_neighbors=k, algorithm="auto")
    neighbors.fit(X_sub)

    distances, indices = neighbors.kneighbors(X_full)

    # Inverse-distance weighting is smoother than an unweighted mean.
    weights = 1.0 / np.maximum(distances, 1e-8)
    propagated = np.sum(weights * pseudo_sub[indices], axis=1) / np.sum(
        weights,
        axis=1,
    )

    return propagated


def _orient_and_scale(pseudotime, marker_axis):
    """Orient early→late using marker axis, then p5–p95 scale to [0,1]."""
    finite = np.isfinite(pseudotime) & np.isfinite(marker_axis)
    if finite.sum() < 10:
        raise ValueError("Too few finite values to orient scPrisma ordering.")

    rho_before = spearmanr(
        pseudotime[finite],
        marker_axis[finite],
    ).correlation

    oriented = pseudotime.copy()
    flipped = bool(np.isfinite(rho_before) and rho_before < 0)
    if flipped:
        oriented = -oriented

    p5, p95 = np.nanpercentile(oriented[finite], [5, 95])
    scaled = np.full_like(oriented, np.nan, dtype=float)
    scaled[finite] = np.clip(
        (oriented[finite] - p5) / max(p95 - p5, 1e-8),
        0.0,
        1.0,
    )

    rho_after = spearmanr(
        scaled[finite],
        marker_axis[finite],
    ).correlation

    return scaled, float(rho_before), float(rho_after), flipped


def run_scprisma_linear_pseudospace(
    adata,
    config,
):
    """Run subsampled scPrisma linear reconstruction and propagate to all cells."""
    from scPrisma.algorithms import get_linear_eig_data, sga_matrix_momentum

    print("=" * 78)
    print("scPrisma linear spectral reconstruction")
    print("=" * 78)

    marker_col = config["marker_axis_col"]
    if marker_col not in adata.obs:
        raise KeyError(
            f"{marker_col!r} is not present in adata.obs. Compute the same "
            "early-versus-late marker axis used to orient DPT before running "
            "this block."
        )

    if "X_harmony" not in adata.obsm:
        raise KeyError("X_harmony is required for full-data kNN propagation.")

    rng = np.random.default_rng(config["random_state"])

    sub_idx = _sample_stratified_indices(
        adata,
        target_size=config["subsample_size"],
        group_key="sample",
        random_state=config["random_state"],
    )

    gene_idx = _select_expression_features(
        adata,
        subsample_idx=sub_idx,
        n_genes=config["n_genes"],
    )

    A, gene_idx = _prepare_expression_matrix(
        adata,
        subsample_idx=sub_idx,
        gene_idx=gene_idx,
    )

    n_cells, n_genes = A.shape
    marker_axis_full = adata.obs[marker_col].to_numpy(dtype=float)
    marker_axis_sub = marker_axis_full[sub_idx]

    print(f"Full object: {adata.n_obs:,} tubules × {adata.n_vars:,} genes")
    print(f"scPrisma subsample: {n_cells:,} tubules × {n_genes:,} genes")
    print(f"Linear alpha: {config['alpha']}")
    print(f"Eigen method: {config['eig_method']!r}")
    print(
        f"Estimated size of one n×n float64 matrix: "
        f"{n_cells * n_cells * 8 / 1024**3:.3f} GiB"
    )

    eig_data = get_linear_eig_data(
        ncells=n_cells,
        alpha=config["alpha"],
        method=config["eig_method"],
        normalize_vectors=config["normalize_vectors"],
    )
    V = _extract_linear_eigenvectors(eig_data, n_cells)
    print(f"Linear spectral template V shape: {V.shape}")

    E0 = _random_bistochastic_matrix(n_cells, rng)

    print("Running sga_matrix_momentum; first call includes Numba compilation...")
    E_opt = sga_matrix_momentum(
        A=A,
        E=E0,
        V=V,
        iterNum=config["iterations"],
        batch_size=min(config["batch_size"], n_genes),
        lr=config["learning_rate"],
        gamma=config["momentum"],
        verbose=True,
    )
    E_opt = np.asarray(E_opt, dtype=np.float64)

    if E_opt.shape != (n_cells, n_cells):
        raise ValueError(
            f"sga_matrix_momentum returned {E_opt.shape}; "
            f"expected {(n_cells, n_cells)}."
        )

    permutation = _hard_permutation(E_opt)
    pseudo_sub, mapping_name, mapping_rho = _choose_permutation_orientation(
        E_opt,
        permutation,
        marker_axis_sub,
    )

    pseudo_sub /= max(float(np.nanmax(pseudo_sub)), 1.0)

    harmony_full = np.asarray(adata.obsm["X_harmony"], dtype=np.float64)
    harmony_sub = harmony_full[sub_idx]

    pseudo_full = _propagate_ordering(
        harmony_full,
        harmony_sub,
        pseudo_sub,
        k=config["knn_propagate"],
    )

    # Preserve exact scPrisma coordinates for cells used in the spectral solve.
    pseudo_full[sub_idx] = pseudo_sub

    scaled, rho_before, rho_after, flipped = _orient_and_scale(
        pseudo_full,
        marker_axis_full,
    )

    adata.obs[config["output_col"]] = scaled
    adata.obs[f"{config['output_col']}_in_scprisma_subsample"] = False
    adata.obs.loc[
        adata.obs_names[sub_idx],
        f"{config['output_col']}_in_scprisma_subsample",
    ] = True

    print()
    print("scPrisma reconstruction summary:")
    print(f"  permutation mapping: {mapping_name}")
    print(f"  subsample mapping rho: {mapping_rho:+.4f}")
    print(f"  full-data marker rho before orientation: {rho_before:+.4f}")
    print(f"  full-data marker rho after orientation:  {rho_after:+.4f}")
    print(f"  direction flipped: {flipped}")
    print(
        f"  {config['output_col']} range: "
        f"{np.nanmin(scaled):.4f}–{np.nanmax(scaled):.4f}"
    )

    if "total_scanpy_dpt" in adata.obs:
        dpt = adata.obs["total_scanpy_dpt"].to_numpy(dtype=float)
        finite = np.isfinite(dpt) & np.isfinite(scaled)
        dpt_rho = spearmanr(dpt[finite], scaled[finite]).correlation
        print(f"  Spearman(scPrisma, DPT): {dpt_rho:+.4f}")
    else:
        dpt_rho = np.nan

    return {
        "subsample_indices": sub_idx,
        "selected_gene_indices": gene_idx,
        "mapping_name": mapping_name,
        "mapping_spearman": mapping_rho,
        "orientation_spearman_before": rho_before,
        "orientation_spearman_after": rho_after,
        "dpt_spearman": float(dpt_rho),
        "output_col": config["output_col"],
    }


def equal_count_bins(values, n_bins=40):
    """Assign approximately equal-count bins along pseudospace."""
    percentile_rank = (
        pd.Series(values)
        .rank(method="average", pct=True)
        .to_numpy()
    )
    bins = np.floor(percentile_rank * n_bins).astype(int)
    return np.clip(bins, 0, n_bins - 1)


def segment_composition(df, pseudospace_col, segment_col, segment_order, n_bins=40):
    """Segment proportions across equal-count scPrisma bins."""
    local = df[[pseudospace_col, segment_col]].dropna().copy()
    local["bin"] = equal_count_bins(local[pseudospace_col], n_bins)

    counts = pd.crosstab(local[segment_col], local["bin"])
    counts = counts.reindex(segment_order, fill_value=0)

    return counts.div(
        counts.sum(axis=0).replace(0, np.nan),
        axis=1,
    ).fillna(0)
