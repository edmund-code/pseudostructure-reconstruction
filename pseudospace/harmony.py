"""Harmony integration helpers: per-condition HVG selection and the rpy2 bridge.

Extracted verbatim from ``6_mouse_only_pseudospace.ipynb`` (Section 1). ``rpy2``
is imported lazily inside :func:`run_harmony_rpy2` (and, transitively,
:func:`_verify_harmony_bridge`) so importing this module does not require
rpy2/R to be available.

:func:`run_harmony_rpy2` bridges to R's ``harmony::RunHarmony`` (harmony >= 2.0;
the older ``harmony::HarmonyMatrix`` entry point this was originally written
against is gone). Call :func:`_verify_harmony_bridge` explicitly (from a smoke
test or a notebook cell -- it is NOT run automatically inside
``run_harmony_rpy2`` since it's too slow for full-size data) to sanity-check
that the bridge is not silently scrambling or mis-orienting the PCA matrix
before trusting a full run.
"""
from __future__ import annotations

import re
from typing import Literal, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


# The bridge below is written for the modern RunHarmony API. Older 1.x releases
# may accept a superficially similar call but do not reproduce the same embedding.
MIN_HARMONY_VERSION = (2, 0, 5)


def _parse_version(version: str) -> tuple[int, ...]:
    """Return a numeric package-version tuple, rejecting ambiguous strings."""
    normalized = re.sub(r"^\[\d+\]\s*", "", version.strip())
    match = re.match(r"^(\d+(?:\.\d+)*)", normalized)
    if not match:
        raise ValueError(f"Cannot parse R harmony version {version!r}.")
    text = match.group(1)
    if "." in text:
        return tuple(int(part) for part in text.split("."))
    return tuple(int(part) for part in normalized.split())


def require_supported_harmony_version(
    version: str,
    minimum: tuple[int, ...] = MIN_HARMONY_VERSION,
) -> None:
    """Fail before integration when R's ``harmony`` package is too old."""
    parsed = _parse_version(version)
    padded = parsed + (0,) * max(0, len(minimum) - len(parsed))
    if padded < minimum:
        required = ".".join(map(str, minimum))
        raise RuntimeError(
            f"R package 'harmony' {version} is unsupported; this workflow requires "
            f"harmony >= {required}. Use the conda environment's pinned r-harmony package "
            "and do not prepend a user-level R library to .libPaths()."
        )


def select_harmony_hvgs_by_condition(
    adata: ad.AnnData,
    group_key: str = 'condition',
    groups: tuple[str, ...] = ('Control', 'IR'),
    mode: Literal['union', 'intersection', 'custom'] = 'intersection',
    min_groups: Optional[int] = None,
    min_mean: float = 0.0125,
    max_mean: float = 3,
    min_disp: float = 0.5,
) -> ad.AnnData:
    """
    Select Harmony HVGs by running HVG selection separately within each condition group,
    then taking genes that are HVG in at least `min_groups` groups.

    `mode` chooses how strict the aggregation is:
      - 'intersection' (default, recommended): a gene must be HVG in EVERY group.
        Yields condition-invariant features -- Harmony integrates on cell-type
        structure shared across conditions, avoiding over-correction of the
        condition-associated genes tested downstream (e.g. in the AKI comparison).
      - 'union': a gene qualifies if it is HVG in at least ONE group. Larger
        feature set, but risks feeding Harmony condition-associated genes and
        partially removing the biological signal of interest.
      - 'custom': use the value passed via `min_groups` directly.

    Backward compat: if `min_groups` is passed explicitly, it overrides `mode`
    (with a printed note). If neither `mode` nor `min_groups` is set, defaults to
    'intersection'.

    Adapted from notebook 1 Part 3 (which used min_groups=2 across 3 groups:
    Human/Control/IR -- effectively "HVG in >=2 of 3", closer to intersection than
    union). For 2 groups (Control/IR), intersection == min_groups=2.

    Provenance note: the ``highly_variable_for_harmony`` gene set selected here
    is the feature set the subsequent ``sc.pp.pca()`` call should use to build
    ``adata.obsm['X_pca']`` -- keep them in sync, since :func:`run_harmony_rpy2`
    can only log PCA provenance from ``adata.uns['pca']['params']``, not verify
    it was built from this exact gene set.
    """
    adata = adata.copy()

    # Resolve mode -> min_groups
    if min_groups is not None:
        resolved_min_groups = int(min_groups)
        resolved_mode = 'custom' if mode == 'custom' else f'custom (min_groups={min_groups} overrides mode={mode!r})'
    elif mode == 'intersection':
        resolved_min_groups = len(groups)
        resolved_mode = 'intersection'
    elif mode == 'union':
        resolved_min_groups = 1
        resolved_mode = 'union'
    elif mode == 'custom':
        raise ValueError("mode='custom' requires min_groups to be set explicitly.")
    else:
        raise ValueError(f"Unknown mode: {mode!r}. Use 'union', 'intersection', or 'custom'.")

    eligible_mask = ~adata.var['exclude_from_harmony_hvg'].values
    eligible_genes = set(adata.var_names[eligible_mask])

    hvg_counts = pd.Series(0, index=adata.var_names, dtype=int)
    per_group_hvgs: dict[str, set] = {}

    print('\n' + '=' * 60)
    print('SELECTING HARMONY HVGs WITHIN EACH CONDITION GROUP')
    print('=' * 60)
    print(f'Groups: {groups}')
    print(f'Mode: {resolved_mode}')
    print(f'Minimum groups required (>= this many groups must call HVG): {resolved_min_groups}')
    print(f'Eligible non-mito/ribo genes: {len(eligible_genes)}')

    for group in groups:
        idx = adata.obs[group_key].astype(str).values == group

        if idx.sum() < 20:
            print(f'  Skipping {group}: only {idx.sum()} tubules')
            per_group_hvgs[group] = set()
            continue

        sub = adata[idx, list(eligible_genes)].copy()
        sub.X = sub.layers['lognorm'].copy()

        unique_samples = sub.obs['sample'].astype(str).nunique()
        batch_key = 'sample' if unique_samples > 1 else None

        print(
            f'\n  {group}: {sub.n_obs} tubules x {sub.n_vars} eligible genes '
            f'({unique_samples} sample(s), batch_key={batch_key})'
        )

        sc.pp.highly_variable_genes(
            sub,
            min_mean=min_mean,
            max_mean=max_mean,
            min_disp=min_disp,
            batch_key=batch_key,
            flavor='seurat',
        )

        group_hvgs = set(sub.var_names[sub.var['highly_variable']])
        per_group_hvgs[group] = group_hvgs
        hvg_counts.loc[list(group_hvgs)] += 1

        print(f'    HVGs in {group}: {len(group_hvgs)}')

    selected = hvg_counts >= resolved_min_groups
    selected = selected & ~adata.var['exclude_from_harmony_hvg']

    adata.var['hvg_condition_count'] = hvg_counts.values
    adata.var['highly_variable_for_harmony'] = selected.values

    union_genes = set.union(*per_group_hvgs.values()) if per_group_hvgs else set()
    nonempty_group_sets = [g for g in per_group_hvgs.values() if g]
    intersection_genes = (
        set.intersection(*nonempty_group_sets)
        if len(nonempty_group_sets) == len(per_group_hvgs) and nonempty_group_sets
        else set()
    )

    n_selected = int(adata.var['highly_variable_for_harmony'].sum())

    print('\nHarmony HVG summary:')
    for group in groups:
        print(f'  HVGs in {group}: {len(per_group_hvgs.get(group, set()))}')
    print(f'  Union (HVG in at least one group): {len(union_genes)}')
    print(f'  Intersection (HVG in every group): {len(intersection_genes)}')
    print(
        f'  Final Harmony feature genes (>= {resolved_min_groups} group(s), non-mito/ribo, mode={resolved_mode}): '
        f'{n_selected}'
    )

    # Warnings scaled to what the mode is doing
    if n_selected < 500:
        print(
            f'\nWARNING: only {n_selected} Harmony features selected. '
            f'Consider lowering min_disp or switching to mode=\'union\' if this is too restrictive.'
        )
    elif resolved_mode == 'union' and n_selected > 5000:
        print(
            f'\nNOTE: {n_selected} Harmony features selected under mode=\'union\'. '
            f"Consider mode='intersection' to avoid feeding condition-associated genes to Harmony "
            f"(intersection would give {len(intersection_genes)} features here)."
        )

    return adata


def sweep_harmony_theta(
    adata: ad.AnnData,
    batch_key: str,
    bio_label_key: str,
    thetas,
    n_pcs: int = 50,
    lambda_val: float = 1,
    max_iter: int = 30,
    tau: int = 1000,
    n_neighbors: int = 25,
    subsample: int | None = 10_000,
    random_state: int = 0,
):
    """Run Harmony at several ``theta`` values and report the batch-mixing vs
    bio-conservation trade-off so a value can be chosen deliberately.

    For each theta this runs :func:`run_harmony_rpy2` on ``adata`` (which must
    already carry ``obsm['X_pca']``), then scores:

      - ``batch_mixing`` -- mean over cells of the normalized Shannon entropy of
        the ``batch_key`` composition among each cell's kNN in the corrected
        space. Higher = batches better mixed (over-correction raises this).
      - ``bio_conservation`` -- silhouette of ``bio_label_key`` on the corrected
        embedding. Higher = biological structure better preserved.

    Aggressive theta raises mixing and lowers bio-conservation in lockstep, so the
    right value is the smallest theta that mixes samples without collapsing the
    ``bio_label_key`` structure -- inspect the returned table, do not maximize
    mixing alone.

    Returns a DataFrame with one row per theta. Does NOT mutate ``adata`` (each
    Harmony run is done on a copy); ``obsm['X_harmony']`` is not written back.
    """
    import numpy as _np
    from sklearn.metrics import silhouette_score
    from sklearn.neighbors import NearestNeighbors

    rng = _np.random.default_rng(random_state)
    batch_all = adata.obs[batch_key].astype(str).to_numpy()
    bio_all = adata.obs[bio_label_key].astype(str).to_numpy()

    if subsample is not None and adata.n_obs > subsample:
        eval_idx = _np.sort(rng.choice(adata.n_obs, size=subsample, replace=False))
    else:
        eval_idx = _np.arange(adata.n_obs)

    def _mean_neighbor_entropy(embedding, labels):
        nn = NearestNeighbors(n_neighbors=min(n_neighbors, embedding.shape[0] - 1))
        nn.fit(embedding)
        _, idx = nn.kneighbors(embedding)
        uniq = _np.unique(labels)
        log_k = _np.log(len(uniq)) if len(uniq) > 1 else 1.0
        entropies = []
        for row in idx:
            counts = pd.Series(labels[row]).value_counts().to_numpy(dtype=float)
            p = counts / counts.sum()
            entropies.append(-_np.sum(p * _np.log(p)) / log_k)
        return float(_np.mean(entropies))

    rows = []
    for theta in thetas:
        result = run_harmony_rpy2(
            adata.copy(), batch_key=batch_key, n_pcs=n_pcs, theta=float(theta),
            lambda_val=lambda_val, max_iter=max_iter, tau=tau,
        )
        emb = _np.asarray(result.obsm['X_harmony'])[eval_idx]
        batch_mixing = _mean_neighbor_entropy(emb, batch_all[eval_idx])
        bio_labels = bio_all[eval_idx]
        bio_conservation = (
            float(silhouette_score(emb, bio_labels))
            if pd.Series(bio_labels).nunique() > 1 else _np.nan
        )
        rows.append({
            'theta': float(theta),
            'batch_mixing_entropy': batch_mixing,
            'bio_conservation_silhouette': bio_conservation,
            'n_eval_cells': int(len(eval_idx)),
        })
        print(f'  theta={theta}: batch_mixing={batch_mixing:.3f}, '
              f'bio_conservation={bio_conservation:.3f}')

    return pd.DataFrame(rows)


def _per_column_correlation(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Pearson correlation between matched columns of two same-shape arrays.

    Used as a provenance/scramble check: if Harmony's output has (near-)zero
    or wildly inconsistent per-column correlation with the input PCA it was
    given, that is a strong signal the matrix was mis-oriented or its row
    order was scrambled somewhere in the R round-trip, rather than
    legitimately batch-corrected (which perturbs but does not decorrelate
    each PC axis outright).
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError(f'Shape mismatch for per-column correlation: {a.shape} vs {b.shape}')
    a_centered = a - a.mean(axis=0, keepdims=True)
    b_centered = b - b.mean(axis=0, keepdims=True)
    numerator = (a_centered * b_centered).sum(axis=0)
    denom = np.sqrt((a_centered ** 2).sum(axis=0) * (b_centered ** 2).sum(axis=0))
    with np.errstate(invalid='ignore', divide='ignore'):
        corr = numerator / denom
    return corr


def _assert_r_matrix_orientation(ro, r_matrix, py_array: np.ndarray, n_checks: int = 8) -> None:
    """Spot-check that an R matrix has the same (row, col) orientation as ``py_array``.

    Samples the four corners plus a handful of random interior points and
    compares ``r_matrix[i+1, j+1]`` (R is 1-based) against ``py_array[i, j]``.
    Raises immediately on any mismatch -- this is the direct defense against a
    silent transpose/scramble during the numpy -> R matrix conversion.
    """
    rng = np.random.default_rng(0)
    n_rows, n_cols = py_array.shape
    check_pts = {(0, 0), (0, n_cols - 1), (n_rows - 1, 0), (n_rows - 1, n_cols - 1)}
    max_pts = min(n_checks, n_rows * n_cols)
    while len(check_pts) < max_pts:
        check_pts.add((int(rng.integers(0, n_rows)), int(rng.integers(0, n_cols))))

    get_elem = ro.r('function(m, i, j) m[i, j]')
    for (i, j) in check_pts:
        r_val = float(get_elem(r_matrix, i + 1, j + 1)[0])
        py_val = float(py_array[i, j])
        if not np.isclose(r_val, py_val, atol=1e-8, rtol=1e-6):
            raise ValueError(
                f'R matrix orientation check FAILED at (row={i}, col={j}): '
                f'R[{i + 1}, {j + 1}]={r_val!r}, Python[{i}, {j}]={py_val!r}. The PCA matrix '
                f'appears to have been transposed or scrambled during R conversion.'
            )


def run_harmony_rpy2(
    adata: ad.AnnData,
    batch_key: str,
    n_pcs: int = 40,
    theta: float = 2,
    lambda_val: float = 1,
    max_iter: int = 30,
    tau: int = 1000,  # NOT harmony>=2.0's own default (0) -- preserved intentionally, see docstring.
    ncores: int = 1,
) -> ad.AnnData:
    """Bridges to R's harmony::RunHarmony via rpy2.

    Updated from the original ``harmony::HarmonyMatrix`` bridge (verbatim from
    notebook 1 Part 3) because newer installs of the R ``harmony`` package no
    longer export the legacy ``HarmonyMatrix`` wrapper -- ``RunHarmony`` is the
    current entry point. ``RunHarmony``'s ``data_mat`` argument is cells x PCs
    (no transpose), unlike the old ``HarmonyMatrix`` which wanted PCs x cells.

    Provenance requirements (not enforced beyond what's noted): ``adata.obsm
    ['X_pca']`` must already be computed on the intended HVG feature set (see
    :func:`select_harmony_hvgs_by_condition`), and ``n_pcs`` must match
    whatever was validated in the trusted pipeline -- this function does not
    know the "right" n_pcs, it only refuses to silently use fewer PCs than
    requested (see the availability check below) and logs whatever PCA
    provenance is available in ``adata.uns['pca']['params']``.

    ``tau`` default note: this function's own default is ``tau=1000``,
    inherited unchanged from the original ``harmony::HarmonyMatrix`` pipeline.
    harmony >= 2.0's own ``harmony_options()`` default is ``tau=0`` -- these
    are materially different regularization strengths. That default is
    preserved here deliberately (not silently changed to match the new
    package default); it is surfaced in the diagnostics print below so any
    caller relying on the implicit old value can see it.

    Call :func:`_verify_harmony_bridge` separately (not from inside this
    function -- it's too slow to run on every call against full-size data) to
    sanity-check the bridge is not scrambling/mis-orienting data before
    trusting a real run.
    """
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr

    harmony = importr('harmony')
    base = importr('base')
    utils = importr('utils')

    # --- Provenance / availability guards -----------------------------------
    available_pcs = adata.obsm['X_pca'].shape[1]
    if n_pcs > available_pcs:
        raise ValueError(
            f'Requested n_pcs={n_pcs} but adata.obsm["X_pca"] only has {available_pcs} PCs '
            f'computed. Silently using fewer PCs than requested would change Harmony\'s '
            f'integration behavior with no visible signal -- re-run PCA with n_comps >= {n_pcs} '
            f'or lower n_pcs explicitly.'
        )

    pca_params = adata.uns.get('pca', {}).get('params') if 'pca' in adata.uns else None

    pca_embeddings = np.ascontiguousarray(adata.obsm['X_pca'][:, :n_pcs])
    n_cells = pca_embeddings.shape[0]
    n_pcs_actual = pca_embeddings.shape[1]

    # --- Metadata alignment guards ------------------------------------------
    if len(adata.obs_names) != n_cells:
        raise ValueError(
            f'adata.obs_names length ({len(adata.obs_names)}) does not match '
            f'adata.obsm["X_pca"] row count ({n_cells}); AnnData is misaligned.'
        )
    n_missing_batch = int(adata.obs[batch_key].isna().sum())
    if n_missing_batch > 0:
        raise ValueError(
            f"batch_key='{batch_key}' has {n_missing_batch} missing (NaN) value(s); "
            f'Harmony requires a complete batch label for every cell.'
        )
    batch_labels = adata.obs[batch_key].astype(str).values
    assert len(batch_labels) == len(adata.obs_names) == n_cells, (
        f'batch_labels ({len(batch_labels)}), adata.obs_names ({len(adata.obs_names)}), and '
        f'PCA row count ({n_cells}) must all match in length and order.'
    )

    # --- Diagnostics ---------------------------------------------------------
    try:
        harmony_version = str(utils.packageVersion('harmony')[0])
    except Exception as exc:
        raise RuntimeError("Could not determine the installed R package 'harmony' version.") from exc
    require_supported_harmony_version(harmony_version)

    try:
        exported_functions = sorted(str(x) for x in base.getNamespaceExports('harmony'))
    except Exception:
        exported_functions = []

    batch_value_counts = pd.Series(batch_labels).value_counts()

    print('=' * 60)
    print('HARMONY (rpy2) DIAGNOSTICS')
    print('=' * 60)
    print(f'  R harmony package version: {harmony_version}')
    print(f'  Exported harmony functions: {exported_functions}')
    print(f'  PCA matrix shape (Python, cells x PCs): {pca_embeddings.shape}')
    if pca_params is not None:
        print(f'  PCA provenance (adata.uns["pca"]["params"]): {pca_params}')
    else:
        print('  PCA provenance: adata.uns["pca"]["params"] not found -- cannot confirm which '
              'feature set/params produced adata.obsm["X_pca"].')
    print(f'  Batch key: {batch_key}')
    print(f'  Batch value counts:\n{batch_value_counts.to_string()}')

    if 'RunHarmony' not in exported_functions:
        raise AttributeError(
            f"Installed harmony package (version {harmony_version}) does not export "
            f"'RunHarmony'. Exported functions: {exported_functions}"
        )

    # Inspect the installed RunHarmony signature so we only pass arguments it
    # actually accepts (e.g. some versions don't expose `do_pca`, and harmony
    # >= 2.0 moved `tau` behind `.options=harmony_options(tau=...)`).
    def _r_formal_names(expr: str) -> set[str]:
        try:
            names = ro.r(
                f'tryCatch(names(formals({expr})), error=function(e) character(0))'
            )
            return set(names)
        except Exception:
            return set()

    accepted_args = _r_formal_names('harmony:::RunHarmony.default')
    if not accepted_args or accepted_args == {'...'}:
        accepted_args = _r_formal_names('harmony::RunHarmony')
    print(f'  RunHarmony.default accepted arguments: {sorted(accepted_args)}')

    print(f'  Requested: theta={theta}, lambda={lambda_val}, max_iter={max_iter}, tau={tau} '
          f'(harmony>=2.0 package default for tau is 0 -- see docstring note).')
    print(f'  Requested ncores={ncores} (pinned to one core for reproducibility).')

    # Resolve tau routing BEFORE entering the localconverter context below.
    # harmony_options() must be built outside that context: calling it *inside*
    # a numpy2ri/pandas2ri localconverter auto-converts its return value to a
    # plain Python list, stripping the R "harmony_options" S3 class -- RunHarmony
    # then rejects it with "'.options' must be created from harmony_options()!".
    # Built here (default converter only), the class survives, and the object
    # is passed into the context below as an already-built R object (untouched
    # by conversion since it's not itself being converted, just passed through).
    if 'tau' in accepted_args:
        tau_route = 'top-level tau= argument'
        r_harmony_options = None
    elif '.options' in accepted_args:
        tau_route = '.options=harmony_options(tau=...)'
        r_harmony_options = harmony.harmony_options(tau=float(tau))
    else:
        raise RuntimeError(
            f'Requested tau={tau} but installed RunHarmony (harmony {harmony_version}) exposes '
            f"neither a top-level 'tau' formal nor '.options'. Refusing to silently drop the "
            f'requested tau -- accepted arguments: {sorted(accepted_args)}'
        )

    with localconverter(ro.default_converter + numpy2ri.converter + pandas2ri.converter):
        # Let the active numpy2ri converter build the R matrix directly from the
        # contiguous cells x PCs array -- no manual flatten/byrow bookkeeping,
        # which is exactly the kind of hand-rolled invariant (order='C' +
        # byrow=True must agree) that can silently transpose or scramble data
        # on a future edit while still producing a plausible-looking result.
        r_pcs = ro.conversion.py2rpy(pca_embeddings)

        r_pcs_nrow, r_pcs_ncol = int(base.nrow(r_pcs)[0]), int(base.ncol(r_pcs)[0])
        if r_pcs_nrow != n_cells or r_pcs_ncol != n_pcs_actual:
            raise ValueError(
                f'R PCA matrix has dimensions ({r_pcs_nrow}, {r_pcs_ncol}), expected '
                f'({n_cells}, {n_pcs_actual}) [cells x PCs]. Conversion did not preserve shape.'
            )
        # Direct element-wise proof the conversion preserved orientation (not just shape).
        _assert_r_matrix_orientation(ro, r_pcs, pca_embeddings)

        r_meta = ro.DataFrame({batch_key: ro.StrVector(batch_labels)})
        r_meta.rownames = ro.StrVector(list(adata.obs_names))

        r_meta_nrow = int(base.nrow(r_meta)[0])
        if r_meta_nrow != n_cells:
            raise ValueError(f'R metadata has {r_meta_nrow} rows, expected {n_cells}.')
        r_meta_rownames = list(r_meta.rownames)
        if r_meta_rownames != list(adata.obs_names):
            raise ValueError('R metadata row names do not match adata.obs_names order.')

        print(f'  R PCA matrix dimensions (cells x PCs): {r_pcs_nrow} x {r_pcs_ncol} -- '
              f'orientation verified element-wise.')
        print(f'  R metadata dimensions: {r_meta_nrow} x {int(base.ncol(r_meta)[0])} -- '
              f'row names verified to match adata.obs_names order.')

        if r_harmony_options is not None:
            is_valid_options = bool(
                ro.r('function(x) inherits(x, "harmony_options")')(r_harmony_options)[0]
            )
            if not is_valid_options:
                raise RuntimeError(
                    "Internal error: the harmony_options() object lost its R 'harmony_options' "
                    "S3 class before reaching RunHarmony (this is exactly the localconverter "
                    "class-stripping bug described above) -- refusing to call RunHarmony with a "
                    "corrupted .options object."
                )

        # NOTE: RunHarmony's only top-level R formal is `(...)` -- it dispatches
        # via UseMethod, so rpy2's usual reserved-word translation (lambda_ ->
        # lambda) never triggers (that only fires when the wrapped function's
        # own formals contain the clashing name). The literal string "lambda_"
        # would otherwise pass straight through to R's internal arg-checker and
        # be rejected, so the R name "lambda" is set directly via dict key.
        run_harmony_kwargs = {
            'vars_use': ro.StrVector([batch_key]),
            'theta': float(theta),
            'lambda': float(lambda_val),
            'max_iter': int(max_iter),
            'ncores': int(ncores),
            'verbose': True,
        }
        if 'do_pca' in accepted_args:
            run_harmony_kwargs['do_pca'] = False
        else:
            print("  Note: installed RunHarmony has no 'do_pca' arg; omitting it "
                  "(PCA embeddings are already being passed in directly).")

        if tau_route == 'top-level tau= argument':
            run_harmony_kwargs['tau'] = float(tau)
        else:
            run_harmony_kwargs['.options'] = r_harmony_options

        if 'return_object' in accepted_args:
            run_harmony_kwargs['return_object'] = False

        effective_params_log = {
            'vars_use': [batch_key],
            'theta': float(theta),
            'lambda': float(lambda_val),
            'max_iter': int(max_iter),
            'ncores': int(ncores),
            'do_pca': run_harmony_kwargs.get('do_pca', '<omitted, arg not accepted>'),
            'tau_requested': float(tau),
            'tau_route': tau_route,
            'return_object': run_harmony_kwargs.get('return_object', '<omitted, arg not accepted>'),
        }
        print(f'  Effective RunHarmony call parameters: {effective_params_log}')

        harmony_result = harmony.RunHarmony(r_pcs, r_meta, **run_harmony_kwargs)

        harmony_array = np.array(harmony_result)

    result_shape = harmony_array.shape
    print(f'  Harmony result shape (as returned by R): {result_shape}')

    # Output orientation: prefer explicit, provenance-based resolution (Harmony
    # returns cells in the same order as meta_data) over shape guessing. A
    # square matrix makes shape alone ambiguous -- refuse to guess rather than
    # silently transpose.
    if n_cells == n_pcs_actual:
        raise ValueError(
            f'Cannot safely determine Harmony output orientation: n_cells == n_pcs_actual == '
            f'{n_cells} (square result). Shape alone cannot disambiguate cells x PCs from PCs x '
            f'cells in this degenerate case, and this function will not guess. Re-run with a '
            f'different n_pcs (so n_cells != n_pcs_actual) or verify orientation manually.'
        )
    if result_shape == (n_cells, n_pcs_actual):
        harmony_embeddings = harmony_array
    elif result_shape == (n_pcs_actual, n_cells):
        harmony_embeddings = harmony_array.T
        print('  Note: Harmony returned PCs x cells; transposed to cells x PCs.')
    else:
        raise ValueError(
            f'Unexpected Harmony result shape {result_shape}; expected either '
            f'({n_cells}, {n_pcs_actual}) [cells x PCs] or '
            f'({n_pcs_actual}, {n_cells}) [PCs x cells].'
        )

    # Provenance check, independent of the shape-based resolution above: confirm
    # the accepted orientation's columns still correlate with the input PCA.
    # Legitimate batch correction perturbs values but should not decorrelate
    # every PC axis outright -- near-zero correlation indicates scrambling.
    col_corrs = _per_column_correlation(pca_embeddings, harmony_embeddings)
    mean_corr = float(np.nanmean(col_corrs))
    print(f'  Per-column correlation(input PCA, Harmony output): mean={mean_corr:.3f}, '
          f'min={float(np.nanmin(col_corrs)):.3f}, max={float(np.nanmax(col_corrs)):.3f}')
    if mean_corr < 0.1:
        raise ValueError(
            f'Harmony output barely correlates with the input PCA per-column '
            f'(mean r={mean_corr:.3f} across {len(col_corrs)} PCs). This indicates a cell-'
            f'ordering/orientation problem rather than legitimate batch correction -- refusing '
            f'to accept the result.'
        )

    print(f'  Harmony done. Final embedding shape: {harmony_embeddings.shape}')

    assert harmony_embeddings.shape == (n_cells, n_pcs_actual), (
        f'Shape mismatch: expected ({n_cells}, {n_pcs_actual}), got {harmony_embeddings.shape}'
    )

    adata.obsm['X_pca_harmony'] = harmony_embeddings
    # Also written under the legacy key: sweep_harmony_theta() and every
    # downstream notebook/module (trajectory.py, scprisma_pseudospace.py,
    # notebooks 2/6) read 'X_harmony', not 'X_pca_harmony'.
    adata.obsm['X_harmony'] = harmony_embeddings

    return adata


def _verify_harmony_bridge(
    pca_embeddings: Optional[np.ndarray] = None,
    n_cells: int = 200,
    n_pcs: int = 15,
    random_state: int = 0,
    corr_min: float = 0.98,
    relative_shift_max: float = 0.25,
) -> dict:
    """No-real-batch-effect, theta=0 identity self-check for the RunHarmony rpy2 bridge.

    Splits cells into two batch labels drawn from the *same* distribution (no
    injected offset) and runs with ``theta=0``. R's ``RunHarmony`` internally
    builds a design matrix via ``contrasts<-`` over ``vars_use`` and errors
    with "contrasts can be applied only to factors with 2 or more levels" for
    a genuinely single-level batch column, so two labels are required -- but
    since neither the labels' underlying distribution nor theta gives Harmony
    any real signal or pressure to correct, the returned embedding should
    still be very close to the input: per-column correlation near 1 and a
    small per-cell shift relative to the input's norm. A low correlation or
    large shift means the bridge (matrix construction, orientation
    resolution, or R round-trip) is scrambling or mis-orienting data rather
    than passing it through -- exactly the failure mode a full run with real
    batches/theta would hide behind a plausible-looking, converged, wrong
    result.

    NOT called automatically inside :func:`run_harmony_rpy2` -- it's a
    dedicated identity run, too slow to repeat on every call against
    full-size data. Call it explicitly from the smoke test or a notebook
    cell before trusting a real integration run.

    Returns a dict with ``mean_column_correlation`` and ``mean_relative_shift``.
    Raises ``ValueError`` if either check fails.
    """
    rng = np.random.default_rng(random_state)
    if pca_embeddings is None:
        pca_embeddings = rng.normal(size=(n_cells, n_pcs)).astype(np.float64)
    else:
        pca_embeddings = np.asarray(pca_embeddings, dtype=np.float64)
        n_cells, n_pcs = pca_embeddings.shape

    adata = ad.AnnData(X=np.zeros((n_cells, 1), dtype=np.float32))
    adata.obsm['X_pca'] = pca_embeddings
    # Two labels (R needs >=2 factor levels), same distribution -- no real
    # batch effect, so theta=0 should leave the embedding essentially alone.
    adata.obs['batch'] = (['batch_1'] * (n_cells // 2)) + (['batch_2'] * (n_cells - n_cells // 2))
    adata.obs_names = [f'cell_{i}' for i in range(n_cells)]

    result = run_harmony_rpy2(
        adata, batch_key='batch', n_pcs=n_pcs, theta=0.0, lambda_val=1.0,
        max_iter=10, tau=0,
    )
    output = np.asarray(result.obsm['X_harmony'])

    col_corrs = _per_column_correlation(pca_embeddings, output)
    mean_corr = float(np.nanmean(col_corrs))

    shift = np.linalg.norm(output - pca_embeddings, axis=1)
    input_norm = np.linalg.norm(pca_embeddings, axis=1)
    relative_shift = float(np.mean(shift / np.clip(input_norm, 1e-8, None)))

    print(f'  [_verify_harmony_bridge] mean per-column correlation: {mean_corr:.4f} '
          f'(threshold: >= {corr_min})')
    print(f'  [_verify_harmony_bridge] mean relative per-cell shift: {relative_shift:.4f} '
          f'(threshold: <= {relative_shift_max})')

    if mean_corr < corr_min:
        raise ValueError(
            f'_verify_harmony_bridge FAILED: no-real-batch-effect theta=0 identity check gives mean '
            f'per-column correlation {mean_corr:.4f} < {corr_min}. The bridge is scrambling or '
            f'mis-orienting the PCA matrix.'
        )
    if relative_shift > relative_shift_max:
        raise ValueError(
            f'_verify_harmony_bridge FAILED: no-real-batch-effect theta=0 identity check gives mean '
            f'relative per-cell shift {relative_shift:.4f} > {relative_shift_max}, larger than '
            f'expected for a no-op correction. The bridge is scrambling or mis-orienting the '
            f'PCA matrix.'
        )

    return {'mean_column_correlation': mean_corr, 'mean_relative_shift': relative_shift}
