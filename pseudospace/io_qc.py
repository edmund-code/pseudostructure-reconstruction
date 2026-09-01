"""Sample loading and QC for mouse-only tubule AnnData objects.

Extracted verbatim from ``6_mouse_only_pseudospace.ipynb`` (Section 1). These
helpers take an AnnData plus explicit parameters and do not read notebook
globals, so they are reusable across notebooks 1-6.
"""
from __future__ import annotations

import os
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


def calculate_qc_metrics(adata: ad.AnnData, species: str = 'mouse') -> ad.AnnData:
    """QC metrics with species-appropriate mito/ribo regex. Verbatim from notebook 1 Part 3."""
    adata = adata.copy()

    if species == 'mouse':
        mito_pattern = '^mt-'
        ribo_pattern = '^Rpl|^Rps'
    else:
        mito_pattern = '^MT-'
        ribo_pattern = '^RPL|^RPS'

    adata.var['mt'] = adata.var_names.str.contains(mito_pattern, regex=True)
    adata.var['ribo'] = adata.var_names.str.contains(ribo_pattern, regex=True)

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt', 'ribo'],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    print(
        f'    Mito %: mean={adata.obs["pct_counts_mt"].mean():.2f}, '
        f'max={adata.obs["pct_counts_mt"].max():.2f}'
    )
    print(
        f'    Ribo %: mean={adata.obs["pct_counts_ribo"].mean():.2f}, '
        f'max={adata.obs["pct_counts_ribo"].max():.2f}'
    )

    return adata


def annotate_mito_ribo_mouse_symbols(adata: ad.AnnData) -> ad.AnnData:
    """
    Annotate mitochondrial and ribosomal genes in mouse-symbol space.

    These genes are KEPT in the final object but excluded from HVG/PCA/Harmony via
    exclude_from_harmony_hvg. Verbatim from notebook 1 Part 3.
    """
    adata = adata.copy()

    mito_pattern = '^mt-'
    ribo_patterns = ['^Rpl', '^Rps', '^Mrpl', '^Mrps']

    adata.var['mt'] = adata.var_names.str.contains(mito_pattern, regex=True)

    ribo = pd.Series(False, index=adata.var_names)
    for pattern in ribo_patterns:
        ribo = ribo | adata.var_names.str.contains(pattern, regex=True)

    adata.var['ribo'] = ribo.values
    adata.var['exclude_from_harmony_hvg'] = adata.var['mt'] | adata.var['ribo']

    print('\nGene annotations:')
    print(f'  Mito genes kept in final object: {int(adata.var["mt"].sum())}')
    print(f'  Ribo genes kept in final object: {int(adata.var["ribo"].sum())}')
    print(
        f'  Genes excluded only from Harmony HVG selection: '
        f'{int(adata.var["exclude_from_harmony_hvg"].sum())}'
    )

    return adata


def sample_name_from_path(path: Path) -> str:
    return path.name.replace('_tubule_by_gene_caleb.h5ad', '')


def load_and_process_mouse_samples(data_dir: Path, files: list[str]) -> dict[str, ad.AnnData]:
    """
    Load each mouse sample and annotate obs metadata. No ortholog conversion -- everything
    stays in native mouse gene symbols. Adapted from notebook 1 Part 3's
    load_and_process_samples with the human branch and ortholog conversion removed.
    """
    adatas = {}

    for f in files:
        name = sample_name_from_path(Path(f))
        print(f'\nProcessing: {name}')

        adata = sc.read_h5ad(os.path.join(data_dir, f))

        if 'counts' in adata.layers:
            adata.X = adata.layers['counts'].copy()

        adata = calculate_qc_metrics(adata, species='mouse')

        adata.obs['original_obs_name'] = adata.obs_names.astype(str).copy()
        adata.obs_names = pd.Index([f'{name}_{obs}' for obs in adata.obs_names])

        adata.obs['sample'] = name
        adata.obs['species'] = 'mouse'
        adata.obs['condition'] = 'Control' if 'Ctrl' in name else 'IR'

        print(f'    {adata.n_obs} tubules x {adata.n_vars} genes (native mouse gene space)')
        adatas[name] = adata

    return adatas


def combine_adatas(adatas: dict[str, ad.AnnData]) -> ad.AnnData:
    """Inner join across the 4 mouse samples (same reference -> effectively a union)."""
    adata = sc.concat(list(adatas.values()), join='inner')

    if not adata.obs_names.is_unique:
        adata.obs_names_make_unique()
    if not adata.var_names.is_unique:
        adata.var_names_make_unique()

    return adata
