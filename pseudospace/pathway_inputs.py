"""Recover the full ortholog panel without changing notebook 03's PT coordinates."""
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .cross_species import map_human_to_mouse_space


def rebuild_pt_expression(pt_obs, source_paths, orthologs, *, target_sum=1e4):
    """Load raw counts for exactly the saved structures, in their original order.

    Keep every accepted pair for the availability audit. Normalize against all
    orthologs measured in EVERY specimen, before any expression/detection filter.
    Human symbol mapping reuses the same implementation as notebook 03.
    """
    columns = ['mouse_symbol', 'human_symbol']
    if (not set(columns + ['mapping_status']).issubset(orthologs) or orthologs.empty or
            orthologs[columns].isna().any().any() or
            any(orthologs[c].astype(str).str.upper().duplicated().any() for c in columns) or
            any(orthologs[c].astype(str).str.strip().eq('').any() for c in columns)):
        raise ValueError('Expected a nonempty accepted one-to-one ortholog map.')
    if (not {'sample', 'comparison_species', 'shared_pseudospace'}.issubset(pt_obs) or
            not pt_obs.index.is_unique or pt_obs.empty or not np.isfinite(target_sum) or target_sum <= 0):
        raise ValueError('Need uniquely identified PT observations and a positive normalization target.')
    if set(pt_obs['sample'].astype(str)) != set(source_paths):
        raise ValueError('Source files must match the saved specimens exactly.')
    if (not pt_obs.comparison_species.isin(['mouse', 'human']).all() or
            not np.isfinite(pt_obs.shared_pseudospace).all()):
        raise ValueError('Expected human/mouse species labels and a finite saved coordinate.')
    for path in source_paths.values():
        if not path.is_file():
            raise FileNotFoundError(f'Missing original specimen matrix: {path}')

    symbols = pd.Index(orthologs.mouse_symbol.astype(str), name='gene')
    availability, matrices, row_ids = {}, [], []
    for sample, path in source_paths.items():
        observed = pt_obs[pt_obs['sample'].astype(str).eq(sample)]
        if observed.comparison_species.nunique() != 1:
            raise ValueError('A specimen cannot contain multiple species.')
        # Notebook 03 prefixes the original structure identifier with the sample name.
        prefix = sample + '_'
        if not observed.index.astype(str).str.startswith(prefix).all():
            raise ValueError('Saved structure identifiers must have their specimen prefix.')
        original_ids = observed.index.astype(str).str[len(prefix):]
        native = ad.read_h5ad(path)
        if not native.obs_names.is_unique or not native.var_names.is_unique:
            raise ValueError(f'Duplicate source identifiers in {sample}.')
        missing = original_ids.difference(native.obs_names)
        if len(missing):
            raise ValueError(f'{sample}: {len(missing)} saved PT structures are missing from the source.')
        native = native[original_ids].copy()
        raw = sparse.csr_matrix(native.layers.get('counts', native.X), dtype=float)
        if (not np.isfinite(raw.data).all() or (raw.data < 0).any() or
                not np.allclose(raw.data, np.round(raw.data), rtol=0, atol=1e-8)):
            raise ValueError(f'{sample}: expected nonnegative raw integer counts.')
        if observed.comparison_species.iloc[0] == 'human':
            native = ad.AnnData(X=raw, var=pd.DataFrame(index=native.var_names))
            mapped, _ = map_human_to_mouse_space(native, orthologs)
            matrix = sparse.csr_matrix(mapped[:, symbols].X)
            measured = mapped.var.loc[symbols, 'measured_in_source_input'].to_numpy(bool)
        else:
            idx = native.var_names.get_indexer(symbols)
            measured = idx >= 0
            transform = sparse.csr_matrix((np.ones(measured.sum()),
                (idx[measured], np.flatnonzero(measured))), shape=(native.n_vars, len(symbols)))
            matrix = raw @ transform
        matrices.append(matrix)
        row_ids.extend(observed.index)
        availability['measured_' + sample] = measured

    order = pd.Index(row_ids).get_indexer(pt_obs.index)
    counts = sparse.vstack(matrices, format='csr')[order]
    var = pd.DataFrame(availability, index=symbols)
    var['human_symbol'] = orthologs.set_index('mouse_symbol').loc[symbols, 'human_symbol'].to_numpy()
    var['measured_in_both_inputs'] = np.logical_and.reduce(list(availability.values()))
    if not var.measured_in_both_inputs.any():
        raise ValueError('No orthologs are measured in every specimen.')
    library_size = np.asarray(counts[:, var.measured_in_both_inputs.to_numpy()].sum(axis=1)).ravel()
    if (library_size <= 0).any():
        raise ValueError('A retained PT structure has zero counts in the shared measured panel.')
    lognorm = counts.multiply((target_sum / library_size)[:, None]).tocsr()
    lognorm.data = np.log1p(lognorm.data)
    result = ad.AnnData(X=lognorm, obs=pt_obs.copy(), var=var)
    result.layers['counts'] = counts
    result.layers['lognorm'] = lognorm.copy()
    result.obs['ortholog_library_size'] = library_size
    result.uns['expression_reconstruction'] = {
        'source': 'original specimen counts; counts layer when present, otherwise X',
        'normalization': 'log1p(count / all-measured-ortholog library size * target_sum)',
        'target_sum': float(target_sum), 'n_accepted_pairs': len(symbols),
        'n_measured_in_every_specimen': int(var.measured_in_both_inputs.sum()),
        'coordinate': 'saved PT structures and shared_pseudospace; unchanged',
    }
    pd.testing.assert_frame_equal(result.obs[pt_obs.columns], pt_obs)
    return result
