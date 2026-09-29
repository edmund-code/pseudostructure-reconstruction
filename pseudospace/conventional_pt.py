"""Specimen-level PT pseudobulk contrasts used by notebook 12."""

from __future__ import annotations

import numpy as np
import pandas as pd


def deseq2_contrast(counts: pd.DataFrame, metadata: pd.DataFrame, *, design: str,
                    contrast: tuple[str, str, str]) -> pd.DataFrame:
    """Fit a signed DESeq2 contrast to integer specimen-level pseudobulk counts.

    The caller defines independent units and the gene filter. In particular, this
    routine does not turn sections from one donor into biological replicates.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    if (counts.empty or counts.shape[0] < 4 or not counts.index.is_unique
            or not counts.columns.is_unique or not counts.index.equals(metadata.index)):
        raise ValueError('Need aligned pseudobulk counts and metadata with unique samples and genes.')
    values = counts.to_numpy()
    if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
        raise ValueError('DESeq2 requires finite, nonnegative integer raw counts.')
    if (counts.sum(axis=0) == 0).any():
        raise ValueError('Remove zero-total genes before DESeq2 fitting.')
    factor, numerator, denominator = contrast
    if factor not in metadata or set(metadata[factor].astype(str)) != {numerator, denominator}:
        raise ValueError('Contrast levels must both occur in the metadata.')

    dds = DeseqDataSet(counts=counts.astype(np.int64), metadata=metadata, design=design,
                       refit_cooks=False, quiet=True, n_cpus=2)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=list(contrast), quiet=True, n_cpus=2)
    stats.summary()
    result = stats.results_df.rename_axis('gene').reset_index()
    result['contrast'] = f'{numerator}_vs_{denominator}'
    return result
