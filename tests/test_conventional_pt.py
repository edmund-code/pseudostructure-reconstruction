import numpy as np
import pandas as pd
import pytest

pytest.importorskip('pydeseq2')

from pseudospace.conventional_pt import deseq2_contrast


def test_signed_pseudobulk_contrast_and_count_guard():
    rng = np.random.default_rng(7)
    names = ['mouse_1', 'mouse_2', 'human_1', 'human_2']
    counts = pd.DataFrame(rng.negative_binomial(12, .25, (4, 120)), index=names,
                          columns=[f'g{i}' for i in range(120)])
    counts.loc[names[2:], counts.columns[:20]] *= 5
    metadata = pd.DataFrame({'species': ['mouse', 'mouse', 'human', 'human']}, index=names)
    result = deseq2_contrast(counts, metadata, design='~species',
                             contrast=('species', 'human', 'mouse'))
    assert result.set_index('gene').loc[counts.columns[:20], 'stat'].median() > 0
    assert result.set_index('gene').loc[counts.columns[:20], 'log2FoldChange'].median() > 1
    assert result.gene.is_unique and len(result) == counts.shape[1]
    with pytest.raises(ValueError, match='integer raw counts'):
        deseq2_contrast(counts.astype(float) + .5, metadata, design='~species',
                        contrast=('species', 'human', 'mouse'))
