"""Specimen-level summaries: balanced curves, pseudobulk profiles and coordinate agreement."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pseudospace.specimen import coordinate_agreement, pseudobulk_profiles, specimen_balanced_curves


def test_balanced_curves_give_every_specimen_equal_weight():
    big = np.array([[10.0, 10.0, 10.0]])
    small = np.array([[0.0, 0.0, 0.0]])
    balanced = specimen_balanced_curves({'big': big, 'small': small})
    assert np.allclose(balanced, 5.0), 'a 1000x larger specimen must not dominate'


def test_balanced_curves_ignore_specimens_without_support():
    full = np.array([[1.0, 1.0, np.nan]])
    partial = np.array([[np.nan, 3.0, np.nan]])
    balanced = specimen_balanced_curves({'full': full, 'partial': partial})
    assert balanced.shape == (1, 3)                 # (features, grid), as elsewhere
    assert balanced[0, 0] == pytest.approx(1.0)
    assert balanced[0, 1] == pytest.approx(2.0)
    assert np.isnan(balanced[0, 2])


def test_pseudobulk_sums_counts_per_specimen_and_bin():
    counts = np.array([
        [1.0, 10.0], [2.0, 20.0], [4.0, 40.0],
        [8.0, 80.0], [16.0, 160.0], [32.0, 320.0],
    ])
    specimen = ['A', 'A', 'A', 'B', 'B', 'B']
    pseudospace = [0.05, 0.06, 0.07, 0.85, 0.90, 0.95]
    frame = pseudobulk_profiles(counts, specimen, pseudospace, n_bins=10, gene_names=['g1', 'g2'],
                                min_structures_per_bin=2)
    first = frame[(frame['specimen'] == 'A') & (frame['pseudospace_bin'] == 0)].iloc[0]
    assert first['n_structures'] == 3
    assert first['count_g1'] == pytest.approx(7.0)
    assert first['count_g2'] == pytest.approx(70.0)
    second = frame[(frame['specimen'] == 'B') & (frame['pseudospace_bin'] == 9)].iloc[0]
    assert second['n_structures'] == 2
    assert second['count_g1'] == pytest.approx(48.0)
    assert second['count_g2'] == pytest.approx(480.0)
    assert len(frame) == 2, 'bins below the minimum structure count are not returned'


def test_pseudobulk_drops_bins_below_the_structure_minimum():
    counts = np.ones((4, 1))
    frame = pseudobulk_profiles(counts, ['A'] * 4, [0.05, 0.06, 0.07, 0.95], n_bins=10,
                                min_structures_per_bin=3)
    assert len(frame) == 1 and frame.iloc[0]['n_structures'] == 3


def test_coordinate_agreement_reports_rank_concordance():
    reference = np.linspace(0, 1, 20)
    same = coordinate_agreement(reference, reference)
    assert same['spearman'] == pytest.approx(1.0)
    assert same['concordant_pairs'] == pytest.approx(1.0)
    flipped = coordinate_agreement(reference, -reference)
    assert flipped['spearman'] == pytest.approx(-1.0)
    assert flipped['concordant_pairs'] == pytest.approx(0.0)
    with_nan = coordinate_agreement(reference, np.where(reference < 0.5, np.nan, reference))
    assert with_nan['n_compared'] == 10
