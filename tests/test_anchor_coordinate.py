import numpy as np
import pandas as pd
import pytest

from pseudospace.anchor_coordinate import (bulk_segment_profiles, census_segment_profiles,
                                           conserved_anchor_mask, positive_control_masks,
                                           pt_specific_mask, segment_contrast)


def test_bulk_profiles_average_replicates_on_log_scale():
    table = pd.DataFrame({'x_PTS1_1': [1, 3], 'x_PTS1_2': [3, 3], 'x_PTS2_1': [7, 0],
                          'x_PTS3_1': [15, 1]}, index=['a', 'b'])
    frame = bulk_segment_profiles(table, '_PT{seg}_')
    assert frame.loc['a', 'S1'] == pytest.approx((1 + 2) / 2)
    assert frame.loc['a', 'S3'] == pytest.approx(4)
    assert frame.loc['a', 'early'] == pytest.approx((1.5 + 3) / 2)


def _census():
    meta = pd.DataFrame({
        'donor_id': ['d1', 'd1', 'd2', 'd2', 'd3'], 'sex': ['male', 'male', 'female', 'female', 'male'],
        'cell_type': ['pct', 's3', 'pct', 's3', 'pct'], 'development_stage': ['adult'] * 5,
        'n_cells': [100, 100, 100, 100, 100], 'total_counts_all_genes': [1e6, 1e6, 2e6, 2e6, 1e6]},
        index=list('abcde'))
    counts = pd.DataFrame({'g1': [100, 300, 200, 600, 50], 'g2': [10, 10, 20, 20, 5]}, index=list('abcde'))
    return counts, meta


def test_census_profiles_keep_complete_donors_and_normalize_by_all_gene_totals():
    counts, meta = _census()
    frame = census_segment_profiles(counts, meta, {'pct': 'early', 's3': 'S3'})
    assert set(frame.columns.get_level_values('donor')) == {'d1', 'd2'}  # d3 lacks S3
    assert frame[('d1', 'male', 'early')]['g1'] == pytest.approx(np.log2(100 + 1))
    contrast, level = segment_contrast(frame)
    assert contrast['g1'] == pytest.approx(np.mean([np.log2(301 / 101), np.log2(301 / 101)]))
    male, _ = segment_contrast(frame, sex='male')
    assert male['g2'] == pytest.approx(0)


def test_anchor_mask_needs_same_sign_magnitude_and_expression_everywhere():
    genes = ['up', 'down', 'mixed', 'weak', 'quiet', 'ineligible']
    contrasts = pd.DataFrame({'h': [1, -2, 1, .2, 2, 3], 'm': [3, -1, -1, 2, 2, 3]}, index=genes)
    levels = pd.DataFrame({'h': [2, 2, 2, 2, .5, 2], 'm': [2, 2, 2, 2, 2, 2]}, index=genes)
    eligible = [True, True, True, True, True, False]
    mask = conserved_anchor_mask(contrasts, levels, eligible, threshold=.5)
    assert mask[mask].index.tolist() == ['up', 'down']


def test_pt_specific_mask_requires_every_species():
    mouse = pd.DataFrame({'PT': [1, 1], 'AL': [1.5, 5]}, index=['a', 'b'])
    human = pd.DataFrame({'PT': [1, 1], 'AL': [3, 1]}, index=['a', 'b'])
    mask = pt_specific_mask({'mouse': mouse, 'human': human})
    assert mask.tolist() == [False, False]
    assert pt_specific_mask({'mouse': mouse}).tolist() == [True, False]


def test_positive_controls_split_opposite_and_flat():
    idx = ['opp', 'flat', 'same', 'disagree']
    human = pd.Series([-1, .1, 1, -1], index=idx)
    mouse = pd.Series([1, 2, 1, 1], index=idx)
    check = pd.Series([2, 3, 1, -1], index=idx)
    level = pd.Series([2, 2, 2, 2], index=idx)
    opposite, flat = positive_control_masks(human, mouse, check, level, [True] * 4)
    assert opposite[opposite].index.tolist() == ['opp']
    assert flat[flat].index.tolist() == ['flat']
