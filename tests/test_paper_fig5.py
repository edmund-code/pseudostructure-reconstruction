"""Synthetic checks for the Fig. 5 screen (analysis/paper/fig5_screen.py)."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis/paper'))
pytest.importorskip('pandas')


def test_mouse_quantile_maps_human_units_through_the_mice():
    from fig5_screen import mouse_quantile

    mean = np.array([0.1, 0.2, 0.3, 0.4, 0.15, 0.9])
    sample = np.array(['Ctrl1A2', 'Ctrl1A4', 'Ctrl1A2', 'Ctrl1A4', 'HUK1_COR1', 'HUK1_MED1'])
    x = mouse_quantile(mean, sample, np.ones(6, bool))
    assert np.allclose(x[:4], [0.125, 0.375, 0.625, 0.875])
    assert x[4] == pytest.approx(0.25) and x[5] == pytest.approx(0.875)


def test_end_contrast_reads_only_inside_the_section_range():
    from fig5_screen import GRID, end_contrast

    e = end_contrast(2.0 * GRID, 0.3, 0.7, (0.0, 1.0))
    s1, s3 = GRID[GRID <= 0.24].mean(), GRID[GRID >= 0.76].mean()
    assert e == pytest.approx(2.0 * (s3 - s1))
    assert np.isnan(end_contrast(2.0 * GRID, 0.3, 0.7, (0.0, 0.7)))      # no S3 core inside the range


def test_reversal_needs_both_species_and_is_blind_under_relabeling():
    from fig5_screen import SECTIONS, SPLITS, category_r

    # gene 0: mouse rises, human falls; gene 1: one human section disagrees; gene 2: the sections disagree within
    # each species, in the pattern of relabeling 1 (rows: Ctrl1A2, Ctrl1A4, COR1, MED1; columns: genes)
    E = dict(zip(SECTIONS, np.array([[1.0, 1.0, 1.0], [1.0, 1.0, -1.0], [-1.0, 1.0, 1.0], [-1.0, -1.0, -1.0]])))
    R = {m: v.copy() for m, v in E.items()}
    ok = np.ones(3, bool)
    calls = category_r(E, R, ok, *SPLITS['species'])
    assert calls.tolist() == [True, False, False]
    # gene 2 reverses between the relabeling-1 arms {COR1, Ctrl1A2} and {MED1, Ctrl1A4}
    assert category_r(E, R, ok, *SPLITS['relabeling 1']).tolist() == [False, False, True]
    R['Ctrl1A2'] = -R['Ctrl1A2']                                            # raw counts disagree with the curve
    assert not category_r(E, R, ok, *SPLITS['species'])[0]


def test_steeper_in_human_needs_every_pair():
    from fig5_screen import SPLITS, category_h

    E = {'Ctrl1A2': np.array([0.3, 0.3]), 'Ctrl1A4': np.array([0.4, 1.0]),
         'HUK1_COR1': np.array([1.5, 1.5]), 'HUK1_MED1': np.array([1.2, 1.2])}
    R = {m: v.copy() for m, v in E.items()}
    assert category_h(E, R, np.ones(2, bool), *SPLITS['species']).tolist() == [True, False]   # 1.2 - 1.0 < 0.5


def test_flattening_residual_removes_a_uniform_attenuation():
    from fig5_screen import flattening_residuals

    rng = np.random.default_rng(0)
    m = rng.normal(0, 1, 500)
    E = {'Ctrl1A2': m, 'Ctrl1A4': m + rng.normal(0, 0.01, 500), 'HUK1_COR1': 0.4 * m, 'HUK1_MED1': 0.4 * m}
    E['HUK1_COR1'][0] = E['HUK1_MED1'][0] = -0.4 * m[0]                    # one reversed gene
    d, pairs, lines = flattening_residuals(E, np.ones(500, bool), 'species')
    assert all(abs(s - 0.4) < 0.02 for _, s in lines)
    assert abs(d[0]) > 0.5 * abs(m[0]) and np.median(np.abs(d[1:])) < 0.02
