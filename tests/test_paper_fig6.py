"""Synthetic checks for the Fig. 6 screen (analysis/paper/fig6_screen.py)."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis/paper'))
pytest.importorskip('pandas')


def test_control_quantile_maps_injured_units_through_the_controls():
    from fig6_screen import control_quantile

    mean = np.array([0.1, 0.2, 0.3, 0.4, 0.15, 0.35, 0.9])
    sample = np.array(['Ctrl1A2', 'Ctrl1A4', 'Ctrl1A2', 'Ctrl1A4', 'IR2A2', 'IR2A4', 'IR2A2'])
    x = control_quantile(mean, sample, np.ones(7, bool))
    assert np.allclose(x[:4], [0.125, 0.375, 0.625, 0.875])          # Fig. 4's rank scale
    assert x[4] == pytest.approx(0.25) and x[5] == pytest.approx(0.75)  # interpolated between control ranks
    assert x[6] == pytest.approx(0.875)                               # beyond the last control unit: clipped


def test_within_core_reads_a_slope_inside_each_segment():
    from fig6_screen import GRID, within_core

    w = within_core(2.0 * GRID, 0.3, 0.7)                             # slope 2 everywhere
    lengths = np.array([0.3 - 0.06, 0.4 - 0.12, 0.3 - 0.06])
    assert np.allclose(w, 2.0 * lengths * 2 / 3, atol=0.02)           # last-third mean minus first-third mean


def test_category_a_needs_the_condition_and_not_the_relabelings():
    from fig6_screen import MICE, category_a

    # gene 0: IR changes inside S3, controls do not -> called; gene 1: the A2 sections change in both conditions
    # (the relabeling-1 factor) -> not specific; gene 2: IR change too small
    W = {m: np.zeros((3, 3)) for m in MICE}
    W['IR2A2'][0, 2] = W['IR2A4'][0, 2] = 1.5
    W['IR2A2'][1, 2] = W['Ctrl1A2'][1, 2] = 1.5
    W['IR2A2'][2, 2] = W['IR2A4'][2, 2] = 0.6
    L = {m: v.copy() for m, v in W.items()}
    calls = category_a(W, L, np.ones((3, 3), bool))
    assert calls[0, 2] and not calls[1, 2] and not calls[2, 2] and calls.sum() == 1
    assert not category_a(W, L, np.zeros((3, 3), bool)).any()         # the support rule


def test_set_screen_ignores_unsupported_genes():
    import fig6_screen as f6

    rng = np.random.default_rng(0)
    g = 400
    gene_w = {}
    for split in f6.SPLITS:
        base = rng.normal(0, 0.1, (g, 3))
        if split == 'condition':
            base[:20, 1] += 1.0                                       # a set that changes inside S2
        base[25, :] = np.nan                                          # an unsupported gene
        gene_w[split] = base
        for k in range(4):
            gene_w[f'{split}|{k}'] = base + rng.normal(0, 0.01, (g, 3))
    members = {'up in S2': np.r_[np.arange(20), 25], 'random': np.arange(100, 130)}
    df = f6.set_screen(gene_w, members, np.random.default_rng(1))
    row = df[(df.split == 'condition') & (df.set == 'up in S2') & (df.segment == 2)].iloc[0]
    assert row.n == 20 and row.call
    assert not df[(df.split == 'condition') & (df.set == 'random')].call.any()
