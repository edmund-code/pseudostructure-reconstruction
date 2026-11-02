"""Shared definitions for the Fig. 3 validation analyses (segment order, robustness, baselines, depth).

Reads the package runs under --results-root (outside Git). Protocol: docs/paper/fig3-validation-experiments.md.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

RUNS = {'aki': 'package_pseudospace/runs/aki', 'hvm': 'package_pseudospace/runs/hvm'}
# the six sections of Fig. 3b: the four mouse sections from the control/AKI fit, the two human sections
# from the human-vs-mouse fit (the controls also sit in that fit; they are scored once, in the AKI fit)
SECTIONS = {'Ctrl1A2': 'aki', 'Ctrl1A4': 'aki', 'IR2A2': 'aki', 'IR2A4': 'aki', 'HUK1_COR1': 'hvm', 'HUK1_MED1': 'hvm'}

# fine stage-03 labels -> the seven displayed groups, in nephron order (thin limbs would sit between PTS3 and AL;
# neither run has any). Podocytes are not tubule and get no rank.
GROUP_ORDER = ['PTS1', 'PTS2', 'PTS3', 'TL', 'AL', 'DCT', 'CNT', 'CD']
GROUP = {'PT-S1': 'PTS1', 'PT-S2': 'PTS2', 'PT-S3': 'PTS3', 'DTL1': 'TL', 'DTL2': 'TL', 'DTL3': 'TL', 'ATL': 'TL',
         'mTAL': 'AL', 'cTAL': 'AL', 'Macula-densa': 'AL', 'DCT1': 'DCT', 'DCT2': 'DCT', 'CNT': 'CNT',
         'CCD': 'CD', 'OMCD': 'CD', 'IMCD': 'CD'}
RANK = {g: i for i, g in enumerate(GROUP_ORDER)}


def load_positions(results_root, run: str) -> pd.DataFrame:
    """One row per unit on the curve: sample, fine label, displayed group, its rank, and t (pseudostructure)."""
    d = pd.read_parquet(Path(results_root) / RUNS[run] / '04_pseudospace/pseudospace.parquet')
    d = d.set_index('obs_name') if 'obs_name' in d.columns else d
    out = pd.DataFrame({'sample': d['sample'].astype(str), 'segment_class': d['segment_class'].astype(str),
                        't': d['pseudospace'].astype(float)}, index=d.index)
    if 'dpt_pseudotime' in d:
        out['dpt'] = d['dpt_pseudotime'].astype(float)
    out['group'] = out.segment_class.map(GROUP)
    out['rank'] = out.group.map(RANK)
    return out


def section_rho(coord: pd.Series, table: pd.DataFrame, sections=None) -> pd.DataFrame:
    """Unit-level Spearman rho of a coordinate against the nephron-order rank of each unit's label, per section."""
    from scipy.stats import spearmanr

    rows = []
    for s in sections or sorted(table['sample'].unique()):
        k = (table['sample'] == s) & table['rank'].notna() & np.isfinite(coord.reindex(table.index))
        rows.append({'sample': s, 'n': int(k.sum()),
                     'rho': float(spearmanr(coord.reindex(table.index)[k], table.loc[k, 'rank']).statistic)})
    return pd.DataFrame(rows)
