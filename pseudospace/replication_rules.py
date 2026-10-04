"""Pre-specified replication rules for positional species claims, and a grouped CSR reader.

Notebook 38 checks each R5 claim ("mouse falls from S1 to S3, human stays flat") against independent
segment-resolved datasets. The rules live here so they are fixed in code before any result is seen
and can be tested on synthetic values.

Conventions: every *arm* is a late-minus-early contrast (S3 - S1) on that dataset's own log scale.
A claim per arm is '+', '-', '0' (flat relative to the other species) or '' (not claimed).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

DEAD_ZONE = 0.1      # |contrast| below this is not a sign in either direction (log units)
FLAT_RATIO = 0.5     # a '0' arm must be at most half the other species' contrast, which must itself be signed
VERDICTS = ('headline', 'supporting', 'not replicated', 'not testable')


def csr_group_sums(data, indices, indptr, rows, codes, n_groups, n_cols, *, block_rows=20_000):
    """Sum selected CSR rows into groups without loading the whole matrix.

    ``data``, ``indices`` and ``indptr`` may be h5py datasets (read block by block). ``rows`` are row
    numbers and ``codes`` their group (0..n_groups-1). Returns ``(sums, totals, counts)``: the
    groups x columns sums, each group's sum over **all** columns, and the number of rows per group.
    """
    rows, codes = np.asarray(rows, np.int64), np.asarray(codes, np.int64)
    if rows.shape != codes.shape or rows.ndim != 1:
        raise ValueError('rows and codes must be 1-D and aligned.')
    if len(np.unique(rows)) != len(rows):
        raise ValueError('rows must be unique.')
    if len(codes) and (codes.min() < 0 or codes.max() >= n_groups):
        raise ValueError('codes must lie in [0, n_groups).')
    ptr = np.asarray(indptr[:], np.int64)
    if len(rows) and (rows.min() < 0 or rows.max() >= len(ptr) - 1):
        raise ValueError('row out of range.')
    order = np.argsort(rows)
    rows, codes = rows[order], codes[order]
    sums = np.zeros(n_groups * n_cols)
    totals = np.zeros(n_groups)
    counts = np.bincount(codes, minlength=n_groups)
    n_rows = len(ptr) - 1
    for first in range(0, n_rows, block_rows):
        last = min(first + block_rows, n_rows)
        lo, hi = np.searchsorted(rows, [first, last])
        if lo == hi:
            continue
        block_data = np.asarray(data[ptr[first]:ptr[last]], float)
        block_index = np.asarray(indices[ptr[first]:ptr[last]], np.int64)
        for row, code in zip(rows[lo:hi], codes[lo:hi]):
            a, b = ptr[row] - ptr[first], ptr[row + 1] - ptr[first]
            values = block_data[a:b]
            np.add.at(sums, code * n_cols + block_index[a:b], values)
            totals[code] += values.sum()
    return sums.reshape(n_groups, n_cols), totals, counts


def arm_passes(claim, value, other=None, *, dead_zone=DEAD_ZONE, flat_ratio=FLAT_RATIO):
    """True/False for one claimed arm; None when unclaimed or not measurable (NaN).

    '+' and '-' need the sign with magnitude at least ``dead_zone``. '0' needs ``|value| <=
    flat_ratio * |other|`` where ``other`` (the other species' contrast) is itself at least
    ``dead_zone`` in magnitude; two flat species do not make a species difference.
    """
    if claim in ('', None) or value is None or not np.isfinite(value):
        return None
    if claim == '+':
        return bool(value >= dead_zone)
    if claim == '-':
        return bool(value <= -dead_zone)
    if claim == '0':
        if other is None or not np.isfinite(other):
            return None
        return bool(abs(other) >= dead_zone and abs(value) <= flat_ratio * abs(other))
    raise ValueError(f'Unknown claim {claim!r}')


def gene_verdict(tests, *, critical, probe_balanced, detection_ok):
    """Classify one gene from its named test results (True/False/None).

    - any failed critical test, or two or more failed tests -> 'not replicated'
    - no test evaluated -> 'not testable'
    - no failure, every critical test evaluated, balanced probes, detection rule met -> 'headline'
    - otherwise (one non-critical failure, an unevaluated critical test, unbalanced probes or the
      detection rule unmet) -> 'supporting'
    """
    evaluated = {k: v for k, v in tests.items() if v is not None}
    if not evaluated:
        return 'not testable'
    failed = [k for k, v in evaluated.items() if v is False]
    if any(k in critical for k in failed) or len(failed) >= 2:
        return 'not replicated'
    critical_unevaluated = any(tests.get(k) is None for k in critical if k in tests)
    if not failed and not critical_unevaluated and probe_balanced and detection_ok:
        return 'headline'
    return 'supporting'


def story_verdict(verdicts, leads, probe_balanced, *, min_pass=0.6, can_lead=True):
    """Story stays as a lead only if every lead gene is 'headline' and at least ``min_pass`` of its
    probe-balanced members are 'headline' or 'supporting'. With the member rule alone it is
    'supporting'; otherwise 'not replicated'. ``verdicts`` and ``probe_balanced`` are gene-indexed.
    """
    verdicts, probe_balanced = pd.Series(verdicts), pd.Series(probe_balanced).reindex(pd.Series(verdicts).index)
    balanced = verdicts[probe_balanced.fillna(False).astype(bool)]
    if balanced.empty:
        return 'not testable', float('nan')
    share = float(balanced.isin(['headline', 'supporting']).mean())
    leads_ok = bool(leads) and all(verdicts.get(g) == 'headline' for g in leads)
    if share >= min_pass and leads_ok and can_lead:
        return 'lead', share
    if share >= min_pass:
        return 'supporting', share
    return 'not replicated', share
