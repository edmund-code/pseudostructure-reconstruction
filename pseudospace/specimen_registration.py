"""Specimen-level offsets for projection onto a frozen mean atlas."""
from __future__ import annotations

import numpy as np

from .canonical_mean_atlas import _distances, _grid, _matrix


def fit_specimen_offset(atlas, calibration_Y, restricted=True, variance_fraction=0.99,
                        penalty=4.0, iterations=5):
    """Fit a shrunk specimen offset, optionally excluding atlas-path directions.

    The source atlas curve is frozen. Calibration rows determine the offset and
    nearest curve points, but no projected query rows enter this fit.
    """
    mean = _matrix(atlas["mean"], "atlas mean")
    grid = _grid(atlas["grid"])
    if len(grid) != len(mean):
        raise ValueError("atlas grid must match the atlas mean")
    y = _matrix(calibration_Y, "calibration_Y")
    if y.shape[1] != mean.shape[1]:
        raise ValueError("calibration_Y feature count does not match atlas")
    if not isinstance(restricted, (bool, np.bool_)):
        raise ValueError("restricted must be boolean")
    variance_fraction = float(variance_fraction)
    if not np.isfinite(variance_fraction) or not 0 < variance_fraction <= 1:
        raise ValueError("variance_fraction must be finite and in (0, 1]")
    penalty = float(penalty)
    if not np.isfinite(penalty) or penalty < 0:
        raise ValueError("penalty must be finite and non-negative")
    if not isinstance(iterations, (int, np.integer)) or isinstance(iterations, (bool, np.bool_)) or iterations < 1:
        raise ValueError("iterations must be a positive integer")

    centered = mean - mean.mean(axis=0, keepdims=True)
    _, singular, vh = np.linalg.svd(centered, full_matrices=False)
    variation = singular * singular
    total = float(variation.sum())
    if total == 0:
        rank, retained = 0, 1.0
    else:
        cumulative = np.cumsum(variation) / total
        cumulative[-1] = 1.0
        rank = int(np.searchsorted(cumulative, variance_fraction, side="left") + 1)
        rank = min(rank, len(singular))
        retained = float(cumulative[rank - 1])
    basis = vh[:rank].T
    if not restricted:
        projector = np.eye(mean.shape[1])
    elif rank == mean.shape[1]:
        projector = np.zeros((mean.shape[1], mean.shape[1]))
    else:
        projector = np.eye(mean.shape[1]) - basis @ basis.T

    offset = np.zeros(mean.shape[1], dtype=float)
    history = []
    previous_positions = None
    for _ in range(iterations):
        positions = np.argmin(_distances(y - offset, {"mean": mean}), axis=1)
        residual = y - mean[positions]
        proposed = projector @ residual.mean(axis=0) / (1.0 + penalty)
        history.append({
            "offset_norm": float(np.linalg.norm(proposed)),
            "position_changes": int(np.count_nonzero(positions != previous_positions))
            if previous_positions is not None else 0,
        })
        offset = proposed
        previous_positions = positions
    return offset, {
        "rank": rank,
        "retained_fraction": retained,
        "offset_norm": float(np.linalg.norm(offset)),
        "history": history,
        "calibration_n": int(len(y)),
    }
