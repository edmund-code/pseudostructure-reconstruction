"""Ask what mechanism could produce the human arm's shallower positional gradients.

The human proximal-tubule arm of this cohort carries systematically smaller gradient
amplitudes than the mouse arm, and the coordinate validation already records it
(``criteria_hvm.csv`` C12, "A3: no attenuation flag", not met). Three mechanisms could
do that: a blurred or mis-registered human axis, a per-gene difference in capture
between the two species' probe panels, or genuinely shallower human zonation.

This script tests the first one and nothing else. A positional smear is a convolution,
so it has two signatures that neither of the other mechanisms has: it can only *reduce*
an amplitude, never raise one, and it must attenuate a rough curve more than a smooth
one. Both are measurable here. The script therefore

  1. defines a clean gene set (converged fits, no separation, amplitudes in range, a
     mean-rate floor in both species), because a separated fit reports amplitudes in the
     tens of thousands of log units and would dominate every summary;
  2. reports the observed amplitude-ratio distribution and the share of genes that are
     *sharper* in human, including the share that are sharper in all four specimen
     pairings;
  3. convolves the mouse curves with Gaussians spanning the committed positional-noise
     range and reads the same two signatures off the simulation;
  4. measures the roughness dependence with the ratio's denominator taken from one mouse
     and the roughness from the other, so that the two do not share estimation noise;
  5. checks whether the flatness reproduces across the two human sections, stratified by
     how attenuated the gene is.

It fits nothing and writes nothing. Run it against a results root that holds the frozen
package pathway outputs:

    python analysis/scripts/attenuation_mechanism_human_vs_mouse.py --results-root <dir>

A "ponytail": the fitted curves are splines with three interior knots, so both arms are
already smoothed by the basis and the simulation cannot probe widths below roughly one
knot spacing. The two signatures it relies on are a sign and a tail share, and neither
depends on resolving small widths; a sharper test would need the per-unit counts.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.stats import spearmanr

# The species gene-model cache, by its content-addressed name under the pathway workflow.
GENE_MODEL = (
    "pt_pathways_package_pseudospace/stage_cache/"
    "gene_model_human_vs_mouse_species__06ba95871e922d1a.npz"
)
GENE_STATS = "pt_pathways_package_pseudospace/tables/gene_statistics_human_vs_mouse.csv"
BINNED_RATES = "pt_pathways_package_pseudospace/tables/s10_mechanics_human_vs_mouse.npz"
BINNED_GENES = "pt_pathways_package_pseudospace/tables/s10_mechanics_genes_human_vs_mouse.csv"

MOUSE = ("Ctrl1A2", "Ctrl1A4")
HUMAN = ("HUK1_COR1", "HUK1_MED1")

# The committed human fold-refit positional noise, in axis units, plus a width above it.
SMEAR_WIDTHS = (0.05, 0.10, 0.163, 0.221, 0.250)


def _roots() -> Path:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", default=os.environ.get("PSEUDOSPACE_RESULTS_ROOT", "results"))
    args, _ = parser.parse_known_args()
    return Path(args.results_root)


def amplitude(curves: np.ndarray) -> np.ndarray:
    """P90 minus P10 of a log-rate curve: a spread that ignores the per-gene offset."""
    return np.percentile(curves, 90, axis=1) - np.percentile(curves, 10, axis=1)


def roughness(curves: np.ndarray, amp: np.ndarray) -> np.ndarray:
    """Total variation divided by range: 1.0 for a monotone curve, larger for a wiggly one."""
    return np.abs(np.diff(curves, axis=1)).sum(axis=1) / np.maximum(amp, 1e-9)


def main() -> None:
    root = _roots()
    model = np.load(root / GENE_MODEL)
    stats = pd.read_csv(root / GENE_STATS, header=[0, 1], index_col=0)
    stats = stats[stats.index.notna()]

    # Guard: the curve matrix and the statistics table must be the same genes in the same
    # order. Without this every number below could be silently off by a row.
    drift = np.abs(stats[("condition", "T_spatial")].to_numpy() - model["stat|T_spatial"]).max()
    assert drift < 1e-9, f"gene order disagrees between the cache and the table ({drift:.3g})"

    genes = np.asarray(stats.index)
    grid = model["grid"]
    dx = float(grid[1] - grid[0])
    mouse_curves, human_curves = model["reference"], model["case"]
    per_specimen = {name: model[f"curve|{name}"] for name in MOUSE + HUMAN}

    separated = np.zeros(len(genes), dtype=bool)
    converged = np.ones(len(genes), dtype=bool)
    for name in MOUSE + HUMAN:
        separated |= model[f"separated|{name}"].astype(bool)
        converged &= model[f"converged|{name}"].astype(bool)

    # A mean-rate floor from the model-free binned rates, so the clean set does not depend
    # on the fit it is used to judge.
    binned = np.load(root / BINNED_RATES, allow_pickle=True)
    binned_genes = pd.read_csv(root / BINNED_GENES)["gene"].tolist()
    assert binned_genes == list(genes), "the binned-rate table is a different gene order"
    specimens = [str(s) for s in binned["obs_specimens"]]
    rates = binned["obs_rate_per10k"]
    rate_of = lambda names: np.vstack(
        [rates[[i for i, s in enumerate(specimens) if s == n]].mean(axis=0) for n in names]
    ).mean(axis=0)

    amp_m, amp_h = amplitude(mouse_curves), amplitude(human_curves)
    amp = {name: amplitude(curves) for name, curves in per_specimen.items()}

    clean = (
        converged
        & ~separated
        & (amp_m > 0.3)
        & (amp_m < 10)
        & (amp_h < 10)
        & (rate_of(MOUSE) > 0.5)
        & (rate_of(HUMAN) > 0.5)
    )
    print(f"separated in at least one specimen: {separated.sum()}")
    print(f"clean genes: {clean.sum()} of {len(genes)}")

    ratio = np.log2(amp_h / amp_m)
    print("\n1. observed amplitude ratio, human over mouse, clean genes")
    for q in (5, 25, 50, 75, 95):
        value = np.percentile(ratio[clean], q)
        print(f"   p{q:<3d} log2 {value:+.3f}   ratio {2 ** value:.3f}")
    sharper_any = ratio[clean] > 0
    both = lambda f: f(amp[HUMAN[0]], amp[MOUSE[0]]) & f(amp[HUMAN[0]], amp[MOUSE[1]]) & f(
        amp[HUMAN[1]], amp[MOUSE[0]]
    ) & f(amp[HUMAN[1]], amp[MOUSE[1]])
    sharper_all = clean & both(np.greater)
    flatter_all = clean & both(np.less)
    print(f"   sharper in human at all: {sharper_any.mean():.3f}")
    print(f"   sharper in all four specimen pairings: {sharper_all.sum()} ({sharper_all.sum() / clean.sum():.3f})")
    print(f"   flatter in all four: {flatter_all.sum()} ({flatter_all.sum() / clean.sum():.3f})")

    print("\n2. a positional smear of the mouse curves, for comparison")
    print("   width   median ratio   frac sharper   rho(log2 ratio, roughness)")
    rough_m = roughness(mouse_curves, amp_m)
    for width in SMEAR_WIDTHS:
        smeared = gaussian_filter1d(mouse_curves, width / dx, axis=1, mode="nearest")
        simulated = np.log2(amplitude(smeared) / amp_m)
        print(
            f"   {width:.3f}   {2 ** np.median(simulated[clean]):12.3f}   "
            f"{(simulated[clean] > 0).mean():12.3f}   {spearmanr(simulated[clean], rough_m[clean]).statistic:+.3f}"
        )
    print(
        f"   OBSERVED {2 ** np.median(ratio[clean]):11.3f}   {sharper_any.mean():12.3f}   "
        f"{spearmanr(ratio[clean], rough_m[clean]).statistic:+.3f}"
    )
    print("   a smear can only lower an amplitude, and must attenuate rough curves more,")
    print("   so a positive roughness dependence and a large sharper-in-human share exclude it.")

    # The ratio's denominator shares estimation noise with any mouse-side stratifier, which
    # would induce a negative correlation on its own. Take the denominator from one mouse and
    # the stratifier from the other.
    split = np.log2(amp_h / amp[MOUSE[1]])
    rough_split = roughness(per_specimen[MOUSE[0]], amp[MOUSE[0]])
    print("\n3. the same dependences with independent noise in the stratifier")
    print(f"   rho(log2 ratio, roughness of the other mouse) = {spearmanr(split[clean], rough_split[clean]).statistic:+.3f}")
    print(f"   rho(log2 ratio, amplitude of the other mouse)  = {spearmanr(split[clean], amp[MOUSE[0]][clean]).statistic:+.3f}")

    print("\n4. does the flatness reproduce across the two human sections?")
    print("   stratum            n     rho(COR1, MED1)   median |log2 COR1/MED1|   both below both mice")
    edges = [(0, 25, "most attenuated"), (25, 50, "p25-50"), (50, 75, "p50-75"), (75, 100, "least attenuated")]
    for low, high, name in edges:
        band = clean & (ratio >= np.percentile(ratio[clean], low)) & (ratio < np.percentile(ratio[clean], high))
        agree = spearmanr(amp[HUMAN[0]][band], amp[HUMAN[1]][band]).statistic
        spread = np.median(np.abs(np.log2(amp[HUMAN[0]][band] / amp[HUMAN[1]][band])))
        print(f"   {name:<18s} {band.sum():4d}   {agree:+.3f}            {spread:.3f}                     {flatter_all[band].mean():.3f}")

    order = np.argsort(-np.where(sharper_all, amp_h / amp_m, 0))
    print("\n5. the genes sharper in human in all four pairings, top 20")
    for i in order[: 20]:
        print(
            f"   {genes[i]:<12s} mouse {amp_m[i]:.2f} -> human {amp_h[i]:.2f}  x{amp_h[i] / amp_m[i]:.2f}"
            f"   ({MOUSE[0]} {amp[MOUSE[0]][i]:.2f}, {MOUSE[1]} {amp[MOUSE[1]][i]:.2f},"
            f" {HUMAN[0]} {amp[HUMAN[0]][i]:.2f}, {HUMAN[1]} {amp[HUMAN[1]][i]:.2f})"
        )


if __name__ == "__main__":
    main()
