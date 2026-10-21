"""Ask whether the sharper-in-human gene set is measuring distal-tubule admixture.

A set of genes carries steeper axial gradients in the human arm than in the mouse arm
(``attenuation_mechanism_human_vs_mouse.py`` reports it: 121 genes sharper in human in
all four specimen pairings). The cheapest attack on that set does not touch the axis at
all — it attacks membership. If a human proximal-tubule unit collects more neighbouring
distal tubule than a mouse unit does, and if that admixture is itself distributed along
the axis, then a distal-identity gene acquires an apparent proximal gradient that is
steeper in human for a reason that has nothing to do with proximal-tubule biology. The
set's own head is the reason to take this seriously: it contains Umod, the most abundant
transcript of the thick ascending limb.

The test has to be a measurement rather than a citation, because an external segment
atlas would answer a different question (where a gene sits in *that* cohort) and cannot
see our units. So this script asks the question with our own objects only:

  1. it builds a per-specimen, per-bin distal-identity score from the model-free binned
     rates, using thick-limb, distal-convoluted and collecting-duct markers and nothing
     from the fitted curves;
  2. it reports whether that score is itself axially structured, and whether it is more
     structured in the human arm — which is the precondition the attack needs;
  3. it correlates every gene's binned profile with the score, so that a gene whose
     gradient is admixture rides the score and a gene whose gradient is proximal does not;
  4. it re-counts the sharper-in-human set after dropping the riders, and names them.

It fits nothing and writes nothing.

    python analysis/scripts/sharpening_set_admixture_human_vs_mouse.py --results-root <dir>

A "ponytail": a per-bin score cannot separate a gene expressed in admixed distal tissue
from a proximal gene that happens to co-vary with where the admixture sits, because both
ride the same 20 numbers. The per-unit version of this test is strictly better and needs
the unit-level matrices; what is here is the cheap screen that decides whether the
expensive one is necessary.
"""

from __future__ import annotations

import argparse
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ConstantInputWarning, spearmanr

GENE_MODEL = (
    "pt_pathways_package_pseudospace/stage_cache/"
    "gene_model_human_vs_mouse_species__06ba95871e922d1a.npz"
)
GENE_STATS = "pt_pathways_package_pseudospace/tables/gene_statistics_human_vs_mouse.csv"
BINNED_RATES = "pt_pathways_package_pseudospace/tables/s10_mechanics_human_vs_mouse.npz"
BINNED_GENES = "pt_pathways_package_pseudospace/tables/s10_mechanics_genes_human_vs_mouse.csv"

MOUSE = ("Ctrl1A2", "Ctrl1A4")
HUMAN = ("HUK1_COR1", "HUK1_MED1")

# Identity genes for the tubule segments that physically abut the proximal tubule. None of
# these is a proximal-tubule gene in either species, so axial structure in their pooled
# signal inside proximal-tubule units is admixture by construction.
DISTAL_MARKERS = (
    "Umod", "Slc12a1", "Kcnj1",            # thick ascending limb
    "Slc12a3", "Pvalb", "Calb1",            # distal convoluted tubule
    "Aqp2", "Aqp3", "Scnn1g",               # collecting duct, principal
    "Atp6v1g3", "Atp6v0d2", "Foxi1", "Slc26a4",  # collecting duct, intercalated
)
# A proximal-tubule identity panel, as the contrast: the score built from these should be
# axially structured in both arms, because that is the axis the units actually span.
PROXIMAL_MARKERS = ("Slc4a4", "Lrp2", "Gatm", "Miox", "Slc22a6")

# Correlation with the distal score, required in both sections of a species, to call a gene
# a rider. The sign is not cosmetic: admixture ADDS a distal gene's signal where the distal
# tissue is, so the admixture hypothesis predicts a positive correlation. A gene strongly
# anti-correlated with the score is high where the admixture is absent, which is what a
# genuinely proximal gene does; it is counted separately rather than discarded.
RIDER_RHO = 0.5


def _results_root() -> Path:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-root", default=os.environ.get("PSEUDOSPACE_RESULTS_ROOT", "results")
    )
    args, _ = parser.parse_known_args()
    return Path(args.results_root)


def amplitude(curves: np.ndarray) -> np.ndarray:
    """P90 minus P10 of a log-rate curve: a spread that ignores the per-gene offset."""
    return np.percentile(curves, 90, axis=1) - np.percentile(curves, 10, axis=1)


def main() -> None:
    # A gene that is zero in every bin of a specimen has no defined rank correlation; it
    # comes back nan, which fails every threshold below, which is the right answer.
    warnings.simplefilter("ignore", ConstantInputWarning)
    root = _results_root()
    model = np.load(root / GENE_MODEL)
    stats = pd.read_csv(root / GENE_STATS, header=[0, 1], index_col=0)
    stats = stats[stats.index.notna()]

    drift = np.abs(stats[("condition", "T_spatial")].to_numpy() - model["stat|T_spatial"]).max()
    assert drift < 1e-9, f"gene order disagrees between the cache and the table ({drift:.3g})"

    genes = np.asarray(stats.index)
    index = {g: i for i, g in enumerate(genes)}

    binned = np.load(root / BINNED_RATES, allow_pickle=True)
    assert pd.read_csv(root / BINNED_GENES)["gene"].tolist() == list(genes), "binned gene order"
    obs_specimens = np.asarray([str(s) for s in binned["obs_specimens"]])
    obs_bins = np.asarray(binned["obs_bins"])
    rates = np.asarray(binned["obs_rate_per10k"])          # (80, n_genes)
    n_units = np.asarray(binned["obs_n"])

    # Per specimen: a (20, n_genes) profile, bins in axis order.
    profile, units = {}, {}
    for name in MOUSE + HUMAN:
        rows = np.where(obs_specimens == name)[0]
        order = rows[np.argsort(obs_bins[rows])]
        profile[name] = np.log1p(rates[order])
        units[name] = n_units[order]
    n_bins = profile[MOUSE[0]].shape[0]

    def panel_score(name: str, panel) -> tuple[np.ndarray, list[str]]:
        """Mean of per-bin z-scored log1p rates over the panel genes that are present."""
        used, columns = [], []
        for gene in panel:
            if gene not in index:
                continue
            values = profile[name][:, index[gene]]
            if rates[obs_specimens == name][:, index[gene]].mean() < 0.1 or values.std() < 1e-9:
                continue
            used.append((values - values.mean()) / values.std())
            columns.append(gene)
        return np.mean(used, axis=0), columns

    print("1. is the distal-identity signal inside the units axially structured at all?")
    print("   specimen     markers used   rho(score, bin)   distal mass per 10k   proximal rho")
    distal, proximal = {}, {}
    for name in MOUSE + HUMAN:
        distal[name], used_d = panel_score(name, DISTAL_MARKERS)
        proximal[name], _ = panel_score(name, PROXIMAL_MARKERS)
        mass = rates[obs_specimens == name][:, [index[g] for g in used_d]].sum(axis=1).mean()
        print(
            f"   {name:<12s} {len(used_d):>12d}   {spearmanr(distal[name], np.arange(n_bins)).statistic:+.3f}"
            f"            {mass:>12.2f}   {spearmanr(proximal[name], np.arange(n_bins)).statistic:+.3f}"
        )
    print(f"   distal panel kept in all four: {sorted(set.intersection(*[set(panel_score(n, DISTAL_MARKERS)[1]) for n in MOUSE + HUMAN]))}")
    print("   the attack needs this score to be both axial and more axial in human.")

    print("\n   agreement of the distal score between the two sections of a species")
    for arm, label in ((MOUSE, "mouse"), (HUMAN, "human")):
        print(f"   {label}: rho({arm[0]}, {arm[1]}) = {spearmanr(distal[arm[0]], distal[arm[1]]).statistic:+.3f}")

    # Per gene, per specimen: how much does its profile ride the distal score?
    rho_distal = {}
    for name in MOUSE + HUMAN:
        columns = profile[name]
        rho_distal[name] = np.array(
            [spearmanr(columns[:, j], distal[name]).statistic for j in range(len(genes))]
        )

    # The clean set and the sharper-in-human set, rebuilt exactly as the attenuation script
    # defines them so the two scripts are talking about the same genes.
    separated = np.zeros(len(genes), dtype=bool)
    converged = np.ones(len(genes), dtype=bool)
    for name in MOUSE + HUMAN:
        separated |= model[f"separated|{name}"].astype(bool)
        converged &= model[f"converged|{name}"].astype(bool)
    amp = {name: amplitude(model[f"curve|{name}"]) for name in MOUSE + HUMAN}
    amp_m, amp_h = amplitude(model["reference"]), amplitude(model["case"])
    rate_of = lambda names: np.vstack(
        [rates[obs_specimens == n].mean(axis=0) for n in names]
    ).mean(axis=0)
    clean = (
        converged & ~separated
        & (amp_m > 0.3) & (amp_m < 10) & (amp_h < 10)
        & (rate_of(MOUSE) > 0.5) & (rate_of(HUMAN) > 0.5)
    )
    both = lambda f: (
        f(amp[HUMAN[0]], amp[MOUSE[0]]) & f(amp[HUMAN[0]], amp[MOUSE[1]])
        & f(amp[HUMAN[1]], amp[MOUSE[0]]) & f(amp[HUMAN[1]], amp[MOUSE[1]])
    )
    sharper = clean & both(np.greater)
    print(f"\n2. the set under test: {sharper.sum()} genes sharper in human in all four pairings")

    up = lambda arm: (rho_distal[arm[0]] > RIDER_RHO) & (rho_distal[arm[1]] > RIDER_RHO)
    down = lambda arm: (rho_distal[arm[0]] < -RIDER_RHO) & (rho_distal[arm[1]] < -RIDER_RHO)
    rider_h, anti_h = up(HUMAN), down(HUMAN)
    rider_m, anti_m = up(MOUSE), down(MOUSE)
    print(f"   positive riders on the human distal score (rho > +{RIDER_RHO} in both sections)")
    print("   — the sign the admixture hypothesis predicts:")
    print(f"     among the sharper set   {(sharper & rider_h).sum():4d} of {sharper.sum()}   ({(sharper & rider_h).sum() / sharper.sum():.3f})")
    print(f"     among the clean set     {(clean & rider_h).sum():4d} of {clean.sum()}   ({(clean & rider_h).sum() / clean.sum():.3f})")
    print(f"   anti-riders (rho < -{RIDER_RHO}), high where the admixture is absent:")
    print(f"     among the sharper set   {(sharper & anti_h).sum():4d} ({(sharper & anti_h).sum() / sharper.sum():.3f})"
          f"   among the clean set {(clean & anti_h).sum():4d} ({(clean & anti_h).sum() / clean.sum():.3f})")
    print(f"   the same on the mouse score: positive {(sharper & rider_m).sum():4d} of {sharper.sum()}"
          f"   clean {(clean & rider_m).sum():4d} of {clean.sum()}")
    print(f"   surviving set after dropping positive human riders: {(sharper & ~rider_h).sum()}")
    print(f"   surviving set after dropping riders of either sign:  {(sharper & ~rider_h & ~anti_h).sum()}")

    print("\n3. does sharpening track the distal score across the clean set?")
    ratio = np.log2(amp_h / amp_m)
    for name in HUMAN + MOUSE:
        print(
            f"   rho(log2 human/mouse amplitude, rho_distal in {name}) = "
            f"{spearmanr(ratio[clean], rho_distal[name][clean]).statistic:+.3f}"
        )
    print("   a set created by admixture would show a strong positive dependence here.")

    print("\n4. the head of the set, with its admixture score")
    order = np.argsort(-np.where(sharper, amp_h / amp_m, 0))
    print("   gene          x sharper   rho_distal human   rho_distal mouse   rate h   rate m   rider")
    rate_h, rate_m = rate_of(HUMAN), rate_of(MOUSE)
    for i in order[:25]:
        rd_h = 0.5 * (rho_distal[HUMAN[0]][i] + rho_distal[HUMAN[1]][i])
        rd_m = 0.5 * (rho_distal[MOUSE[0]][i] + rho_distal[MOUSE[1]][i])
        mark = "RIDER+" if rider_h[i] else ("anti-" if anti_h[i] else "")
        print(
            f"   {genes[i]:<12s} {amp_h[i] / amp_m[i]:>9.2f}   {rd_h:>16.3f}   {rd_m:>16.3f}"
            f"   {rate_h[i]:>6.1f}   {rate_m[i]:>6.1f}   {mark}"
        )

    print("\n5. the head of the set that survives the admixture screen")
    survivor = sharper & ~rider_h
    s_order = np.argsort(-np.where(survivor, amp_h / amp_m, 0))
    print("   gene          x sharper   rho_distal human   rate h   rate m")
    for i in s_order[:20]:
        rd_h = 0.5 * (rho_distal[HUMAN[0]][i] + rho_distal[HUMAN[1]][i])
        print(
            f"   {genes[i]:<12s} {amp_h[i] / amp_m[i]:>9.2f}   {rd_h:>16.3f}"
            f"   {rate_h[i]:>6.1f}   {rate_m[i]:>6.1f}"
        )

    print("\n6. the genes the attack named, wherever they sit")
    for gene in ("Umod", "Ppp1r1a", "Mob1b", "Ddr1", "Dpp4", "Ggact", "Tmem184a", "Oxct1", "Pter", "Abhd14a", "Slc36a1"):
        if gene not in index:
            print(f"   {gene:<10s} not in the universe")
            continue
        i = index[gene]
        rd_h = 0.5 * (rho_distal[HUMAN[0]][i] + rho_distal[HUMAN[1]][i])
        rd_m = 0.5 * (rho_distal[MOUSE[0]][i] + rho_distal[MOUSE[1]][i])
        flags = []
        if sharper[i]:
            flags.append("sharper-in-human")
        if rider_h[i]:
            flags.append("POSITIVE RIDER")
        elif anti_h[i]:
            flags.append("anti-rider")
        if not clean[i]:
            flags.append("not clean")
        print(
            f"   {gene:<10s} x{amp_h[i] / amp_m[i]:>5.2f}   rho_distal h {rd_h:+.3f}  m {rd_m:+.3f}"
            f"   rate h {rate_h[i]:>7.2f}  m {rate_m[i]:>7.2f}   {', '.join(flags)}"
        )


if __name__ == "__main__":
    main()
