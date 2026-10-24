"""Compare unweighted and equal-weighted S1->S2 / S2->S3 transitions on the frozen PT coordinate.

Notebook 72 places each specimen's segment transitions with an unweighted logistic regression of label on
position (``logistic_transitions``). Its 0.5 point moves with the two labels' counts, and the human sections
carry 7-9 PT-S1 units per PT-S2 unit, so the human S1->S2 transitions land deep in S2 and HUK1_MED1's cross
its S2->S3 transition (pt-package-issues.md, issues 6 and 26). This script refits the same inputs with each
label weighted by 1 / its count, checks that the unweighted fit still reproduces 72's table, and prints both,
plus a mouse control in which S2 is thinned to the human ratio. It fits nothing else and writes nothing:

    python analysis/scripts/balanced_segment_transitions.py --results-root <dir>
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pseudospace.continuous_descriptors import logistic_transitions  # noqa: E402
from pseudospace.package_inputs import support_positions  # noqa: E402

FROZEN = "pt_package_coordinate_validation/tables/transitions.csv"
UNITS = "package_pseudospace/{ds}/pt_units.h5ad"
THIN_RATIOS = (1, 3, 7, 9)   # S1:S2 after thinning mouse S2; human is 5.7-9.1


def _roots() -> Path:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", default=os.environ.get("PSEUDOSPACE_RESULTS_ROOT", "results"))
    args, _ = parser.parse_known_args()
    return Path(args.results_root)


def inputs(results: Path, ds: str):
    """Notebook 72's positions, labels and specimens: common-support units, min-max scaled."""
    obs = ad.read_h5ad(results / UNITS.format(ds=ds), backed="r").obs
    specimen, segment = obs["sample"].astype(str).to_numpy(), obs.segment_class.astype(str).to_numpy()
    mask, position = support_positions(obs.pseudospace.to_numpy(float), specimen, segment)
    return position, segment[mask], specimen[mask]


def main() -> None:
    results = _roots()
    frozen = pd.read_csv(results / FROZEN)
    tables = []
    for ds in ("hvm", "aki"):
        position, segment, specimen = inputs(results, ds)
        plain = logistic_transitions(position, segment, specimen)
        check = plain.merge(frozen[frozen.dataset == ds], on=["specimen", "transition"])
        assert len(check) == 8 and np.allclose(check.position_x, check.position_y, atol=1e-6), "inputs differ from 72"
        fair = logistic_transitions(position, segment, specimen, balanced=True)
        counts = pd.crosstab(specimen, segment)
        table = check.rename(columns={"position_y": "unweighted (72)"}).drop(columns="position_x").merge(
            fair.rename(columns={"position": "equal-weighted"}), on=["specimen", "transition"])
        table["S1:S2:S3 units"] = table.specimen.map(
            lambda n: ":".join(str(counts.loc[n].get(s, 0)) for s in ("PT-S1", "PT-S2", "PT-S3")))
        tables.append(table)
    table = pd.concat(tables, ignore_index=True)
    table["shift"] = table["equal-weighted"] - table["unweighted (72)"]
    pd.set_option("display.width", 200)
    print(table.round(3).to_string(index=False))
    print("\nGroup means")
    print(table.groupby(["dataset", "group", "transition"])[["unweighted (72)", "equal-weighted"]].mean().round(3))
    print("\nPooled over each dataset's four specimens (the cuts notebooks 79-84 use)")
    print(table.groupby(["dataset", "transition"])[["unweighted (72)", "equal-weighted"]].mean().round(3))

    # control: a mouse with balanced labels, S2 thinned to the human ratio
    position, segment, specimen = inputs(results, "hvm")
    rng = np.random.default_rng(0)
    print("\nMouse S1->S2 with S2 thinned")
    for name in ("Ctrl1A2", "Ctrl1A4"):
        s1 = np.flatnonzero((specimen == name) & (segment == "PT-S1"))
        s2 = np.flatnonzero((specimen == name) & (segment == "PT-S2"))
        for ratio in THIN_RATIOS:
            rows = np.r_[s1, rng.choice(s2, len(s1) // ratio, replace=False)]
            pair = [logistic_transitions(position[rows], segment[rows], specimen[rows], ("PT-S1", "PT-S2"),
                                         balanced=b).position.iloc[0] for b in (False, True)]
            print(f"  {name} S1:S2 {ratio}:1  unweighted {pair[0]:.3f}  equal-weighted {pair[1]:.3f}")


if __name__ == "__main__":
    main()
