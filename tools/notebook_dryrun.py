"""Execute every code cell of a workflow notebook against a synthetic dataset in a throwaway directory.

This is a test harness, not a notebook run: it builds a synthetic AnnData satisfying the notebook's data
contract, points PSEUDOSPACE_DATA_ROOT / PSEUDOSPACE_RESULTS_ROOT at a temporary tree, and executes
the notebook's own cell sources in order under the Agg backend. Nothing is written into the
repository and no private data is touched. Real-data scale, values and figure quality are not tested.

Notebook 06 is the default target; ``--notebook`` runs any other one. The synthetic tree carries the
artifacts several notebooks read (the PT object, 05's gene-rewiring tables, and the ortholog map and
pathway libraries 07 resolves membership through), so a notebook that needs none of them simply
ignores the extras.

Synthetic genes are built in regimes, so known groups land in known phenotype classes and the
pathway enrichment comes out non-empty - that is what exercises the pathway figures, which empty
synthetic results kept skipping. One regime is a displacement regime: there is no shifted phenotype
class any more, so it exercises the supporting-displacement path instead - the flag is set and the gene
is classed by what its residual shape does.
"""
import json, os, pathlib, shutil, sys, tempfile, traceback
from contextlib import contextmanager

import numpy as np
import pandas as pd

REPO = pathlib.Path(__file__).resolve().parent.parent
DEFAULT_NOTEBOOK = "analysis/notebooks/06_human_mouse_spatial_rewiring.ipynb"


def _notebook_path() -> pathlib.Path:
    """Notebook to run, from --notebook or the default, resolved against the repo root."""
    args = sys.argv[1:]
    if "--notebook" in args:
        return (REPO / args[args.index("--notebook") + 1]).resolve()
    return REPO / DEFAULT_NOTEBOOK

os.environ["MPLCONFIGDIR"] = "/tmp/mpl_dry"
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad

@contextmanager
def _noop(directory):
    """Stand-in for TemporaryDirectory when the harness must keep its outputs."""
    directory.mkdir(parents=True, exist_ok=True)
    yield directory


LANDMARKS = ["Slc5a2", "Slc5a12", "Gatm", "Lrp2", "Cubn", "Slc34a1",
             "Slc22a6", "Slc13a3", "Cyp2e1", "Slc22a7", "Cyp7b1", "Slc7a13", "Slc6a18", "Acsm3"]
LITERATURE = ["Vcam1", "Spp1", "Prom1", "Vim", "S100a6", "Anxa4"]
SAMPLES = {"Ctrl1A2": "mouse", "Ctrl1A4": "mouse", "HUK1_COR1": "human", "HUK1_MED1": "human"}
N_PER_SAMPLE = 150


def build_tree(root: pathlib.Path):
    """Write a synthetic PT object, pathway membership and 05's artifacts under `root`."""
    results = root / "results"
    upstream = results / "human_vs_healthy_mouse"
    for sub in ("curves", "gene_rewiring", "stage_cache"):
        (upstream / sub).mkdir(parents=True, exist_ok=True)
    (root / "data").mkdir(parents=True, exist_ok=True)

    genes = LANDMARKS + LITERATURE + [f"Gene{index}" for index in range(200)]
    n_genes = len(genes)
    head = len(LANDMARKS) + len(LITERATURE)

    rng = np.random.default_rng(0)
    regimes = np.array(["mixed"] * n_genes, dtype=object)
    regimes[head:head + 60] = "conserved"
    regimes[head + 60:head + 120] = "complex"
    regimes[head + 120:head + 180] = "shifted"
    regimes[head + 180:] = "amplitude"

    amplitude = rng.uniform(0.2, 0.6, size=n_genes)
    amplitude[:len(LANDMARKS)] = 0.5                  # patterned, so the figure-1A QC panel has content
    centre = rng.uniform(0.15, 0.85, size=n_genes)
    displacement = np.zeros(n_genes)
    scale = np.ones(n_genes)
    n_complex = int((regimes == "complex").sum())
    n_shifted = int((regimes == "shifted").sum())
    n_amplitude = int((regimes == "amplitude").sum())
    # Beyond the +-0.20 search cap, so the bounded phase search cannot rescue them -> complex.
    displacement[regimes == "complex"] = rng.choice([0.30, -0.30], size=n_complex)
    displacement[regimes == "shifted"] = rng.choice([0.12, -0.12], size=n_shifted)
    scale[regimes == "amplitude"] = rng.choice([0.3, 2.0], size=n_amplitude)

    # --- structures ---------------------------------------------------------------------
    positions, species, sample = [], [], []
    for name, kind in SAMPLES.items():
        for value in np.sort(rng.uniform(0.02, 1.0, size=N_PER_SAMPLE)):
            positions.append(value); species.append(kind); sample.append(name)
    position = np.asarray(positions)
    species = np.asarray(species)
    sample = np.asarray(sample)
    human = species == "human"

    values = np.zeros((position.size, n_genes))
    for index in range(n_genes):
        centre_here = np.where(human, centre[index] + displacement[index], centre[index])
        amplitude_here = np.where(human, amplitude[index] * scale[index], amplitude[index])
        bump = np.exp(-((position - centre_here) ** 2) / (2 * 0.08 ** 2))
        values[:, index] = amplitude_here * bump + rng.normal(0, 0.05, size=position.size)
    values = np.clip(values, 0.0, None)
    values[rng.random(values.shape) < 0.15] = 0.0      # zero inflation -> realistic detection

    obs = pd.DataFrame({
        "comparison_species": species, "sample": sample, "total_scanpy_dpt": position,
        "region": np.where(human, "cortex", "kidney"),
    }, index=[f"structure_{index}" for index in range(position.size)])
    var = pd.DataFrame({"measured_in_both_inputs": np.ones(n_genes, dtype=bool)}, index=genes)
    adata = ad.AnnData(X=values, obs=obs, var=var)
    adata.layers["lognorm"] = values.copy()
    adata.layers["counts"] = np.rint(values * 10).astype(np.int64)
    adata.write_h5ad(upstream / "cross_species_pt_dpt.h5ad")

    # --- pathway membership, built from the regimes so enrichment is non-empty -----------
    def group(name):
        return [genes[index] for index in np.flatnonzero(regimes == name)]

    random_choices = np.random.default_rng(7).choice(genes, size=(2, 40), replace=False)
    pathway_rows = [
        ("Reactome_2022", "Peroxisome R-HSA-900", group("conserved")[:30]),
        ("Reactome_2022", "Fatty acid metabolism R-HSA-901", group("complex")[:30]),
        ("Reactome_2022", "Oxidative phosphorylation R-HSA-902", group("shifted")[:30]),
        ("MSigDB_Hallmark_2020", "Hypoxia", group("amplitude")[:30]),
        ("MSigDB_Hallmark_2020", "TNF-alpha Signaling via NF-kB", list(random_choices[0])),
        ("MSigDB_Hallmark_2020", "Bile acid metabolism", list(random_choices[1])),
        ("KEGG_2019_Mouse", "Tiny set", ["Gene1"]),
    ]
    pd.DataFrame([
        {"library": library, "pathway": pathway, "n_requested": len(members),
         "n_with_ortholog": len(members), "n_assayed": len(members), "n_tested": len(members),
         "n_genes_present": len(members), "genes_present": repr(list(members)),
         "genes_assayed": repr(list(members)), "retained": True, "exclusion_reason": ""}
        for library, pathway, members in pathway_rows
    ]).to_csv(upstream / "curves" / "pathway_membership_coverage.csv", index=False)

    # --- 07's inputs: the accepted ortholog map and the pathway JSON libraries --------------------
    # 07 resolves pathway membership through the ortholog map itself instead of reading a coverage
    # table, so both artifacts have to exist. The sets are built from the regimes above, so 07's
    # enrichment lands on known groups instead of on nothing.
    pd.DataFrame({
        "human_symbol": genes, "mouse_symbol": genes, "mapping_status": "synthetic_identity",
    }).to_csv(upstream / "ortholog_map_used.csv", index=False)

    gene_sets_dir = root / "data" / "mouse_vs_human" / "pathway_gene_sets"
    gene_sets_dir.mkdir(parents=True, exist_ok=True)
    early = [gene for gene, centre_value in zip(genes, centre) if centre_value < 0.35]
    libraries = {
        "Reactome_2022.json": {
            "Early PT transport R-HSA-1": early[:30],
            "Fatty acid metabolism R-HSA-901": group("complex")[:30],
            "Oxidative phosphorylation R-HSA-902": group("shifted")[:30],
        },
        "MSigDB_Hallmark_2020.json": {
            "Hypoxia": group("amplitude")[:30],
            "TNF-alpha Signaling via NF-kB": group("conserved")[:30],
        },
        # 'Tiny set' is below the minimum member count, so it must be reported as excluded rather
        # than scored: that is the path an under-sized library entry takes.
        "KEGG_2019_Mouse.json": {
            "Bile acid metabolism": group("mixed")[:30],
            "Tiny set": ["Gene1"],
        },
    }
    for filename, sets in libraries.items():
        (gene_sets_dir / filename).write_text(json.dumps(sets))

    # --- 05's artifacts, read by the cross-check cells only -----------------------------
    grid = np.linspace(0.02, 1.0, 101)
    prior_mouse = rng.normal(size=(n_genes, grid.size)) * 0.3
    np.savez_compressed(
        upstream / "gene_rewiring" / "balanced_gene_curves.npz",
        genes=np.asarray(genes), grid=grid, grid_unit=(grid - grid[0]) / (grid[-1] - grid[0]),
        mouse=prior_mouse, human=prior_mouse + rng.normal(0, 0.05, size=(n_genes, grid.size)),
    )
    pd.DataFrame({
        "gene": genes, "mouse_amplitude": rng.uniform(0.05, 1.0, size=n_genes),
        "human_amplitude": rng.uniform(0.05, 1.0, size=n_genes),
        "level_effect_human_minus_mouse": rng.normal(size=n_genes),
    }).to_csv(upstream / "gene_rewiring" / "gene_curve_atlas.csv", index=False)
    pd.DataFrame({
        "gene": genes, "mean_log2_cpm_mouse": rng.normal(size=n_genes),
        "mean_log2_cpm_human": rng.normal(size=n_genes),
        "log2fc_human_vs_mouse": rng.normal(size=n_genes),
    }).to_csv(upstream / "gene_rewiring" / "whole_PT_gene_logFC.csv", index=False)
    return upstream


def main():
    # --keep writes into a persistent directory so the produced tables can be inspected.
    keep = "--keep" in sys.argv
    quiet = "--verbose" not in sys.argv
    notebook_path = _notebook_path()
    persistent = pathlib.Path("/tmp/dryrun_keep")
    if keep and persistent.exists():
        shutil.rmtree(persistent)
    with (tempfile.TemporaryDirectory() if not keep else _noop(persistent)) as temporary:
        root = pathlib.Path(persistent if keep else temporary)
        upstream = build_tree(root)
        os.environ["PSEUDOSPACE_DATA_ROOT"] = str(root / "data")
        os.environ["PSEUDOSPACE_RESULTS_ROOT"] = str(root / "results")

        notebook = json.load(open(notebook_path))
        if not quiet:
            print(f"# dry-run of {notebook_path.relative_to(REPO)}"
                  f"  ({len(notebook['cells'])} cells)")
        sources = [("".join(cell["source"]), index) for index, cell in enumerate(notebook["cells"])
                   if cell["cell_type"] == "code"]

        namespace = {"__name__": "__main__"}
        failures = []
        for source, index in sources:
            try:
                exec(compile(source, f"cell_{index}", "exec"), namespace)
            except Exception:
                failures.append((index, source.split("\n")[0][:70], traceback.format_exc()))
                print(f"!!! FAILED cell {index}: {source.split(chr(10))[0][:70]}")
            else:
                if not quiet:
                    print(f"    ok  cell {index}")

        print("\n" + "=" * 90)
        if failures:
            print(f"{len(failures)} cell(s) failed:")
            for index, first_line, trace in failures:
                print(f"\n--- cell {index}: {first_line}")
                for line in trace.rstrip().splitlines()[-12:]:
                    print(f"    {line}")
        else:
            print("every code cell executed without error")
        figures = sorted((upstream / "spatial_rewiring/figures").glob("*.png"))
        tables = sorted((upstream / "spatial_rewiring/tables").glob("*.csv"))
        # Every notebook's own outputs, so a --notebook run other than 06 also reports what it wrote.
        # `stage_cache` payloads are compute artifacts, not results.
        produced = sorted(
            path.relative_to(root) for path in (root / "results").rglob("*")
            if path.is_file() and path.suffix in (".csv", ".png") and "stage_cache" not in path.parts
        )
        print(f"\ndry-run: {len(sources) - len(failures)}/{len(sources)} cells ok, "
              f"{len(figures)} figures, {len(tables)} tables")
        if not quiet:
            for path in figures + tables:
                print(f"  {path.name}")
            for path in produced:
                print(f"  {path}")
        return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
