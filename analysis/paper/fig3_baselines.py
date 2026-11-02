"""Fig. 3e trajectory baselines, scored the same way as pseudostructure t.

Methods: stage-04 DPT (already computed), scFates `tl.tree` method `ppt`, Slingshot and GASTON. The first three
use the fit's pass-2 Harmony embedding (all 50 dimensions) and start at the stage-04 root unit. GASTON is spatial:
it runs per section on unit centroids plus GLM-PCs of the unit counts, and its isodepth is oriented so the units
in the root's coarse Leiden cluster (stage 04's root rule) have the lower median. Each coordinate is scored with
`fig3_common.section_rho`: unit-level Spearman rho against the nephron-order rank of the label, per section.
Protocol: docs/paper/fig3-validation-experiments.md, section 3e.

Reads the package runs under --results-root (outside Git); writes <results-root>/paper_figures/fig3/baselines/:
coordinates_<run>.parquet, baseline_rho.csv, run_meta.json and methods_notes.md. Finished method coordinates are
kept under baselines/work/ and reused unless --force.

External environments, passed by flag or environment variable (never installed into kidney-pseudospace):
  --rscript        an Rscript with slingshot and DelayedMatrixStats (conda-forge + bioconda: r-base 4.4,
                   bioconductor-slingshot, bioconductor-delayedmatrixstats); runs fig3_slingshot.R
  --gaston-python  a Python with gaston-spatial, glmpca and CPU torch; it re-invokes this file with --gaston-worker
Threads: FIG3_THREADS (default 16) caps BLAS/OpenMP/numba threads and GASTON restart workers.

  python analysis/paper/fig3_baselines.py --rscript <r-env>/bin/Rscript --gaston-python <gaston-env>/bin/python
"""
from __future__ import annotations

import os

THREADS = int(os.environ.get('FIG3_THREADS', '16'))
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ.setdefault(_v, str(THREADS))

import argparse  # noqa: E402
import contextlib  # noqa: E402
import io  # noqa: E402
import json  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

SEED = 0  # the project's scFates seed (notebook 13)
# scFates ppt: the scFates "Explore sigma" tutorial's PC-space recipe; 50 nodes, as the package's elastic curve
SCFATES = {'Nodes': 50, 'sigmas': [1000, 100, 10, 1, 0.1, 0.01], 'ppt_lambda': 100, 'ppt_nsteps': 50,
           'n_map': 1, 'seed': SEED}
SLINGSHOT = {'approx_points': 150}
# GASTON: its own Visium HD tutorial (GLM-PCA settings, NN defaults, 30 restarts); the lowest final loss wins
GASTON = {'num_dims': 8, 'num_genes': 5000, 'penalty': 100, 'glmpca_max_iter': 30, 'glmpca_eps': 1e-4,
          'isodepth_arch': [20, 20], 'expression_arch': [20, 20], 'epochs': 10000, 'optimizer': 'adam',
          'restarts': 30}
METHODS = ('scfates', 'slingshot', 'gaston')


def _repo() -> Path:
    return next(p for p in [Path(__file__).resolve(), *Path(__file__).resolve().parents] if (p / 'pseudospace').is_dir())


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--data-root', default=os.environ.get('PSEUDOSPACE_DATA_ROOT'))
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT'))
    p.add_argument('--runs', nargs='+', default=['aki', 'hvm'])
    p.add_argument('--methods', nargs='+', default=list(METHODS), choices=METHODS)
    p.add_argument('--rscript', default=os.environ.get('FIG3_RSCRIPT', 'Rscript'))
    p.add_argument('--gaston-python', default=os.environ.get('FIG3_GASTON_PYTHON'))
    p.add_argument('--force', action='store_true', help='refit methods whose coordinates are already in work/')
    p.add_argument('--gaston-worker', help=argparse.SUPPRESS)  # internal: GLM-PCA + NN restarts on one section folder
    return p.parse_args(argv)


# ---------------------------------------------------------------- GASTON worker (runs in the GASTON environment)

def _gaston_restart(task):
    folder, seed = task
    import torch
    from gaston import neural_net
    torch.set_num_threads(1)
    d = np.load(Path(folder) / 'nn_input.npz')
    S, A = torch.tensor(d['S']), torch.tensor(d['A'])
    t0 = time.time()
    mod, _ = neural_net.train(S, A, S_hidden_list=GASTON['isodepth_arch'], A_hidden_list=GASTON['expression_arch'],
                              epochs=GASTON['epochs'], checkpoint=GASTON['epochs'], optim=GASTON['optimizer'],
                              seed=seed, device='cpu')
    with torch.no_grad():
        loss = float(neural_net.get_loss(mod, S, A))  # process_NN_output.process_files' selection criterion
        iso = mod.spatial_embedding(S).numpy().ravel()
    return seed, loss, iso, time.time() - t0


def gaston_worker(folder: Path):
    """GLM-PCA of one section's counts, then GASTON restarts; keeps the isodepth of the lowest-loss model."""
    import multiprocessing as mp
    from glmpca import glmpca
    from gaston import neural_net
    d = np.load(folder / 'input.npz')
    counts, coords = d['counts'].astype(float), d['coords']
    t0 = time.time()
    res = glmpca.glmpca(counts.T, GASTON['num_dims'], fam='poi', penalty=GASTON['penalty'], verbose=False,
                        ctl={'maxIter': GASTON['glmpca_max_iter'], 'eps': GASTON['glmpca_eps'], 'optimizeTheta': True})
    A = np.asarray(res['factors'], dtype=float)
    glmpca_s = time.time() - t0
    if not np.isfinite(A).all():
        raise ValueError(f'{folder.name}: GLM-PCA gave non-finite factors')
    S_t, A_t = neural_net.load_rescale_input_data(coords, A)  # z-scores both, as the tutorial does
    np.savez(folder / 'nn_input.npz', S=S_t.numpy(), A=A_t.numpy())
    t1 = time.time()
    with mp.get_context('spawn').Pool(min(THREADS, GASTON['restarts'])) as pool:
        out = pool.map(_gaston_restart, [(str(folder), s) for s in range(GASTON['restarts'])])
    seeds, losses, isos, secs = zip(*out)
    best = int(np.argmin(losses))
    np.savez(folder / 'output.npz', isodepth=isos[best], glmpcs=A, losses=np.array(losses), seeds=np.array(seeds),
             restart_seconds=np.array(secs), best_seed=seeds[best], glmpca_seconds=glmpca_s,
             nn_seconds=time.time() - t1, glmpca_iters=len(res.get('dev', [])))


# ---------------------------------------------------------------- fits (kidney-pseudospace environment)

def fit_scfates(X, root_idx):
    import anndata as ad
    import pandas as pd
    import scFates as scf
    obj = ad.AnnData(obs=pd.DataFrame(index=[str(i) for i in range(len(X))]), obsm={'X_harmony': X})
    c = SCFATES
    sigma = scf.tl.explore_sigma(obj, Nodes=c['Nodes'], use_rep='X_harmony', sigmas=c['sigmas'], seed=c['seed'])
    log = io.StringIO()
    with contextlib.redirect_stdout(log):  # simpleppt prints whether the fit converged within ppt_nsteps
        scf.tl.tree(obj, Nodes=c['Nodes'], use_rep='X_harmony', method='ppt', ppt_sigma=sigma,
                    ppt_lambda=c['ppt_lambda'], ppt_nsteps=c['ppt_nsteps'], seed=c['seed'])
    print(log.getvalue(), end='')
    status = [ln.strip() for ln in log.getvalue().splitlines() if 'converged' in ln]
    g = obj.uns['graph']
    F = np.asarray(g['F'], dtype=float)  # dims x nodes
    node = int(np.argmin(np.linalg.norm(F - X[root_idx][:, None], axis=0)))
    scf.tl.root(obj, node)
    scf.tl.pseudotime(obj, n_jobs=2, n_map=c['n_map'], seed=c['seed'])
    t = obj.obs['t'].to_numpy(float)
    if not np.isfinite(t).all():
        raise ValueError('scFates gave non-finite pseudotime')
    tips = [int(v) for v in np.atleast_1d(g['tips'])]
    meta = {'sigma': float(sigma), 'ppt_status': status[-1] if status else 'unknown', 'n_nodes_kept': int(F.shape[1]), 'root_node': node, 'root_node_is_tip': node in tips,
            'n_tips': len(tips), 'n_forks': len(np.atleast_1d(g['forks'])), 'n_segments': int(obj.obs['seg'].nunique())}
    return t, meta


def fit_slingshot(X, units, clusters, start, folder: Path, rscript: str):
    import pandas as pd
    folder.mkdir(parents=True, exist_ok=True)
    np.ascontiguousarray(X, dtype='<f8').tofile(folder / 'embedding.f64')
    (folder / 'dims.txt').write_text(f'{X.shape[0]}\n{X.shape[1]}\n')
    (folder / 'units.txt').write_text('\n'.join(units) + '\n')
    (folder / 'clusters.txt').write_text('\n'.join(clusters) + '\n')
    (folder / 'start.txt').write_text(f'{start}\n')
    script = Path(__file__).with_name('fig3_slingshot.R')
    subprocess.run([rscript, '--vanilla', str(script), str(folder), str(SLINGSHOT['approx_points'])], check=True)
    out = pd.read_csv(folder / 'slingshot.csv').set_index('unit')
    if list(out.index.astype(str)) != list(units):
        raise ValueError('Slingshot output units do not match the input order')
    meta = {'start_cluster': start, 'lineages': (folder / 'lineages.txt').read_text().splitlines(),
            'n_lineages': int(sum(c.startswith('pt_') for c in out.columns)), 'r_session': (folder / 'session.txt').read_text().splitlines(),
            'pseudotime': 'slingAvgPseudotime = rowSums(slingPseudotime(na=FALSE) * slingCurveWeights(as.probs=TRUE))'}
    return out['slingshot'].to_numpy(float), meta


def fit_gaston_section(counts, coords, folder: Path, gaston_python: str):
    from scipy import sparse
    folder.mkdir(parents=True, exist_ok=True)
    tot = np.asarray(counts.sum(axis=0)).ravel()
    genes = np.sort(np.argsort(tot, kind='stable')[-GASTON['num_genes']:])  # top genes by section UMI, as the tutorial
    Y = counts[:, genes]
    Y = Y.toarray() if sparse.issparse(Y) else np.asarray(Y)
    if not np.array_equal(Y, np.round(Y)) or Y.min() < 0:
        raise ValueError(f'{folder.name}: layers["counts"] must hold non-negative integer counts')
    np.savez(folder / 'input.npz', counts=Y.astype(np.int32), coords=np.asarray(coords, dtype=float))
    # -I: the GASTON environment never imports from the script directory or the working directory
    subprocess.run([gaston_python, '-I', str(Path(__file__).resolve()), '--gaston-worker', str(folder)], check=True,
                   cwd=folder)
    (folder / 'input.npz').unlink()  # the dense count block is large and is rebuilt from pass2.h5ad on a rerun
    o = np.load(folder / 'output.npz')
    meta = {'n_units': int(Y.shape[0]), 'n_genes': int(Y.shape[1]), 'best_seed': int(o['best_seed']),
            'best_loss': float(o['losses'].min()), 'loss_range': [float(o['losses'].min()), float(o['losses'].max())],
            'glmpca_seconds': round(float(o['glmpca_seconds']), 1), 'nn_seconds': round(float(o['nn_seconds']), 1),
            'restart_seconds_median': round(float(np.median(o['restart_seconds'])), 1)}
    return o['isodepth'].astype(float), meta


def _versions(gaston_python: str | None):
    from importlib.metadata import version
    v = {k: version(k) for k in ('scFates', 'simpleppt', 'anndata', 'numpy', 'scipy', 'pandas')}
    v['python'] = sys.version.split()[0]
    if gaston_python:
        code = ("from importlib.metadata import version as v; import sys, json; print(json.dumps({k: v(k) for k in "
                "('gaston-spatial','glmpca','torch','numpy','scikit-learn','statsmodels')} | {'python': sys.version.split()[0]}))")
        v['gaston_env'] = json.loads(subprocess.run([gaston_python, '-I', '-c', code], check=True, capture_output=True,
                                                    text=True).stdout)
    return v


def run_one(run, results_root: Path, out: Path, args, meta):
    import anndata as ad
    import pandas as pd
    from fig3_common import RUNS, SECTIONS, load_positions, section_rho

    base = results_root / RUNS[run]
    table = load_positions(results_root, run)
    a = ad.read_h5ad(base / '03_integration/pass2.h5ad', backed='r')
    root = json.loads((base / '04_pseudospace/manifest.json').read_text())['summary']['root_unit']
    if not table.index.equals(pd.Index(a.obs_names)):
        raise ValueError(f'{run}: pseudospace.parquet units differ from pass2.h5ad')
    if root not in table.index:
        raise ValueError(f'{run}: root unit {root} not on the curve')
    obs = a.obs
    X = np.asarray(a.obsm['X_harmony'], dtype=float)
    if X.shape[1] != 50 or not np.isfinite(X).all():
        raise ValueError(f'{run}: X_harmony must be finite with 50 dimensions, got {X.shape}')
    root_idx = int(table.index.get_loc(root))
    units = [str(u) for u in table.index]
    coords = table[['sample', 't']].copy()
    if 'dpt' in table:
        coords['dpt'] = table['dpt']
    m = meta.setdefault(run, {'root_unit': root, 'root_sample': str(obs.loc[root, 'sample']),
                              'root_leiden_pass2': str(obs.loc[root, 'leiden_pass2']),
                              'root_leiden_coarse': str(obs.loc[root, 'leiden_coarse']), 'n_units': len(units)})
    work = out / 'work'
    work.mkdir(parents=True, exist_ok=True)

    for method in METHODS:
        cache = work / f'{run}_{method}.parquet'
        if cache.exists() and not (args.force and method in args.methods):
            prev = pd.read_parquet(cache)
            coords[method] = prev[method].reindex(coords.index)
            info_path = cache.with_suffix('.json')  # fit details travel with the cached coordinate
            if info_path.exists():
                m[method] = json.loads(info_path.read_text())
            continue
        if method not in args.methods:
            coords[method] = np.nan
            continue
        t0 = time.time()
        if method == 'scfates':
            values, info = fit_scfates(X, root_idx)
        elif method == 'slingshot':
            values, info = fit_slingshot(X, units, obs['leiden_pass2'].astype(str).tolist(),
                                         str(obs.loc[root, 'leiden_pass2']), work / f'slingshot_{run}', args.rscript)
        else:
            if not args.gaston_python:
                raise SystemExit('--gaston-python (or FIG3_GASTON_PYTHON) is required for GASTON')
            counts = a.layers['counts']
            values = np.full(len(units), np.nan)
            info = {'orientation_rule': "per section, units with leiden_coarse == the root unit's leiden_coarse "
                                        "(stage 04's root rule) get the lower median isodepth"}
            root_coarse = str(obs.loc[root, 'leiden_coarse'])
            for s in sorted(table['sample'].unique()):
                k = np.flatnonzero((table['sample'] == s).to_numpy())
                iso, sec = fit_gaston_section(counts[k], obs.iloc[k][['centroid_x', 'centroid_y']].to_numpy(float),
                                              work / f'gaston_{run}_{s}', args.gaston_python)
                in_root = obs.iloc[k]['leiden_coarse'].astype(str).eq(root_coarse).to_numpy()
                if not 0 < in_root.sum() < len(k):
                    raise ValueError(f'{run} {s}: the root coarse cluster cannot orient this section')
                flip = bool(np.median(iso[in_root]) > np.median(iso[~in_root]))
                values[k] = -iso if flip else iso
                info[s] = sec | {'n_root_cluster_units': int(in_root.sum()), 'flipped': flip}
                print(f'[gaston] {run} {s}: best seed {sec["best_seed"]}, flipped {flip}', flush=True)
        info['seconds'] = round(time.time() - t0, 1)
        info['threads'] = THREADS  # BLAS/OpenMP/numba threads; GASTON restart workers
        m[method] = info
        coords[method] = values
        coords[[method]].to_parquet(cache)
        cache.with_suffix('.json').write_text(json.dumps(info, indent=1, default=str))
        print(f'[{method}] {run} done in {info["seconds"]} s', flush=True)

    coords = coords[['t', 'dpt', 'scfates', 'slingshot', 'gaston', 'sample']]
    coords.index.name = 'obs_name'
    coords.to_parquet(out / f'coordinates_{run}.parquet')
    sections = [s for s, r in SECTIONS.items() if r == run]
    rows = []
    for method in ('t', 'dpt', 'scfates', 'slingshot', 'gaston'):
        if coords[method].notna().any():
            rows.append(section_rho(coords[method], table, sections).assign(run=run, method=method))
    return pd.concat(rows, ignore_index=True)[['run', 'method', 'sample', 'n', 'rho']]


def write_notes(out: Path, meta, rho):
    v = meta['versions']
    lines = ['Fig. 3e trajectory baselines: methods notes (generated by analysis/paper/fig3_baselines.py)', '',
             'Protocol: docs/paper/fig3-validation-experiments.md, section 3e.',
             'Score: unit-level Spearman rho between each coordinate and the nephron-order rank of the unit label, per '
             'section (fig3_common.section_rho). Ctrl1A2, Ctrl1A4, IR2A2, IR2A4 are scored in the aki fit; HUK1_COR1 and '
             'HUK1_MED1 in the hvm fit.', '',
             'Shared inputs: <run>/03_integration/pass2.h5ad (read backed), obsm X_harmony, all 50 dimensions; root = '
             '<run>/04_pseudospace/manifest.json summary.root_unit.', '',
             f"Versions (kidney-pseudospace): {json.dumps({k: x for k, x in v.items() if k != 'gaston_env'})}",
             f"Versions (GASTON env): {json.dumps(v.get('gaston_env'))}",
             'Environments: scFates runs in kidney-pseudospace; Slingshot in a separate conda env (conda-forge + bioconda, '
             'r-base 4.4, bioconductor-slingshot, bioconductor-delayedmatrixstats); GASTON in a separate venv '
             '(gaston-spatial, glmpca, statsmodels, CPU torch). Slingshot R session: see the slingshot r_session entries.', '',
             't: pseudostructure (stage-04 elastic curve), column pseudospace of 04_pseudospace/pseudospace.parquet.',
             'dpt: stage-04 DPT baseline, column dpt_pseudotime of the same file; scored as is.', '',
             f"scfates: scFates tl.explore_sigma(Nodes={SCFATES['Nodes']}, use_rep='X_harmony', sigmas={SCFATES['sigmas']}, "
             f"seed={SEED}) picks ppt_sigma without labels; then tl.tree(method='ppt', Nodes={SCFATES['Nodes']}, "
             f"use_rep='X_harmony', ppt_sigma=<explored>, ppt_lambda={SCFATES['ppt_lambda']}, ppt_nsteps={SCFATES['ppt_nsteps']}, "
             f"seed={SEED}, metric euclidean, device cpu). This is the scFates Explore-sigma tutorial's PC-space recipe; "
             'the project had no earlier ppt settings (notebook 13 used tl.curve), its seed 0 is kept. tl.root at the '
             f"principal point nearest the root unit in X_harmony (Euclidean); tl.pseudotime(n_jobs=2, n_map={SCFATES['n_map']}, "
             f'seed={SEED}), column t. A branched tree is used as given.',
             f"slingshot: slingshot(X_harmony, clusterLabels=obs leiden_pass2, start.clus=<root unit's leiden_pass2>, "
             f"approx_points={SLINGSHOT['approx_points']}, all other arguments default) in R (analysis/paper/fig3_slingshot.R); "
             'coordinate = slingAvgPseudotime, i.e. rowSums(slingPseudotime(na=FALSE) * slingCurveWeights(as.probs=TRUE)).',
             f"gaston: per section; GLM-PCA (glmpca.glmpca, fam='poi', {GASTON['num_dims']} factors, penalty {GASTON['penalty']}, "
             f"maxIter {GASTON['glmpca_max_iter']}, eps {GASTON['glmpca_eps']}, optimizeTheta True) on layers['counts'] of the "
             f"top {GASTON['num_genes']} genes by section UMI total (GASTON Visium HD tutorial); spatial input = obs centroid_x, "
             f"centroid_y; both z-scored by gaston.neural_net.load_rescale_input_data; neural_net.train(S_hidden_list="
             f"{GASTON['isodepth_arch']}, A_hidden_list={GASTON['expression_arch']}, epochs={GASTON['epochs']}, optim="
             f"'{GASTON['optimizer']}', lr 1e-3, full batch, ReLU, device cpu) for seeds 0..{GASTON['restarts'] - 1}; the "
             'model with the lowest final MSE (process_NN_output.process_files criterion) gives the isodepth, '
             "spatial_embedding(S). Orientation: within each section, units with leiden_coarse equal to the root unit's "
             "leiden_coarse (stage 04's root rule, coarse cluster 0) get the lower median isodepth; no labels and no t "
             'are used.', '', 'Per-run fit details:']
    for run in (r for r in ('aki', 'hvm') if r in meta):
        lines.append(f'{run}: {json.dumps(meta[run], default=str)}')
    lines += ['', 'Wall-clock seconds per fit (whole run; GASTON is the sum over its sections):']
    for run in (r for r in ('aki', 'hvm') if r in meta):
        lines.append(f"{run}: " + ', '.join(f"{k} {meta[run][k].get('seconds')}" for k in METHODS if k in meta[run]))
    lines += ['', 'rho per section:', rho.pivot_table(index='method', columns='sample', values='rho').round(3).to_string()]
    (out / 'methods_notes.md').write_text('\n'.join(lines) + '\n')


def main(argv=None):
    args = parse_args(argv)
    if args.gaston_worker:
        gaston_worker(Path(args.gaston_worker))
        return
    import pandas as pd
    repo = _repo()
    sys.path.insert(0, str(repo / 'analysis/paper'))
    results_root = Path(args.results_root or repo / 'results')
    out = results_root / 'paper_figures/fig3/baselines'
    out.mkdir(parents=True, exist_ok=True)
    meta_path = out / 'run_meta.json'
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    meta['versions'] = _versions(args.gaston_python)
    rho = pd.concat([run_one(r, results_root, out, args, meta) for r in args.runs], ignore_index=True)
    prev = out / 'baseline_rho.csv'
    if prev.exists():  # keep rows of runs not refitted in this call
        old = pd.read_csv(prev)
        rho = pd.concat([old[~old['run'].isin(args.runs)], rho], ignore_index=True)
    rho.to_csv(prev, index=False)
    meta_path.write_text(json.dumps(meta, indent=1, default=str))
    write_notes(out, meta, rho)
    print(rho.pivot_table(index='method', columns='sample', values='rho').round(3).to_string())


if __name__ == '__main__':
    main()
