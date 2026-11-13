"""Fig. 3d robustness: how stable is the pseudostructure coordinate t when its inputs are perturbed?

Each rerun perturbs the stage-02 QC units of a reference package run, then reruns the frozen stage 03 (two-pass
Harmony, Leiden, automatic labels) and stage 04 (elastic principal curve) with the run's own configuration, the
R harmony backend and a fixed BLAS/OpenMP thread count, as `analysis/scripts/position_half_sampling.py` does.
Protocol: docs/paper/fig3-validation-experiments.md, section 3d.

Perturbations (fractions 10, 20, 30, 50, 70, 100 %; seeds 0, 1, 2):
  tubules    a random subset of the QC units, stratified by section (the same fraction of every section)
  reads      binomial thinning of every unit's raw counts; stage 02's lognorm (library size to 1e4, log1p over all
             genes) is recomputed from the thinned counts, and stage 03 renormalises after its gene filter as usual
  genes      a random subset of the genes offered to the integration HVG selection. The hook is the QC stage's
             var column `exclude_from_harmony_hvg`: genes outside the subset are flagged, so they stay in the
             object (the HVG mean bins are formed over the offered genes only, as with any exclusion) and every
             later step that reads all genes (gene filter, clustering feature space, panel scores, the stage-04
             marker axis and its control genes) still sees them. The 2,000 clustering-feature HVGs, which only
             feed marker ranking and panel scores, are NOT restricted: the package has no hook for them.
  gene_half  five complementary splits of the same reference pool into halves A and B; each half is a rerun
Scores (on the units each rerun shares with the reference): Spearman rho against the reference t, rho against
the reference labels' nephron-order rank per section and pooled, the rerun's own stage-04 `segment_order_rho`,
and whether stage 04 oriented it the other way (rho < 0; recorded, never corrected).

Reads private data and writes under --results-root (outside Git). Resumable: finished reruns are skipped.

  python analysis/paper/fig3_robustness.py --run aki --jobs 8          # every rerun not yet done
  python analysis/paper/fig3_robustness.py --run aki --only reads_f100_s0
  python analysis/paper/fig3_robustness.py --run aki --score           # robustness.csv, gene_halves.csv
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd


def _repo() -> Path:
    return next(p for p in [Path(__file__).resolve(), *Path(__file__).resolve().parents] if (p / 'pseudospace').is_dir())


REPO = _repo()
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from fig3_common import RUNS, load_positions, section_rho  # noqa: E402

OUT = 'paper_figures/fig3/robustness'
FRACTIONS = (0.1, 0.2, 0.3, 0.5, 0.7, 1.0)
SEEDS = (0, 1, 2)
N_SPLITS = 5
CONDITIONS = ('tubules', 'reads', 'genes')
ARM_ID = {'tubules': 1, 'reads': 2, 'genes': 3, 'gene_half': 4, 'exclude': 5}   # keeps the arms' random streams
#                                                                                     independent
EXCLUDE_COL = 'exclude_from_harmony_hvg'


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--data-root', default=os.environ.get('PSEUDOSPACE_DATA_ROOT', str(REPO / 'data')))
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--run', choices=sorted(RUNS), default='aki')
    p.add_argument('--jobs', type=int, default=8)
    p.add_argument('--threads', type=int, default=4, help='BLAS/OpenMP threads per rerun; FIXED by the protocol '
                   '(Leiden partitions change with the thread count)')
    p.add_argument('--tmp', default='/dev/shm' if Path('/dev/shm').is_dir() else '', help='scratch for stage files')
    p.add_argument('--only', nargs='+', help='run only these reruns (names as in reruns/)')
    p.add_argument('--job', help=argparse.SUPPRESS)        # internal: run one rerun in this process
    p.add_argument('--score', action='store_true', help='score the finished reruns and write the tables')
    return p.parse_args(argv)


# ---------------------------------------------------------------------------------------------------------------
# reruns
# ---------------------------------------------------------------------------------------------------------------
def job_names() -> list[str]:
    names = [f'{c}_f{round(f * 100):03d}_s{s}' for c in CONDITIONS for f in FRACTIONS for s in SEEDS]
    return names + [f'gene_half_s{k}_{h}' for k in range(N_SPLITS) for h in 'AB']


def parse_job(name: str) -> dict:
    if name.startswith('gene_half_'):
        _, _, seed, half = name.split('_')
        return {'condition': 'gene_half', 'fraction': 0.5, 'seed': int(seed[1:]), 'half': half}
    condition, fraction, seed = name.split('_')
    assert condition in CONDITIONS, name
    return {'condition': condition, 'fraction': int(fraction[1:]) / 100, 'seed': int(seed[1:]), 'half': ''}


def perturb(qc, spec: dict, target_sum: float) -> tuple[object, dict]:
    """The QC units of one rerun (a subset copy, or `qc` edited in place), and a record of what was changed."""
    import scanpy as sc
    from scipy import sparse

    cond, f = spec['condition'], spec['fraction']
    rng = np.random.default_rng([ARM_ID[cond], spec['seed']])
    info = {}
    if cond == 'tubules':
        section = qc.obs['sample'].astype(str).to_numpy()
        keep = np.zeros(qc.n_obs, bool)
        for s in sorted(set(section)):
            idx = np.flatnonzero(section == s)
            keep[rng.permutation(idx)[:max(1, round(f * len(idx)))]] = True
        info['n_kept_by_section'] = pd.Series(section[keep]).value_counts().sort_index().to_dict()
        return qc[keep].copy(), info
    if cond == 'reads':                              # reads and genes edit the caller's object in place
        counts = sparse.csr_matrix(qc.layers['counts'], dtype=np.float32, copy=True)
        before = float(counts.sum())
        counts.data = rng.binomial(counts.data.astype(np.int64), f).astype(np.float32)
        counts.eliminate_zeros()
        lib = np.asarray(counts.sum(1)).ravel()
        if (lib == 0).any():
            raise ValueError(f'{int((lib == 0).sum())} units have no counts left after thinning')
        tmp = sc.AnnData(X=counts.copy())            # stage 02 `normalize_log1p`: library size over all genes, log1p
        sc.pp.normalize_total(tmp, target_sum=target_sum)
        sc.pp.log1p(tmp)
        qc.layers['counts'], qc.layers['lognorm'] = counts, tmp.X
        qc.X = tmp.X.copy()
        info.update(read_fraction_realised=float(counts.sum()) / before, median_library=float(np.median(lib)))
        return qc, info
    excluded = qc.var[EXCLUDE_COL].to_numpy(bool)
    if cond == 'exclude':    # a named gene list withdrawn from the pool (Fig. 6 injury-blind axis)
        withdrawn = qc.var_names.isin(spec['genes']) & ~excluded
        qc.var[EXCLUDE_COL] = excluded | withdrawn
        info.update(n_pool=int((~excluded).sum()), n_withdrawn=int(withdrawn.sum()))
        return qc, info
    pool = qc.var_names[~excluded]
    perm = rng.permutation(len(pool))
    if cond == 'genes':
        offered = pool[perm[:round(f * len(pool))]]
    else:   # gene_half: complementary halves of one permutation
        cut = len(pool) // 2
        offered = pool[perm[:cut] if spec['half'] == 'A' else perm[cut:]]
    qc.var[EXCLUDE_COL] = excluded | ~qc.var_names.isin(offered)
    info.update(n_pool=int(len(pool)), n_offered=int(len(offered)))
    return qc, info


def run_job(args, name: str, out: Path, spec: dict | None = None) -> None:
    """One rerun: perturb the QC units, then COAST's two passes and curve (coast_fit.two_pass), as the reference fit."""
    import anndata as ad
    import coast_fit

    t0 = time.time()
    spec = spec or parse_job(name)
    P = coast_fit.panels(args.run)
    qc = ad.read_h5ad(Path(args.results_root) / RUNS[args.run] / '02_qc/units_qc.h5ad')
    n_qc = qc.n_obs
    qc, info = perturb(qc, spec, 1e4)
    try:
        withdrawn = qc.var_names[qc.var[EXCLUDE_COL].to_numpy(bool)]
        n_input = qc.n_obs
        p1, p2, _ = coast_fit.two_pass(qc, P, exclude_genes=withdrawn)
        del qc
        stage04 = coast_fit.stage04_table(p2)
        table = pd.DataFrame({'t': stage04['pseudospace'], 'sample': stage04['sample'],
                              'segment_class': stage04['segment_class'], 'coarse_class': stage04['coarse_class']},
                             index=stage04.index)
        om = coast_fit.order_metrics(stage04, P['segment_order']).iloc[0]
        meta = {'job': name, **spec, **info, 'method': 'COAST', 'n_qc_units': n_qc, 'n_input_units': n_input,
                'n_units_pass2': int(p2.n_obs), 'n_units': int(len(table)), 'n_genes_pass2': int(p2.n_vars),
                'n_integration_features_pass1': int(p1.var['coast_feature'].sum()),
                'n_integration_features_pass2': int(p2.var['coast_feature'].sum()),
                'segment_order_rho': float(om.segment_order_rho), 'n_inversions': int(om.n_inversions),
                'root_unit': str(stage04['pseudospace'].idxmin()), 'threads': args.threads,
                'seconds': round(time.time() - t0)}
        (out / f'{name}.json').write_text(json.dumps(meta, indent=1, default=str))
        table.to_parquet(out / f'{name}.parquet.part')
        (out / f'{name}.parquet.part').rename(out / f'{name}.parquet')     # the parquet marks a finished rerun
        (out / f'{name}.error.json').unlink(missing_ok=True)
    except Exception:
        (out / f'{name}.error.json').write_text(json.dumps({'job': name, **spec, 'seconds': round(time.time() - t0),
                                                            'error': traceback.format_exc()}, indent=1))
        raise


def launch(args, out: Path) -> None:
    todo = [n for n in (args.only or job_names()) if not (out / f'{n}.parquet').exists()]
    # longest first: the full-size reruns, then the tubule subsets from large to small
    todo.sort(key=lambda n: -parse_job(n)['fraction'] if n.startswith('tubules') else -2)
    env = {**os.environ, **{k: str(args.threads) for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')}}
    base = [sys.executable, str(Path(__file__).resolve()), '--data-root', args.data_root, '--results-root',
            args.results_root, '--run', args.run, '--threads', str(args.threads), '--tmp', args.tmp]
    running: list[tuple[str, subprocess.Popen]] = []
    failed = []
    t0 = time.time()
    print(f'{len(todo)} reruns to do, {args.jobs} at a time', flush=True)
    while todo or running:
        while todo and len(running) < args.jobs:
            n = todo.pop(0)
            log = open(out / f'{n}.log', 'w')
            running.append((n, subprocess.Popen([*base, '--job', n], env=env, stdout=log, stderr=subprocess.STDOUT)))
        time.sleep(5)
        for n, proc in list(running):
            if proc.poll() is not None:
                running.remove((n, proc))
                if proc.returncode:
                    failed.append(n)
                print(f'{n} {"FAILED" if proc.returncode else "done"} ({round(time.time() - t0)} s elapsed)', flush=True)
    print(f'finished in {round(time.time() - t0)} s; failed: {failed}', flush=True)


# ---------------------------------------------------------------------------------------------------------------
# scores
# ---------------------------------------------------------------------------------------------------------------
def _rho(a, b) -> float:
    from scipy.stats import spearmanr
    return float(spearmanr(a, b).statistic)


def label_scores(t: pd.Series, ref: pd.DataFrame, sections) -> dict:
    """rho of a coordinate against the REFERENCE labels' rank: per section and pooled over every section."""
    table = ref.loc[ref.index.intersection(t.index)]
    per = section_rho(t, table, sections).set_index('sample')['rho']
    k = table['rank'].notna() & np.isfinite(t.reindex(table.index))
    return {'rho_labels_pooled': _rho(t.reindex(table.index)[k], table.loc[k, 'rank']),
            **{f'rho_labels_{s}': float(per[s]) for s in sections}}


def score(args, out: Path, root: Path) -> None:
    ref = load_positions(args.results_root, args.run)
    sections = sorted(ref['sample'].unique())
    # the 100 % rerun at the protocol's thread count: rerunning the unperturbed data at 4 threads does not give the
    # reference's Leiden partition (the reference ran at the default thread count), so the score against this
    # rerun separates each perturbation's effect from the thread count's
    full_path = out / 'tubules_f100_s0.parquet'
    full = pd.read_parquet(full_path)['t'] if full_path.exists() else None
    rows = []
    for name in job_names():
        spec = parse_job(name)
        f = out / f'{name}.parquet'
        if not f.exists():
            err = out / f'{name}.error.json'
            if err.exists():
                rows.append({'run': args.run, 'job': name, **spec, 'error': json.loads(err.read_text())['error'].strip().splitlines()[-1]})
            continue
        t = pd.read_parquet(f)['t']
        meta = json.loads((out / f'{name}.json').read_text())
        shared = t.index.intersection(ref.index)
        rho_ref = _rho(t[shared], ref.loc[shared, 't'])
        med = t[shared].groupby(ref.loc[shared, 'group']).median()   # partial flips keep rho > 0: the PT limb
        #                                                              entered from its S3 end is one of them
        rows.append({'run': args.run, 'job': name, **spec, 'n_units': int(len(t)), 'rho_vs_reference': rho_ref,
                     **label_scores(t, ref, sections), 'segment_order_rho': meta['segment_order_rho'],
                     'seconds': meta['seconds'], 'n_shared': int(len(shared)), 'orientation_flipped': bool(rho_ref < 0),
                     'max_abs_dt': float((t[shared] - ref.loc[shared, 't']).abs().max()),
                     'pt_limb_reversed': bool(med.get('PTS1', np.nan) > med.get('PTS3', np.nan)),
                     'rho_vs_full_rerun': (_rho(t[both], full[both]) if full is not None and len(both := t.index.intersection(full.index))
                                           else np.nan),
                     'n_inversions': meta['n_inversions'], 'n_input_units': meta['n_input_units'],
                     'n_integration_features_pass1': meta['n_integration_features_pass1'],
                     'n_clusters_pass1': meta['fingerprint_pass1']['n_clusters'],
                     'fingerprint_pass1': meta['fingerprint_pass1']['membership_sha1'],
                     'fingerprint_pass2': meta['fingerprint_pass2']['membership_sha1'], 'root_unit': meta['root_unit'],
                     'error': ''})
    table = pd.DataFrame(rows)
    table.to_csv(root / 'robustness.csv', index=False)

    halves = []
    for k in range(N_SPLITS):
        a, b = (out / f'gene_half_s{k}_{h}.parquet' for h in 'AB')
        if not (a.exists() and b.exists()):
            continue
        ta, tb = pd.read_parquet(a)['t'], pd.read_parquet(b)['t']
        both = ta.index.intersection(tb.index)
        r = table.set_index(['condition', 'seed', 'half']).loc['gene_half']
        halves.append({'run': args.run, 'split_seed': k, 'n_A': len(ta), 'n_B': len(tb), 'n_shared_AB': len(both),
                       'rho_A_vs_B': _rho(ta[both], tb[both]), 'rho_A_vs_reference': r.loc[(k, 'A'), 'rho_vs_reference'],
                       'rho_B_vs_reference': r.loc[(k, 'B'), 'rho_vs_reference'],
                       'rho_labels_pooled_A': r.loc[(k, 'A'), 'rho_labels_pooled'],
                       'rho_labels_pooled_B': r.loc[(k, 'B'), 'rho_labels_pooled']})
    pd.DataFrame(halves).to_csv(root / 'gene_halves.csv', index=False)
    pd.DataFrame([{'run': args.run, 'n_units': len(ref), **label_scores(ref['t'], ref, sections)}]).to_csv(
        root / 'reference_scores.csv', index=False)

    pd.set_option('display.width', 250)
    done = table[table['error'] == ''] if 'error' in table else table
    print(f'{len(done)} reruns scored, {int((table.get("error", pd.Series(dtype=str)) != "").sum())} failed')
    print('reference rho vs its own labels:', pd.read_csv(root / 'reference_scores.csv').round(3).to_dict('records')[0])
    cols = ['rho_vs_reference', 'rho_vs_full_rerun', 'rho_labels_pooled', *[f'rho_labels_{s}' for s in sections], 'segment_order_rho']
    print(done.groupby(['condition', 'fraction'])[cols].median().round(3).to_string())
    print('flipped orientations:', done.loc[done['orientation_flipped'].astype(bool), 'job'].tolist())
    print('PT limb reversed (median t of PT-S1 above PT-S3):', done.loc[done['pt_limb_reversed'].astype(bool), 'job'].tolist())
    full = done[done['fraction'] == 1.0]
    if len(full):
        print('100 % arms (reproduction check):')
        print(full[['condition', 'seed', 'n_units', 'n_shared', 'rho_vs_reference', 'max_abs_dt', 'fingerprint_pass1',
                    'fingerprint_pass2']].to_string(index=False))
    if halves:
        print(pd.DataFrame(halves).round(4).to_string(index=False))


def main(argv=None):
    args = parse_args(argv)
    root = Path(args.results_root) / OUT / args.run
    out = root / 'reruns'
    out.mkdir(parents=True, exist_ok=True)
    if args.job is not None:
        run_job(args, args.job, out)
    elif args.score:
        score(args, out, root)
    else:
        launch(args, out)


if __name__ == '__main__':
    main()
