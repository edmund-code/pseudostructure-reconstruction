"""Spatially blocked delete-half jackknife of the human-vs-mouse pseudospace.

Each draw keeps a random half of the spatial blocks of every section (stage-02 QC units of the reference run),
reruns the frozen stage 03 (two-pass Harmony, Leiden, automatic labels) and stage 04 (elastic curve), and stores
every kept unit's position and labels. `--combine` puts all draws on the PT scale and writes the units x draws
position table, per-unit summaries and a per-draw table. Protocol: docs/results/pt-position-uncertainty.md.

Reads private data and writes under --results-root (outside Git). Resumable: finished draws are skipped.

  python analysis/scripts/position_half_sampling.py --draws 100 --jobs 8
  python analysis/scripts/position_half_sampling.py --combine
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


def _repo() -> Path:
    return next(p for p in [Path(__file__).resolve(), *Path(__file__).resolve().parents] if (p / 'pseudospace').is_dir())


REPO = _repo()
sys.path.insert(0, str(REPO))


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--data-root', default=os.environ.get('PSEUDOSPACE_DATA_ROOT', str(REPO / 'data')))
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--config', default=str(REPO / 'analysis/package_configs/human_vs_mouse'))
    p.add_argument('--run', default='package_pseudospace/runs/hvm', help='reference run, relative to the results root')
    p.add_argument('--out', default='package_pseudospace/position_half_sampling', help='relative to the results root')
    p.add_argument('--draws', type=int, default=100)
    p.add_argument('--start', type=int, default=0)
    p.add_argument('--jobs', type=int, default=8)
    p.add_argument('--threads', type=int, default=4, help='BLAS/OpenMP threads per draw; FIXED by the protocol: Leiden partitions change with the thread count')
    p.add_argument('--per-block', type=float, default=25.0)
    p.add_argument('--tmp', default='/dev/shm' if Path('/dev/shm').is_dir() else '', help='scratch for per-draw stage files')
    p.add_argument('--draw', type=int, help=argparse.SUPPRESS)       # internal: run one draw in this process
    p.add_argument('--combine', action='store_true')
    p.add_argument('--evaluate', action='store_true', help="the protocol's pre-declared statistics (after --combine)")
    return p.parse_args(argv)


def run_draw(args, seed: int, out: Path) -> None:
    import anndata as ad
    from pseudospace_reconstruction.config import RunPaths, load_config
    from pseudospace_reconstruction.integrate import run_two_pass
    from pseudospace_reconstruction.pipeline import ensure_r_home, run_stage04, stage03_arguments
    from pseudospace.position_uncertainty import half_sample, spatial_blocks

    t0 = time.time()
    ensure_r_home()
    config = load_config(Path(args.config), Path(args.data_root))
    stage03 = stage03_arguments(config)
    qc = ad.read_h5ad(Path(args.results_root) / args.run / '02_qc/units_qc.h5ad')
    section = qc.obs['sample'].astype(str).to_numpy()
    blocks = spatial_blocks(qc.obs[['centroid_x', 'centroid_y']].to_numpy(float), section, args.per_block)
    keep = half_sample(blocks, section, np.random.default_rng(seed))
    tmp = Path(args.tmp) / f'position_half_sampling_{seed:03d}' if args.tmp else out / f'tmp_draw_{seed:03d}'
    try:
        res = run_two_pass(qc[keep].copy(), stage03['panels'], stage03['params'])
        paths = RunPaths(tmp)
        paths.stage('03', create=True)
        p2 = res.pass2   # stage 04 reads obs, the embedding and the lognorm layer (all genes: the marker-axis score
        #                  draws its control genes from every gene), so the count matrices are not written
        ad.AnnData(obs=p2.obs, var=p2.var, obsm={'X_harmony': p2.obsm['X_harmony']},
                   layers={'lognorm': p2.layers['lognorm']}).write(paths.pass2)
        paths.manifest('03').write_text(json.dumps({'stage': '03', 'half_sampling_seed': seed}))
        run_stage04(config, paths)
        ps = pd.read_parquet(tmp / '04_pseudospace/pseudospace.parquet')['pseudospace']
        obs = res.pass2.obs
        table = pd.DataFrame({'pseudospace': ps.reindex(obs.index).to_numpy(), 'sample': obs['sample'].astype(str),
                              'segment_class': obs['segment_class'].astype(str), 'coarse_class': obs['coarse_class'].astype(str)},
                             index=obs.index)
        table.to_parquet(out / f'draw_{seed:03d}.parquet')
        om = pd.read_csv(tmp / '04_pseudospace/order_metrics.csv')
        om = om[(om.coordinate == 'pseudospace') & (om.units == 'all')].iloc[0]
        meta = {'seed': seed, 'n_kept_qc_units': int(keep.sum()), 'n_qc_units': int(len(keep)),
                'n_blocks': int(len(np.unique(blocks))), 'pass1_fingerprint': dict(res.fingerprint_pass1),
                'segment_order_rho': float(om.segment_order_rho), 'n_inversions': int(om.n_inversions),
                'seconds': round(time.time() - t0)}
        (out / f'draw_{seed:03d}.json').write_text(json.dumps(meta, default=str))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def launch(args, out: Path) -> None:
    todo = [s for s in range(args.start, args.start + args.draws) if not (out / f'draw_{s:03d}.parquet').exists()]
    env = {**os.environ, **{k: str(args.threads) for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')}}
    base = [sys.executable, str(Path(__file__).resolve()), '--data-root', args.data_root, '--results-root', args.results_root,
            '--config', args.config, '--run', args.run, '--out', args.out, '--per-block', str(args.per_block), '--tmp', args.tmp]
    running: list[tuple[int, subprocess.Popen]] = []
    failed = []
    while todo or running:
        while todo and len(running) < args.jobs:
            s = todo.pop(0)
            log = open(out / f'draw_{s:03d}.log', 'w')
            running.append((s, subprocess.Popen([*base, '--draw', str(s)], env=env, stdout=log, stderr=subprocess.STDOUT)))
        time.sleep(5)
        for s, proc in list(running):
            if proc.poll() is not None:
                running.remove((s, proc))
                if proc.returncode:
                    failed.append(s)
                print(f'draw {s:03d} {"FAILED" if proc.returncode else "done"}', flush=True)
    print(f'finished; failed draws: {failed}', flush=True)


def combine(args, out: Path) -> None:
    from pseudospace.position_uncertainty import unit_summary

    cols, meta = {}, []
    for f in sorted(out.glob('draw_*.parquet')):
        d = pd.read_parquet(f)
        d = d[d.coarse_class == 'PT']
        cols[f.stem] = pd.Series(_pt_scale(d.pseudospace.to_numpy(float), d['sample'].to_numpy(), d.segment_class.to_numpy()),
                                 index=d.index)
        meta.append(json.loads(f.with_suffix('.json').read_text()))
    positions = pd.DataFrame(cols)
    positions.to_parquet(out / 'pt_positions_by_draw.parquet')
    unit_summary(positions).to_parquet(out / 'pt_position_summary.parquet')
    pd.DataFrame(meta).to_csv(out / 'draws.csv', index=False)
    print(f'{positions.shape[1]} draws, {positions.shape[0]} PT units; written to {out}')


FLIP_GENES = ('Mme', 'Dcxr', 'Ctsh')
MOUSE = ('Ctrl1A2', 'Ctrl1A4')
HUMAN = ('HUK1_COR1', 'HUK1_MED1')


def _pt_scale(raw, spec, seg):
    from pseudospace.package_inputs import common_support_mask
    mask = common_support_mask(raw, spec)
    x = np.where(mask, (raw - raw[mask].min()) / np.ptp(raw[mask]), np.nan)
    return 1 - x if np.nanmedian(x[seg == 'PT-S1']) > np.nanmedian(x[seg == 'PT-S3']) else x


def _statistics(x, spec, seg, cp10k) -> dict:
    """The protocol's pre-declared statistics on one axis (PT units only; x on the PT scale)."""
    from sklearn.linear_model import LogisticRegression

    def boundary(k, first, second):
        k = k & np.isfinite(x) & np.isin(seg, [first, second])
        m = LogisticRegression(C=1e4, class_weight='balanced').fit(x[k, None], (seg[k] == second).astype(int))
        return -m.intercept_[0] / m.coef_[0, 0]

    out = {}
    b = {(s, t): boundary(spec == s, *t) for s in MOUSE + HUMAN for t in (('PT-S1', 'PT-S2'), ('PT-S2', 'PT-S3'))}
    for t, name in ((('PT-S1', 'PT-S2'), 'S1->S2'), (('PT-S2', 'PT-S3'), 'S2->S3')):
        out[f'{name} boundary, human - mouse'] = np.mean([b[(s, t)] for s in HUMAN]) - np.mean([b[(s, t)] for s in MOUSE])
    cut12 = np.mean([b[(s, ('PT-S1', 'PT-S2'))] for s in MOUSE])
    cut23 = np.mean([b[(s, ('PT-S2', 'PT-S3'))] for s in MOUSE])
    for j, gene in enumerate(cp10k.columns):
        for group, names in (('mouse', MOUSE), ('human', HUMAN)):
            k = np.isin(spec, names) & np.isfinite(x)
            y = cp10k.to_numpy()[:, j]
            out[f'{gene} S3/S1 log2FC {group}'] = np.log2((y[k & (x > cut23)].mean() + .05) / (y[k & (x <= cut12)].mean() + .05))
    return out


def evaluate(args, out: Path) -> None:
    """Full-data estimates (reference run) with delete-half jackknife SEs over the draws."""
    import anndata as ad
    from scipy import sparse
    from pseudospace.position_uncertainty import call_share, interval

    ref = ad.read_h5ad(Path(args.results_root) / args.run / '03_integration/pass2.h5ad', backed='r')
    pt = (ref.obs.coarse_class.astype(str) == 'PT').to_numpy()
    # counts for EVERY QC unit: a draw labels a slightly different set of units PT than the reference does
    qc = ad.read_h5ad(Path(args.results_root) / args.run / '02_qc/units_qc.h5ad')
    counts = qc.layers['counts'] if 'counts' in qc.layers else qc.X
    counts = counts.tocsr() if sparse.issparse(counts) else sparse.csr_matrix(counts)
    nuclear = ~np.asarray(qc.var_names.str.startswith('mt-'))
    lib = np.asarray(counts[:, nuclear].sum(1)).ravel()
    cols = [qc.var_names.get_loc(g) for g in FLIP_GENES]
    cp10k = pd.DataFrame(np.asarray(counts[:, cols].todense()) / lib[:, None] * 1e4, index=qc.obs_names, columns=FLIP_GENES)

    def on_axis(table: pd.DataFrame) -> dict:
        table = table[table.coarse_class == 'PT']
        x = _pt_scale(table.pseudospace.to_numpy(float), table['sample'].to_numpy(), table.segment_class.to_numpy())
        return _statistics(x, table['sample'].to_numpy(), table.segment_class.to_numpy(), cp10k.reindex(table.index))

    ro = ref.obs[pt]
    full = on_axis(pd.DataFrame({'pseudospace': pd.read_parquet(Path(args.results_root) / args.run / '04_pseudospace/pseudospace.parquet')
                                 ['pseudospace'].reindex(ro.index), 'sample': ro['sample'].astype(str),
                                 'segment_class': ro.segment_class.astype(str), 'coarse_class': 'PT'}, index=ro.index))
    files = sorted(out.glob('draw_*.parquet'))
    draws = pd.DataFrame([on_axis(pd.read_parquet(f)) for f in files], index=[int(f.stem[-3:]) for f in files])
    draws.to_csv(out / 'statistics_by_draw.csv')
    rows = [{'statistic': k, **interval(v, draws[k].to_numpy())} for k, v in full.items()]
    for gene in FLIP_GENES:
        rows.append({'statistic': f'{gene} reverses direction (share of draws)',
                     'estimate': call_share(list(np.sign(draws[f'{gene} S3/S1 log2FC mouse']) != np.sign(draws[f'{gene} S3/S1 log2FC human'])),
                                            bool), 'n_draws': len(draws)})
    table = pd.DataFrame(rows)
    table.to_csv(out / 'predeclared_statistics.csv', index=False)
    summary = pd.read_parquet(out / 'pt_position_summary.parquet').join(
        pd.DataFrame({'species': np.where(ro['sample'].astype(str).str.startswith('HUK'), 'human', 'mouse')}, index=ro.index), how='inner')
    zone = pd.cut(summary['mean'], [-np.inf, .342, .612, np.inf], labels=['S1 zone', 'S2 zone', 'S3 zone'])
    se = summary.groupby(['species', zone], observed=True)['se'].median().unstack()
    se.to_csv(out / 'unit_se_by_species_zone.csv')
    pd.set_option('display.width', 200)
    print(table.round(3).to_string(index=False))
    print('\nmedian per-unit position SE (PT lengths):')
    print(se.round(3).to_string())

    # --- not pre-declared (added after seeing the per-unit SEs) --------------------------------------------------
    # 1 how much of a unit's spread is the whole draw's scale moving: align each draw to the reference by a
    #   least-squares affine map on the units they share, and recompute per-unit SEs on the residual positions
    # 2 the pre-declared statistics without the draws whose stage-04 order had inversions
    positions = pd.read_parquet(out / 'pt_positions_by_draw.parquet')
    ref_x = summary['mean'].reindex(positions.index) * np.nan
    ref_table = pd.read_parquet(Path(args.results_root) / args.run / '04_pseudospace/pseudospace.parquet')['pseudospace'].reindex(ro.index)
    ref_x = pd.Series(_pt_scale(ref_table.to_numpy(float), ro['sample'].astype(str).to_numpy(), ro.segment_class.astype(str).to_numpy()), index=ro.index)
    aligned = {}
    for c in positions:
        both = positions[c].notna() & ref_x.reindex(positions.index).notna()
        a, b = np.polyfit(positions.loc[both, c], ref_x.reindex(positions.index)[both], 1)
        aligned[c] = a * positions[c] + b
    aligned = pd.DataFrame(aligned)
    from pseudospace.position_uncertainty import unit_summary as _us
    s2 = _us(aligned).join(summary[['species']], how='inner')
    zone2 = pd.cut(ref_x.reindex(s2.index), [-np.inf, .342, .612, np.inf], labels=['S1 zone', 'S2 zone', 'S3 zone'])
    print('\n[post hoc] median per-unit SE after affine alignment of each draw to the reference (scale removed):')
    print(s2.groupby(['species', zone2], observed=True)['se'].median().unstack().round(3).to_string())
    meta = pd.read_csv(out / 'draws.csv').set_index('seed')
    good = [i for i, f in enumerate(sorted(out.glob('draw_*.parquet'))) if meta.loc[int(f.stem[-3:]), 'n_inversions'] == 0]
    rows2 = [{'statistic': k, **interval(v, draws[k].to_numpy()[good])} for k, v in full.items()]
    print(f'\n[post hoc] pre-declared statistics on the {len(good)} draws without order inversions:')
    print(pd.DataFrame(rows2).round(3).to_string(index=False))


def main(argv=None):
    args = parse_args(argv)
    out = Path(args.results_root) / args.out
    out.mkdir(parents=True, exist_ok=True)
    if args.draw is not None:
        run_draw(args, args.draw, out)
    elif args.combine:
        combine(args, out)
    elif args.evaluate:
        evaluate(args, out)
    else:
        launch(args, out)


if __name__ == '__main__':
    main()
