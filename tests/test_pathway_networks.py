import numpy as np
import pandas as pd
import pytest

from pseudospace import pathway_networks as pn

pytest.importorskip('igraph')


def test_segment_cuts_minimise_misassigned_units():
    rng = np.random.default_rng(0)
    x = rng.uniform(0, 1, 3000)
    seg = np.where(x < .3, 'S1', np.where(x < .7, 'S2', 'S3'))
    flip = rng.choice(len(x), 60, replace=False)
    seg[flip] = rng.choice(['S1', 'S2', 'S3'], 60)
    c12, c23 = pn.segment_cuts(x, seg)
    assert abs(c12 - .3) < .03 and abs(c23 - .7) < .03
    # a clean split: every candidate between the last S1 and the first S2 ties, so the cut is their midpoint
    c12, c23 = pn.segment_cuts([0, 1, 2, 5, 6, 9], ['S1', 'S1', 'S1', 'S2', 'S2', 'S3'])
    assert c12 == pytest.approx(3.5) and c23 == pytest.approx(7.5)
    with pytest.raises(ValueError):
        pn.segment_cuts([0, 1, 2], ['S1', 'S1', 'S1'])


def test_gene_shapes_extremes_sides_and_lists():
    grid = np.linspace(0, 1, 101)
    cuts = (.3, .7)
    rising = np.log(2) * (grid - .5) * 2                      # +1 log2 at the end, -1 at the start
    trough = -np.log(2) * np.exp(-((grid - .5) / .05) ** 2)   # relatively low in S2
    falling = -rising + 3.                                    # a level offset must not matter
    delta = np.vstack([rising, trough, falling, rising])
    shape, table = pn.gene_shapes(delta, [True, True, True, False], grid, cuts)
    assert np.allclose(shape[0], (rising - rising.mean()) / np.log(2)) and np.isnan(shape[3]).all()
    assert table.loc[0, ['direction', 'segment', 'side']].tolist() == ['up', 'S3', 1]
    assert table.loc[1, ['direction', 'segment', 'side']].tolist() == ['down', 'S2', 0]
    assert table.loc[2, ['direction', 'segment', 'side']].tolist() == ['up', 'S1', -1]
    assert table.half_range[0] == pytest.approx(1., abs=1e-9) and not table.usable[3]
    lists = pn.position_split_lists(table, [True, True, True, True])
    assert len(lists) == 6
    stacked = np.vstack(list(lists.values()))
    assert stacked.sum(axis=0).tolist() == [1, 1, 1, 0]      # every usable selected gene in exactly one list
    assert lists[('up', 'S3')][0] and lists[('down', 'S2')][1] and lists[('up', 'S1')][2]


def test_driver_jaccard_and_position_score():
    names, J = pn.driver_jaccard({'a': {'x', 'y'}, 'b': {'y', 'z'}, 'c': set()})
    assert names == ['a', 'b', 'c']
    assert J[0, 1] == pytest.approx(1 / 3) and J[0, 0] == 1 and J[2, 2] == 0 and J[0, 2] == 0
    assert pn.position_score([1, 1, 1, 1]) == pytest.approx(2.)
    assert pn.position_score([]) == 0 and pn.score_class(-1.) == 'toward S1' and pn.score_class(.5) == 'mixed'


def test_leiden_programs_groups_cliques_keeps_isolates_and_is_deterministic():
    names = [f'c{i}' for i in range(9)]
    W = np.zeros((9, 9))
    for block in (range(0, 4), range(4, 8)):
        for i in block:
            for j in block:
                W[i, j] = .6
    W[3, 4] = W[4, 3] = .12                                    # a weak bridge above threshold
    W[0, 8] = W[8, 0] = .05                                    # below threshold: c8 is an isolate
    p = pd.Series(np.linspace(.01, .09, 9), index=names)
    labels = pn.leiden_programs(names, W, threshold=.1, best_p=p)
    assert labels[[f'c{i}' for i in range(4)]].nunique() == 1 and labels[[f'c{i}' for i in range(4, 8)]].nunique() == 1
    assert labels['c0'] != labels['c4'] and labels['c8'] == 3   # two programs of four, then the isolate
    assert labels['c0'] == 1                                   # tie on size broken by the smaller best p
    again = pn.leiden_programs(names, W, threshold=.1, best_p=p)
    assert labels.equals(again)
    alone = pn.leiden_programs(names, W, threshold=.9)
    assert sorted(alone) == list(range(1, 10))


def test_program_stats_coherence_replication_and_null():
    rng = np.random.default_rng(1)
    grid = np.linspace(0, 1, 50)
    late, early = grid - .5, .5 - grid
    shape = np.vstack([late + rng.normal(0, .05, 50) for _ in range(30)] + [early + rng.normal(0, .05, 50) for _ in range(30)])
    shape -= shape.mean(axis=1, keepdims=True)
    drivers = {f'L{i}': list(range(i * 5, i * 5 + 5)) for i in range(6)}
    drivers.update({f'E{i}': list(range(30 + i * 5, 30 + i * 5 + 5)) for i in range(6)})
    spec = (shape + rng.normal(0, .05, shape.shape), shape + rng.normal(0, .05, shape.shape))
    stats = pn.ProgramStats(shape=shape, drivers=drivers, specimen_shapes=spec)
    by_shape = pd.Series({c: c[0] for c in drivers})
    table, profiles = stats.per_program(by_shape)
    assert (table.coherence > .95).all() and (table.replication > .95).all() and (table.genes == 30).all()
    coh, rep = stats.summary(by_shape)
    null = stats.null(by_shape, n=200, seed=3)
    assert pn.null_z(coh, null[:, 0]) > 3
    assert pn.null_z(1., np.ones(10)) != pn.null_z(1., np.ones(10))       # no spread -> NaN
    one = pd.Series({c: i for i, c in enumerate(drivers)})
    single, _ = stats.per_program(one)
    assert len(single) == 12 and (single.calls == 1).all()


def test_hierarchies_and_ancestor_naming(tmp_path):
    rel = tmp_path / 'rel.txt'
    rel.write_text('R-HSA-1\tR-HSA-2\nR-HSA-2\tR-HSA-3\nR-HSA-2\tR-HSA-4\nR-HSA-2\tR-HSA-6\nR-HSA-1\tR-HSA-5\n'
                   'R-MMU-1\tR-MMU-2\n')
    paths = tmp_path / 'paths.txt'                             # 100 human pathways, most of them unrelated
    paths.write_text(''.join(f'R-HSA-{i}\tName {i}\tHomo sapiens\n' for i in range(1, 101)) + 'R-MMU-1\tX\tMus musculus\n')
    h = pn.reactome_hierarchy(rel, paths)
    assert h.size == 100 and h.ancestors('R-HSA-3') == {'R-HSA-2', 'R-HSA-1'} and h.descendants('R-HSA-1') == 5
    hs = {'R-HSA-': h}
    picks, covered = pn.name_program(['R-HSA-3', 'R-HSA-4', 'R-HSA-6'], hs)
    assert picks == ['R-HSA-2'] and covered == 3               # 3 x log(100/4) beats a leaf's 1 x log(100) and the root
    picks, _ = pn.name_program(['R-HSA-3', 'R-HSA-4', 'R-HSA-6'], hs, max_descendant_share=.02)
    assert picks == ['R-HSA-3', 'R-HSA-4']                     # the parent is too large: leaves, at most two picks
    picks, covered = pn.name_program(['R-HSA-3', None, None, None, None, None, None, None], hs)
    assert picks == [] and covered == 0                        # one of eight calls is below the 15% share
    keg = tmp_path / 'br.keg'
    keg.write_text('+C\tMap number\n!\nAMetabolism\nB  Lipid metabolism\nC    00061  Fatty acid biosynthesis\n'
                   'C    00062  Fatty acid elongation\nACellular Processes\nB  Transport and catabolism\n'
                   'C    04140  Autophagy - animal\n' + ''.join(f'C    05{i:03d}  Other {i}\n' for i in range(20)))
    k, lookup = pn.kegg_hierarchy(keg)
    assert k.size == 23 and pn.kegg_node('Fatty acid biosynthesis', lookup) == 'K:map00061'
    assert pn.kegg_node('Autophagy', lookup) == 'K:map04140' and pn.kegg_node('Nothing here', lookup) is None
    assert k.ancestors('K:map00062') == {'K:B:Lipid metabolism', 'K:A:Metabolism'}
    picks, _ = pn.name_program(['K:map00061', 'K:map00062'], {'K:': k}, max_descendant_share=1.)
    assert picks == ['K:B:Lipid metabolism']


def test_dominance():
    directions = {'pass': 1, 'jaccard': 1, 'exposure': -1}
    tol = {'jaccard': .05, 'exposure': .05}
    assert pn.dominance({'pass': 1, 'jaccard': .8, 'exposure': .2}, {'pass': 1, 'jaccard': .5, 'exposure': .2}, directions, tol) == 'better'
    assert pn.dominance({'pass': 0, 'jaccard': .8, 'exposure': .2}, {'pass': 1, 'jaccard': .5, 'exposure': .2}, directions, tol) == 'mixed'
    assert pn.dominance({'pass': 1, 'jaccard': .5, 'exposure': .6}, {'pass': 1, 'jaccard': .52, 'exposure': .2}, directions, tol) == 'worse'
    assert pn.dominance({'pass': 1, 'jaccard': np.nan}, {'pass': 1, 'jaccard': .5}, directions, tol) == 'equal'
    assert pn.dominance({}, {}, directions, tol) == 'undefined'
