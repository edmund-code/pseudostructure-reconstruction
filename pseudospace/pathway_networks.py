"""Pathway networks along the PT: scoranet's over-representation networks with positional additions (notebook 93).

Protocol: ``docs/results/pt-pathway-networks.md``. scoranet (GitHub mikemartinez99/scoranet, v1.0.0) links enriched
terms whose DE genes overlap (Jaccard >= 0.1), groups them with Leiden (modularity, resolution 1) and names each group
by GO ancestors. Here the same network is built on any method's pathway calls, with the PT axis added:

- ``segment_cuts``: the S1|S2 and S2|S3 cuts that misassign the fewest labelled units.
- ``gene_shapes``: each gene's positional difference delta(x) (centred, log2), its half-range and its extreme.
- ``position_split_lists``: scoranet's up/down split extended to up/down x segment (row N's six lists).
- ``driver_jaccard``, ``leiden_programs``: the network and its programs (isolates kept as one-call programs).
- ``position_score``: scoranet's direction score with late/early in place of up/down.
- ``ProgramStats``: program coherence R^2 and specimen replication r, with a null of random partitions of the calls.
- ``Hierarchy``, ``reactome_hierarchy``, ``kegg_hierarchy``, ``name_program``: scoranet's ancestor naming on Reactome and
  KEGG.
- ``dominance``: the protocol's verdict rule.
"""
from __future__ import annotations

import random
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

SEGMENTS = ('S1', 'S2', 'S3')
LN2 = np.log(2.)


def segment_cuts(position, segment):
    """(c12, c23): each cut minimises the units on the wrong side (S1 against S2 + S3; S1 + S2 against S3). Candidate
    cuts are the midpoints between consecutive distinct positions; ties go to the midpoint of the tied candidates."""
    x = np.asarray(position, float)
    s = np.asarray(segment).astype(str)
    order = np.argsort(x, kind='stable')
    x, s = x[order], s[order]
    values = np.unique(x)
    if len(values) < 2:
        raise ValueError('Need at least two distinct positions.')
    candidates = (values[:-1] + values[1:]) / 2
    cuts = []
    for early in ({'S1'}, {'S1', 'S2'}):
        is_early = np.isin(s, list(early))
        if is_early.all() or not is_early.any():
            raise ValueError('Each cut needs units on both sides.')
        left = np.searchsorted(x, candidates, side='right')            # units below each candidate
        late_below = np.concatenate([[0], np.cumsum(~is_early)])[left]
        early_above = is_early.sum() - np.concatenate([[0], np.cumsum(is_early)])[left]
        errors = late_below + early_above
        tied = candidates[errors == errors.min()]
        cuts.append(float((tied.min() + tied.max()) / 2))
    if not cuts[0] < cuts[1]:
        raise ValueError(f'Cuts are not ordered: {cuts}')
    return tuple(cuts)


def grid_segments(grid, cuts):
    """Segment of each grid point: S1 below c12, S3 at or above c23, S2 between."""
    g = np.asarray(grid, float)
    return np.where(g < cuts[0], 'S1', np.where(g >= cuts[1], 'S3', 'S2'))


def gene_shapes(delta, usable, grid, cuts):
    """delta: genes x grid natural-log case - reference differences. Returns (shape, table): shape = genes x grid
    (Delta - its grid mean) / ln 2 (NaN rows for unusable genes); table per gene: usable, half_range, direction ('up'
    when max >= -min, else 'down'), extreme (grid index), segment (of the extreme) and side (+1 late: up in S3 or down
    in S1; -1 early: up in S1 or down in S3; 0 otherwise)."""
    D = np.asarray(delta, float)
    ok = np.asarray(usable, bool) & np.isfinite(D).all(axis=1)
    shape = np.full(D.shape, np.nan)
    shape[ok] = (D[ok] - D[ok].mean(axis=1, keepdims=True)) / LN2
    hi, lo = np.nanmax(np.where(ok[:, None], shape, 0.), axis=1), np.nanmin(np.where(ok[:, None], shape, 0.), axis=1)
    up = hi >= -lo
    extreme = np.where(up, np.nanargmax(np.where(ok[:, None], shape, -np.inf), axis=1),
                       np.nanargmin(np.where(ok[:, None], shape, np.inf), axis=1))
    seg = grid_segments(grid, cuts)[extreme]
    side = np.where((up & (seg == 'S3')) | (~up & (seg == 'S1')), 1, np.where((up & (seg == 'S1')) | (~up & (seg == 'S3')), -1, 0))
    table = pd.DataFrame({'usable': ok, 'half_range': np.where(ok, (hi - lo) / 2, np.nan),
                          'direction': np.where(ok, np.where(up, 'up', 'down'), ''), 'extreme': np.where(ok, extreme, -1),
                          'segment': np.where(ok, seg, ''), 'side': np.where(ok, side, 0)})
    return shape, table


def position_split_lists(table, selected):
    """Row N's six lists: {(direction, segment): boolean vector}, over the genes of ``table`` (``gene_shapes``)."""
    sel = np.asarray(selected, bool) & table.usable.to_numpy()
    return {(d, s): sel & table.direction.eq(d).to_numpy() & table.segment.eq(s).to_numpy()
            for d in ('up', 'down') for s in SEGMENTS}


def driver_jaccard(drivers):
    """drivers: {call: set of genes}. Returns (names, Jaccard matrix) in ``drivers`` order (0 where a union is empty)."""
    names = list(drivers)
    genes = sorted(set().union(*drivers.values())) if drivers else []
    index = {g: i for i, g in enumerate(genes)}
    rows = [i for i, n in enumerate(names) for _ in drivers[n]]
    cols = [index[g] for n in names for g in drivers[n]]
    M = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(names), len(genes)))
    inter = (M @ M.T).toarray()
    size = np.asarray(M.sum(axis=1)).ravel()
    union = size[:, None] + size[None, :] - inter
    J = np.divide(inter, union, out=np.zeros_like(inter), where=union > 0)
    return names, J


def leiden_programs(names, weights, *, threshold, resolution=1., seed=1234, best_p=None):
    """scoranet's grouping: an undirected graph with an edge where ``weights`` >= ``threshold`` (weight = the value),
    Leiden on modularity (igraph ``community_leiden``, beta 0.01, iterated until stable, Python ``random`` seeded
    before the call). Isolates are one-call programs. Programs are numbered 1.. by size, ties by the smallest
    ``best_p`` (then first appearance). Returns a Series call -> program."""
    import igraph as ig
    names = list(names)
    n = len(names)
    if n == 0:
        return pd.Series(dtype=int, name='program')
    W = np.asarray(weights, float)
    iu, ju = np.triu_indices(n, 1)
    keep = (W[iu, ju] >= threshold) & (W[iu, ju] > 0)
    graph = ig.Graph(n=n, edges=list(zip(iu[keep].tolist(), ju[keep].tolist())))
    graph.es['weight'] = W[iu, ju][keep].tolist()
    if graph.ecount():
        random.seed(seed)
        membership = graph.community_leiden(objective_function='modularity', weights='weight', resolution=resolution,
                                             beta=.01, n_iterations=-1).membership
    else:
        membership = list(range(n))
    labels = pd.Series(membership, index=names)
    p = pd.Series(np.ones(n) if best_p is None else np.asarray(best_p, float), index=names)
    first = pd.Series(np.arange(n), index=names)
    summary = pd.DataFrame({'size': labels.value_counts(), 'p': p.groupby(labels).min(), 'first': first.groupby(labels).min()})
    order = summary.sort_values(['size', 'p', 'first'], ascending=[False, True, True]).index
    relabel = {old: new for new, old in enumerate(order, start=1)}
    return labels.map(relabel).rename('program')


def profile_correlation(profiles):
    """Pearson correlation between rows (calls x grid); 0 for a constant row."""
    P = np.asarray(profiles, float)
    C = P - P.mean(axis=1, keepdims=True)
    norm = np.sqrt((C ** 2).sum(axis=1))
    ok = norm > 0
    C[ok] /= norm[ok, None]
    C[~ok] = 0.
    return np.clip(C @ C.T, -1., 1.)


def position_score(sides):
    """scoranet's direction score with late (+1) and early (-1) genes: (n_late - n_early) / sqrt(n); 0 when empty."""
    s = np.asarray(sides)
    return float((np.sum(s > 0) - np.sum(s < 0)) / np.sqrt(len(s))) if len(s) else 0.


def score_class(score, threshold=1.):
    return 'toward S3' if score >= threshold else 'toward S1' if score <= -threshold else 'mixed'


@dataclass
class ProgramStats:
    """Coherence and replication of programs of calls, with a null of random partitions.

    shape: genes x grid delta (log2, centred; NaN rows unusable); drivers: {call: gene indices}; specimen_shapes: two
    genes x grid arrays (each case specimen against the reference mean, centred, log2; NaN rows unusable)."""
    shape: np.ndarray
    drivers: dict
    specimen_shapes: tuple = ()
    _calls: list = field(init=False, repr=False)
    _D: sparse.csr_matrix = field(init=False, repr=False)

    def __post_init__(self):
        self._calls = list(self.drivers)
        G = self.shape.shape[0]
        usable = np.isfinite(self.shape).all(axis=1)
        rows = [i for i, c in enumerate(self._calls) for g in self.drivers[c] if usable[g]]
        cols = [g for c in self._calls for g in self.drivers[c] if usable[g]]
        self._D = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(self._calls), G))
        self._X = np.where(usable[:, None], self.shape, 0.)
        self._ss = (self._X ** 2).sum(axis=1)
        if self.specimen_shapes:
            both = np.all([np.isfinite(s).all(axis=1) for s in self.specimen_shapes], axis=0) & usable
            self._rep = [np.where(both[:, None], s, 0.) for s in self.specimen_shapes]
            self._rep_ok = both.astype(float)

    def _unions(self, labels):
        """Programs x genes 0/1 union of the calls' driver genes, for a program label per call."""
        codes, uniq = pd.factorize(pd.Series(labels).reindex(self._calls).to_numpy())
        P = sparse.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniq), len(codes)))
        U = (P @ self._D)
        U.data[:] = 1.
        return uniq, np.bincount(codes, minlength=len(uniq)), U

    def per_program(self, labels):
        """DataFrame per program: calls, genes, coherence R^2, replication r, the profile half-range."""
        uniq, sizes, U = self._unions(labels)
        n = np.asarray(U.sum(axis=1)).ravel()
        S = np.asarray(U @ self._X)
        mean = np.divide(S, n[:, None], out=np.zeros_like(S), where=n[:, None] > 0)
        tot = np.asarray(U @ self._ss).ravel()
        r2 = np.divide(n * (mean ** 2).sum(axis=1), tot, out=np.full(len(n), np.nan), where=tot > 0)
        out = pd.DataFrame({'calls': sizes, 'genes': n.astype(int), 'coherence': r2,
                            'half_range': (mean.max(axis=1) - mean.min(axis=1)) / 2}, index=uniq)
        if self.specimen_shapes:
            m = np.asarray(U @ self._rep_ok).ravel()
            prof = [np.divide(np.asarray(U @ s), m[:, None], out=np.full(S.shape, np.nan), where=m[:, None] > 0) for s in self._rep]
            out['replication'] = [_pearson(a, b) for a, b in zip(*prof)]
            out['replication genes'] = m.astype(int)
        return out, mean

    def summary(self, labels):
        """(call-weighted mean coherence, call-weighted mean replication)."""
        t, _ = self.per_program(labels)
        coh = _weighted(t.coherence, t.calls)
        rep = _weighted(t.replication, t.calls) if 'replication' in t else np.nan
        return coh, rep

    def null(self, labels, *, n=1000, seed=93):
        """(coherence, replication) of ``n`` random partitions of the calls with the same program sizes."""
        lab = pd.Series(labels).reindex(self._calls).to_numpy()
        rng = np.random.default_rng(seed)
        out = np.empty((n, 2))
        for k in range(n):
            out[k] = self.summary(pd.Series(rng.permutation(lab), index=self._calls))
        return out


def _weighted(values, weights):
    v, w = np.asarray(values, float), np.asarray(weights, float)
    ok = np.isfinite(v)
    return float((v[ok] * w[ok]).sum() / w[ok].sum()) if ok.any() else np.nan


def _pearson(a, b):
    if not (np.isfinite(a).all() and np.isfinite(b).all()) or np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def null_z(observed, null):
    """(observed - null mean) / null SD; NaN when the null has no spread."""
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    if not np.isfinite(observed) or len(null) < 2 or np.std(null, ddof=1) == 0:
        return np.nan
    return float((observed - null.mean()) / np.std(null, ddof=1))


@dataclass
class Hierarchy:
    """A pathway hierarchy: parents (child -> set of parents), display names, and the number of pathways it holds.
    Descendants and ancestors are transitive."""
    parents: dict
    names: dict
    size: int

    def __post_init__(self):
        children = {}
        for child, ps in self.parents.items():
            for p in ps:
                children.setdefault(p, set()).add(child)
        self._children = children
        self._anc, self._desc = {}, {}

    def ancestors(self, node):
        if node not in self._anc:
            seen, stack = set(), list(self.parents.get(node, ()))
            while stack:
                p = stack.pop()
                if p not in seen:
                    seen.add(p)
                    stack.extend(self.parents.get(p, ()))
            self._anc[node] = seen
        return self._anc[node]

    def descendants(self, node):
        if node not in self._desc:
            seen, stack = set(), list(self._children.get(node, ()))
            while stack:
                c = stack.pop()
                if c not in seen:
                    seen.add(c)
                    stack.extend(self._children.get(c, ()))
            self._desc[node] = len(seen)
        return self._desc[node]


def reactome_hierarchy(relation_file, pathways_file, *, species='Homo sapiens', prefix='R-HSA-'):
    """Reactome's pathway hierarchy for one species (``ReactomePathwaysRelation.txt``: parent<TAB>child;
    ``ReactomePathways.txt``: id<TAB>name<TAB>species). Nodes are stable ids."""
    names = {}
    for line in Path(pathways_file).read_text().splitlines():
        parts = line.split('\t')
        if len(parts) >= 3 and parts[2] == species:
            names[parts[0]] = parts[1]
    parents = {}
    for line in Path(relation_file).read_text().splitlines():
        parts = line.split('\t')
        if len(parts) == 2 and parts[0].startswith(prefix) and parts[1].startswith(prefix):
            parents.setdefault(parts[1], set()).add(parts[0])
    return Hierarchy(parents=parents, names=names, size=len(names))


def _norm(text):
    return re.sub(r'[^a-z0-9]', '', text.lower())


def kegg_hierarchy(keg_file):
    """KEGG's pathway map hierarchy (``br08901``): maps (C) under subcategories (B) under categories (A). Returns
    (Hierarchy, lookup): nodes 'K:A:<name>', 'K:B:<name>' and 'K:map<number>'; ``lookup`` maps each map name to its
    node. The hierarchy size is the number of maps."""
    parents, names, lookup = {}, {}, {}
    a = b = None
    for line in Path(keg_file).read_text().splitlines():
        if line.startswith('A'):
            a = 'K:A:' + line[1:].strip()
            names[a] = line[1:].strip()
        elif line.startswith('B ') and line[1:].strip():
            b = 'K:B:' + line[1:].strip()
            names[b] = line[1:].strip()
            parents[b] = {a}
        elif line.startswith('C '):
            m = re.match(r'C\s+(\d+)\s+(.*)$', line)
            if m:
                node, name = 'K:map' + m.group(1), re.sub(r'\s*\[.*\]$', '', m.group(2)).strip()
                names[node] = name
                parents[node] = {b}
                lookup.setdefault(name, node)
    n_maps = sum(1 for k in names if k.startswith('K:map'))
    return Hierarchy(parents=parents, names=names, size=n_maps), lookup


def _strip(name):
    return _norm(re.sub(r'\s*-\s*animal\b', '', re.sub(r'\([^)]*\)', '', name)))


def kegg_node(name, lookup):
    """A KEGG set name's map node: exact normalised match, then both names without parenthetical text and
    '- animal'; None when neither matches."""
    for key in (_norm, _strip):
        table = {}
        for k, node in lookup.items():
            table.setdefault(key(k), node)
        if key(name) in table:
            return table[key(name)]
    return None


def name_program(nodes, hierarchies, *, min_share=.15, max_descendant_share=.04, max_picks=2):
    """scoranet's ancestor rule. nodes: one hierarchy node (or None) per call of the program; hierarchies: {prefix:
    Hierarchy} with node ids starting with their prefix. Each call covers its node and the node's ancestors; a
    candidate has at most ``max_descendant_share`` of its hierarchy as descendants; IC = log(size / (descendants + 1)).
    Up to ``max_picks`` greedy picks, each maximising (calls newly covered) x IC and newly covering at least
    ``min_share`` of all the program's calls; ties to the lower id. Returns (pick ids, calls covered)."""
    n = len(nodes)

    def hier(node):
        return next(h for p, h in hierarchies.items() if node.startswith(p))

    covered_by = [({nd} | hier(nd).ancestors(nd)) if nd is not None else set() for nd in nodes]
    candidates = sorted(set().union(*covered_by)) if covered_by else []
    picks, uncovered = [], np.ones(n, bool)
    if not candidates:
        return picks, 0
    covers = np.array([[c in s for c in candidates] for s in covered_by], bool).reshape(n, len(candidates))
    size = np.array([hier(c).descendants(c) for c in candidates], float)
    total = np.array([hier(c).size for c in candidates], float)
    info = np.log(total / (size + 1))
    eligible = size <= max_descendant_share * total
    for _ in range(max_picks):
        gain = covers[uncovered].sum(axis=0)
        ok = eligible & (gain / n >= min_share) & ~np.isin(candidates, picks)
        if not ok.any():
            break
        score = np.where(ok, gain * info, -np.inf)
        best = int(np.argmax(score))                     # first maximum: candidates are sorted, so the lower id
        picks.append(candidates[best])
        uncovered &= ~covers[:, best]
        if not uncovered.any():
            break
    return picks, int((~uncovered).sum())


def dominance(a, b, directions, tolerances):
    """Verdict of a against b over criteria (dicts keyed alike). directions: +1 higher is better, -1 lower is better
    (booleans compare as 0/1); a criterion missing or NaN on either side is skipped. Returns 'better', 'worse',
    'equal' (all within tolerance) or 'mixed'."""
    signs = []
    for k, d in directions.items():
        x, y = a.get(k), b.get(k)
        if x is None or y is None:
            continue
        x, y = float(x), float(y)
        if not (np.isfinite(x) and np.isfinite(y)):
            continue
        diff = d * (x - y)
        signs.append(0 if abs(x - y) <= tolerances.get(k, 0.) else int(np.sign(diff)))
    if not signs:
        return 'undefined'
    if all(s == 0 for s in signs):
        return 'equal'
    if all(s >= 0 for s in signs):
        return 'better'
    if all(s <= 0 for s in signs):
        return 'worse'
    return 'mixed'
