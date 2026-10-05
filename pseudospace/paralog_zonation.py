"""Ortholog families, probe specificity and positional profiles of paralogs (notebook 62).

The cross-species analyses use ``cross_species.build_one_to_one_ortholog_map``. Genes that map to
several partners, or to none, are outside that space. These helpers study them inside each species,
on a shared PT coordinate, without any ortholog mapping:

* ``hcop_orthogroups`` turns HCOP pairs into connected human-mouse components (orthogroups) under a
  declared support rule and classifies every gene as 1:1, 1:many, many:1, many:many or no consensus;
* ``transcript_table`` / ``probe_offtarget`` audit each Visium probe against its species
  transcriptome (bowtie2 alignments) so that paralog claims carry a probe status;
* ``weighted_curves`` fits per-specimen expression curves and nested F-tests for every gene at once;
* ``matched_pair_draws``, ``global_two_sided``, ``separation_call`` and ``nnls_union`` are the
  family-level statistics of stage 2.

Nothing here selects genes by a species contrast.
"""
from __future__ import annotations

import gzip
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import nnls
from scipy.stats import f as f_dist

CLASS_ORDER = ('1:1', '1:many', 'many:1', 'many:many', 'no_consensus', 'no_ortholog')
PSEUDOGENE = re.compile(r'pseudogene', re.IGNORECASE)


# ---------------------------------------------------------------------------------------------
# Ortholog families
# ---------------------------------------------------------------------------------------------

def support_count(value) -> int:
    """Number of databases in an HCOP ``support`` field (0 when missing)."""
    text = '' if value is None else str(value).strip()
    if text in ('', '-', 'nan', 'None'):
        return 0
    return len({part.strip() for part in text.split(',') if part.strip()})


def _key(prefix, ensembl, symbol):
    ensembl, symbol = str(ensembl).strip(), str(symbol).strip()
    if ensembl not in ('', '-', 'nan'):
        return f'{prefix}:{ensembl}'
    if symbol not in ('', '-', 'nan'):
        return f'{prefix}:sym:{symbol}'
    return None


def hcop_orthogroups(table, min_support=3):
    """Connected components of HCOP human-mouse pairs with ``support`` >= ``min_support``.

    ``table`` needs human_ensembl_gene, human_symbol, mouse_ensembl_gene, mouse_symbol and support.
    Genes are keyed by Ensembl id (``H:ENSG...``/``M:ENSMUSG...``), by symbol when the id is
    missing. Returns one row per gene seen in any HCOP row with columns key, species, ensembl,
    symbol, max_support, component (-1 when no kept pair), n_human, n_mouse and ortholog_class
    (``CLASS_ORDER`` without ``no_ortholog``, which only applies to genes absent from HCOP).
    """
    needed = {'human_ensembl_gene', 'human_symbol', 'mouse_ensembl_gene', 'mouse_symbol', 'support'}
    missing = needed - set(table.columns)
    if missing:
        raise ValueError(f'HCOP table lacks {sorted(missing)}')
    rows = table.copy()
    rows['h'] = [_key('H', e, s) for e, s in zip(rows.human_ensembl_gene, rows.human_symbol)]
    rows['m'] = [_key('M', e, s) for e, s in zip(rows.mouse_ensembl_gene, rows.mouse_symbol)]
    rows = rows.dropna(subset=['h', 'm'])
    rows['n_support'] = rows.support.map(support_count)
    pairs = rows.groupby(['h', 'm'], as_index=False).n_support.max()

    genes = pd.concat([
        rows[['h', 'human_ensembl_gene', 'human_symbol', 'n_support']].set_axis(
            ['key', 'ensembl', 'symbol', 'n_support'], axis=1).assign(species='human'),
        rows[['m', 'mouse_ensembl_gene', 'mouse_symbol', 'n_support']].set_axis(
            ['key', 'ensembl', 'symbol', 'n_support'], axis=1).assign(species='mouse')])
    genes = (genes.sort_values('n_support', ascending=False)
             .groupby('key', as_index=False)
             .agg(species=('species', 'first'), ensembl=('ensembl', 'first'), symbol=('symbol', 'first'),
                  max_support=('n_support', 'max')))

    kept = pairs[pairs.n_support >= min_support]
    parent = {k: k for k in genes.key}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for h, m in zip(kept.h, kept.m):
        rh, rm = find(h), find(m)
        if rh != rm:
            parent[max(rh, rm)] = min(rh, rm)
    linked = set(kept.h) | set(kept.m)
    roots = pd.Series({k: find(k) for k in genes.key})
    genes['root'] = genes.key.map(roots).where(genes.key.isin(linked))
    sizes = genes.dropna(subset=['root']).groupby(['root', 'species']).size().unstack(fill_value=0)
    sizes = sizes.reindex(columns=['human', 'mouse'], fill_value=0)
    order = {root: i for i, root in enumerate(sorted(sizes.index))}
    genes['component'] = genes.root.map(order).fillna(-1).astype(int)
    genes['n_human'] = genes.root.map(sizes.human).fillna(0).astype(int)
    genes['n_mouse'] = genes.root.map(sizes.mouse).fillna(0).astype(int)
    genes['ortholog_class'] = [
        'no_consensus' if c < 0 else ('1:1' if (nh, nm) == (1, 1) else '1:many' if nh == 1 else
                                      'many:1' if nm == 1 else 'many:many')
        for c, nh, nm in zip(genes.component, genes.n_human, genes.n_mouse)]
    return genes.drop(columns='root').reset_index(drop=True)


def classify_panel(panel, species, groups):
    """Attach orthogroup columns to a panel (index = symbol, column ``gene_id`` = Ensembl id).

    Matches by Ensembl id first and by symbol second. Panel genes absent from HCOP get class
    ``no_ortholog`` and component -1.
    """
    prefix = 'H' if species == 'human' else 'M'
    sub = groups[groups.species.eq(species)]
    by_id = sub.set_index('key')
    by_symbol = sub.drop_duplicates('symbol').set_index('symbol')
    out = []
    for symbol, gene_id in zip(panel.index.astype(str), panel.gene_id.astype(str)):
        key = f'{prefix}:{gene_id}'
        if key in by_id.index:
            row, how = by_id.loc[key], 'ensembl'
        elif symbol in by_symbol.index:
            row, how = by_symbol.loc[symbol], 'symbol'
        else:
            out.append({'gene': symbol, 'gene_id': gene_id, 'hcop_key': None, 'matched_by': 'none',
                        'max_support': 0, 'component': -1, 'n_human': 0, 'n_mouse': 0,
                        'ortholog_class': 'no_ortholog'})
            continue
        out.append({'gene': symbol, 'gene_id': gene_id, 'hcop_key': row.name if how == 'ensembl' else row.key,
                    'matched_by': how, 'max_support': int(row.max_support), 'component': int(row.component),
                    'n_human': int(row.n_human), 'n_mouse': int(row.n_mouse),
                    'ortholog_class': row.ortholog_class})
    return pd.DataFrame(out).set_index('gene')


# ---------------------------------------------------------------------------------------------
# Probe specificity
# ---------------------------------------------------------------------------------------------

_HEADER = re.compile(r'^>(\S+)\s+\S+\s+(?:chromosome|scaffold|primary_assembly):[^:]+:([^:]+):(\d+):(\d+):(-?1)'
                     r'.*?\sgene:(\S+)\s+gene_biotype:(\S+)\s+transcript_biotype:(\S+)(?:\s+gene_symbol:(\S+))?')


def transcript_table(fasta_paths):
    """Transcript id -> gene id (unversioned), symbol, biotypes and genomic interval, from Ensembl headers."""
    rows = []
    for path in fasta_paths:
        opener = gzip.open if str(path).endswith('.gz') else open
        with opener(path, 'rt') as handle:
            for line in handle:
                if line[0] != '>':
                    continue
                m = _HEADER.match(line)
                if m is None:
                    raise ValueError(f'Unparsed FASTA header in {path}: {line[:120]}')
                tx, chrom, start, end, strand, gene, gene_bt, tx_bt, symbol = m.groups()
                rows.append((tx.split('.')[0], gene.split('.')[0], symbol or '', gene_bt, tx_bt, chrom,
                             int(start), int(end), int(strand)))
    table = pd.DataFrame(rows, columns=['transcript', 'gene_id', 'symbol', 'gene_biotype', 'transcript_biotype',
                                        'chrom', 'start', 'end', 'strand'])
    return table.drop_duplicates('transcript').set_index('transcript')


def gene_intervals(transcripts):
    """Gene id -> chrom, strand, start (min), end (max), symbol and biotype."""
    return transcripts.groupby('gene_id').agg(chrom=('chrom', 'first'), strand=('strand', 'first'),
                                              start=('start', 'min'), end=('end', 'max'),
                                              symbol=('symbol', 'first'), gene_biotype=('gene_biotype', 'first'))


def write_combined_fasta(fasta_paths, out_path):
    """Concatenate (gzipped) FASTA files into one plain FASTA with unversioned transcript ids."""
    with open(out_path, 'w') as out:
        for path in fasta_paths:
            opener = gzip.open if str(path).endswith('.gz') else open
            with opener(path, 'rt') as handle:
                for line in handle:
                    if line[0] == '>':
                        out.write('>' + line[1:].split()[0].split('.')[0] + '\n')
                    else:
                        out.write(line)


def run_bowtie2(index_prefix, probes_fasta, sam_path, *, threads=16, max_mismatches=5):
    """All end-to-end, gapless alignments with <= ``max_mismatches`` mismatches (both strands).

    Exact 10-nt seeds at every offset: any 50-nt alignment with <= 4 mismatches contains an exact
    10-mer, so it is found (within bowtie2's extension effort, -D 30 -R 3); 5-mismatch hits are found
    when a 10-nt gap between mismatches exists. A 1-mismatch seed (-N 1) found the same alignments on
    a 560-probe benchmark at 66x the run time.
    """
    command = ['bowtie2', '-p', str(threads), '-f', '-x', str(index_prefix), '-U', str(probes_fasta),
               '--end-to-end', '-a', '--no-unal', '--no-hd', '--score-min', f'L,-{6 * max_mismatches},0',
               '--mp', '6,6', '--np', '6', '--rdg', '100,100', '--rfg', '100,100', '-N', '0', '-L', '10',
               '-i', 'C,1,0', '-D', '30', '-R', '3', '--reorder', '-S', str(sam_path)]
    subprocess.run(command, check=True, capture_output=True)


def read_alignments(sam_path):
    """Probe, transcript, reverse-strand flag and mismatches (XM tag) for each SAM record."""
    rows = []
    with open(sam_path) as handle:
        for line in handle:
            if line[0] == '@':
                continue
            parts = line.rstrip('\n').split('\t')
            flag = int(parts[1])
            if flag & 4:
                continue
            xm = next((int(p[5:]) for p in parts[11:] if p.startswith('XM:i:')), None)
            if xm is None or 'I' in parts[5] or 'D' in parts[5]:
                continue
            rows.append((parts[0], parts[2], bool(flag & 16), xm))
    return pd.DataFrame(rows, columns=['probe', 'transcript', 'reverse', 'mismatches'])


def probe_offtarget(alignments, probe_gene, transcripts, *, cross=3, near=5, count_pseudogenes=False):
    """Per-probe off-target audit.

    ``probe_gene`` maps probe id -> own Ensembl gene id. The on-target strand is the majority strand
    of own-gene alignments over all probes. Off-target = an alignment on that strand to another gene,
    excluding genes whose interval overlaps the own gene on the same chromosome strand (read-through
    or nested loci), genes with the own gene's symbol (duplicate annotations of one gene, e.g. on an
    unplaced scaffold) and, unless ``count_pseudogenes``, pseudogene biotypes. Returns one row per probe:
    on_target (bool), best_offtarget_mm (NaN when none), offtarget_genes (';'-joined symbols at the
    best mismatch level) and status (``cross`` <= cross, ``near`` <= near, ``unique``, ``unverified``).
    """
    a = alignments.copy()
    a['gene_id'] = a.transcript.map(transcripts.gene_id)
    if a.gene_id.isna().any():
        raise ValueError('Alignments to transcripts missing from the transcript table.')
    a['own'] = a.probe.map(probe_gene)
    own_hits = a[a.gene_id.eq(a.own)]
    if own_hits.empty:
        raise ValueError('No probe aligned to its own gene; check the probe and transcript inputs.')
    on_reverse = bool(own_hits.reverse.mean() > .5)
    a = a[a.reverse.eq(on_reverse)]
    intervals = gene_intervals(transcripts)
    off = a[a.gene_id.ne(a.own)].copy()
    if len(off):
        own_iv = intervals.reindex(off.own.to_numpy())
        off_iv = intervals.reindex(off.gene_id.to_numpy())
        same_locus = ((own_iv.chrom.to_numpy() == off_iv.chrom.to_numpy())
                      & (own_iv.strand.to_numpy() == off_iv.strand.to_numpy())
                      & (own_iv.start.to_numpy() <= off_iv.end.to_numpy())
                      & (off_iv.start.to_numpy() <= own_iv.end.to_numpy()))
        own_symbol = own_iv.symbol.fillna('').to_numpy()
        same_symbol = (own_symbol != '') & (own_symbol == off_iv.symbol.fillna('').to_numpy())
        off = off[~(same_locus | same_symbol)]
        if not count_pseudogenes:
            off = off[~off.gene_id.map(intervals.gene_biotype).astype(str).str.contains(PSEUDOGENE)]
        off['symbol'] = off.gene_id.map(intervals.symbol)
    best = off.groupby('probe').mismatches.min() if len(off) else pd.Series(dtype=float)
    genes_at_best = (off.merge(best.rename('best'), left_on='probe', right_index=True)
                     .query('mismatches == best').groupby('probe').symbol
                     .agg(lambda s: ';'.join(sorted(set(s))))) if len(off) else pd.Series(dtype=str)
    probes = pd.Index(sorted(probe_gene.index if hasattr(probe_gene, 'index') else probe_gene), name='probe')
    out = pd.DataFrame(index=probes)
    out['on_target'] = out.index.isin(set(own_hits[own_hits.reverse.eq(on_reverse)].probe))
    out['best_offtarget_mm'] = best.reindex(out.index)
    out['offtarget_genes'] = genes_at_best.reindex(out.index).fillna('')
    mm = out.best_offtarget_mm
    out['status'] = np.where(~out.on_target, 'unverified',
                             np.where(mm.le(cross), 'cross', np.where(mm.le(near), 'near', 'unique')))
    out.attrs['on_target_reverse'] = on_reverse
    return out


def gene_probe_status(probe_status, probe_gene):
    """Gene-level status: ambiguous > unverified > near > clean (worst probe wins); n probes."""
    table = probe_status.assign(gene_id=pd.Series(probe_gene).reindex(probe_status.index).to_numpy())
    rank = {'unique': 0, 'near': 1, 'unverified': 2, 'cross': 3}
    names = {0: 'clean', 1: 'near', 2: 'unverified', 3: 'ambiguous'}
    worst = table.status.map(rank).groupby(table.gene_id).max().map(names)
    return pd.DataFrame({'probe_status': worst, 'n_probes': table.groupby('gene_id').size(),
                         'best_offtarget_mm': table.groupby('gene_id').best_offtarget_mm.min(),
                         'offtarget_genes': table.groupby('gene_id').offtarget_genes.agg(
                             lambda s: ';'.join(sorted({g for x in s for g in x.split(';') if g})))})


# ---------------------------------------------------------------------------------------------
# Curves and tests
# ---------------------------------------------------------------------------------------------

def weighted_nested_f(expression, weights, full, null):
    """Weighted least-squares nested F-test for every column of ``expression`` at once.

    ``expression`` is structures x genes (sparse or dense), ``weights`` positive per structure,
    ``full`` and ``null`` design matrices with the null nested in the full model. Returns
    (coef_full (p x genes), F, p-value, df1, df2).
    """
    w = np.asarray(weights, float)
    if (w <= 0).any() or not np.isfinite(w).all():
        raise ValueError('Weights must be positive and finite.')
    X1, X0 = np.asarray(full, float), np.asarray(null, float)
    Y = sparse.csr_matrix(expression) if sparse.issparse(expression) else np.asarray(expression, float)
    if sparse.issparse(Y):
        wy = sparse.diags(w) @ Y
        yy = np.asarray(Y.multiply(Y).T @ w).ravel()
        xty1 = np.asarray(wy.T @ X1).T
        xty0 = np.asarray(wy.T @ X0).T
    else:
        wy = Y * w[:, None]
        yy = (Y * Y).T @ w
        xty1, xty0 = X1.T @ wy, X0.T @ wy
    g1 = X1.T @ (X1 * w[:, None])
    g0 = X0.T @ (X0 * w[:, None])
    coef1 = np.linalg.solve(g1, xty1)
    coef0 = np.linalg.solve(g0, xty0)
    rss1 = np.maximum(yy - np.einsum('pg,pg->g', coef1, xty1), 0)
    rss0 = np.maximum(yy - np.einsum('pg,pg->g', coef0, xty0), 0)
    df1 = np.linalg.matrix_rank(X1) - np.linalg.matrix_rank(X0)
    df2 = len(w) - np.linalg.matrix_rank(X1)
    with np.errstate(divide='ignore', invalid='ignore'):
        F = ((rss0 - rss1) / df1) / (rss1 / df2)
    F = np.where(np.isfinite(F), F, 0.)
    return coef1, F, f_dist.sf(F, df1, df2), df1, df2


def log_profile(curve, pseudocount):
    """log2(max(curve, 0) + pseudocount)."""
    return np.log2(np.maximum(np.asarray(curve, float), 0) + pseudocount)


def rowwise_pearson(a, b):
    """Pearson correlation of matching rows of two (genes x grid) arrays (NaN for flat rows)."""
    a = np.asarray(a, float) - np.mean(a, axis=-1, keepdims=True)
    b = np.asarray(b, float) - np.mean(b, axis=-1, keepdims=True)
    den = np.sqrt((a * a).sum(-1) * (b * b).sum(-1))
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(den > 0, (a * b).sum(-1) / den, np.nan)


def rms_shape_distance(a, b):
    """RMS difference of mean-centred profiles (last axis = grid)."""
    a = np.asarray(a, float) - np.mean(a, axis=-1, keepdims=True)
    b = np.asarray(b, float) - np.mean(b, axis=-1, keepdims=True)
    return np.sqrt(np.mean((a - b) ** 2, axis=-1))


def separation_call(profiles, groups, *, min_gap=.5):
    """Two-group shape separation for four (or more) specimen profiles of one unit.

    ``profiles`` is specimens x grid, ``groups`` two labels per specimen (two specimens each). Call
    when every within-group distance is smaller than every between-group distance and the mean
    between-group distance exceeds the mean within-group distance by ``min_gap``. Returns
    (call, gap, mean within, mean between).
    """
    P = np.asarray(profiles, float)
    g = np.asarray(groups)
    n = len(P)
    within, between = [], []
    for i in range(n):
        for j in range(i + 1, n):
            d = float(rms_shape_distance(P[i], P[j]))
            (within if g[i] == g[j] else between).append(d)
    if not within or not between:
        raise ValueError('Need at least two groups with replicated members.')
    gap = float(np.mean(between) - np.mean(within))
    return bool(max(within) < min(between) and gap >= min_gap), gap, float(np.mean(within)), float(np.mean(between))


# ---------------------------------------------------------------------------------------------
# Paralog partitioning and the union fit
# ---------------------------------------------------------------------------------------------

def pair_statistics(peaks, standardized, i, j):
    """Peak separation |peak_i - peak_j| and dissimilarity 1 - r for index arrays i, j.

    ``standardized`` rows are mean-centred, unit-norm profiles, so r is a row dot product.
    """
    i, j = np.asarray(i), np.asarray(j)
    sep = np.abs(np.asarray(peaks)[i] - np.asarray(peaks)[j])
    r = np.einsum('...k,...k->...', standardized[i], standardized[j])
    return sep, 1 - r


def standardize_rows(profiles):
    """Mean-centre and scale each row to unit Euclidean norm (flat rows stay zero)."""
    x = np.asarray(profiles, float) - np.mean(profiles, axis=1, keepdims=True)
    norm = np.linalg.norm(x, axis=1, keepdims=True)
    return np.divide(x, norm, out=np.zeros_like(x), where=norm > 0)


def matched_pair_draws(pairs, strata, family, eligible, n_draws, rng):
    """Random partner pairs matched on stratum, outside the observed pair's family.

    ``pairs`` is an (n, 2) integer array of gene indices; ``strata`` an integer label per gene;
    ``family`` a family id per gene (-1 = none); ``eligible`` a boolean mask. Returns two (n, n_draws)
    index arrays (a, b) with a != b, stratum(a) = stratum(i), stratum(b) = stratum(j), and neither
    a nor b in the family of the observed pair. Raises when a stratum has no candidate.
    """
    pairs = np.asarray(pairs, int)
    strata, family, eligible = np.asarray(strata), np.asarray(family), np.asarray(eligible, bool)
    out_a = np.empty((len(pairs), n_draws), int)
    out_b = np.empty((len(pairs), n_draws), int)
    for k, (i, j) in enumerate(pairs):
        fam = family[i]
        outside = eligible & ((family != fam) | (fam < 0))
        cand_a = np.flatnonzero(outside & (strata == strata[i]))
        cand_b = np.flatnonzero(outside & (strata == strata[j]))
        if len(cand_a) == 0 or len(cand_b) == 0 or (len(cand_a) == 1 and np.array_equal(cand_a, cand_b)):
            raise ValueError(f'No matched candidates for pair {(i, j)}.')
        a = rng.choice(cand_a, n_draws)
        b = rng.choice(cand_b, n_draws)
        clash = a == b
        while clash.any():
            b[clash] = rng.choice(cand_b, int(clash.sum()))
            clash = a == b
        out_a[k], out_b[k] = a, b
    return out_a, out_b


def global_two_sided(observed, null):
    """Two-sided Monte-Carlo p-value of ``observed`` against draws ``null`` (with +1 correction)."""
    null = np.asarray(null, float)
    lo = (1 + np.sum(null <= observed)) / (1 + len(null))
    hi = (1 + np.sum(null >= observed)) / (1 + len(null))
    return float(min(1., 2 * min(lo, hi)))


def unit_mean(profiles):
    """Scale each row (linear-scale curve) to mean 1 after flooring at 0."""
    x = np.maximum(np.asarray(profiles, float), 0)
    m = x.mean(axis=-1, keepdims=True)
    return np.divide(x, m, out=np.zeros_like(x), where=m > 0)


def r_squared(y, fit):
    y, fit = np.asarray(y, float), np.asarray(fit, float)
    tss = np.sum((y - y.mean()) ** 2)
    return float(1 - np.sum((y - fit) ** 2) / tss) if tss > 0 else np.nan


def nnls_union(target, components):
    """Non-negative fit of ``target`` (grid) by rows of ``components`` (k x grid).

    Returns (R^2 of the union fit, coefficients, R^2 of the best single non-negative scaled
    component, index of that component). R^2 is around the target mean.
    """
    y = np.asarray(target, float)
    X = np.asarray(components, float)
    coef, _ = nnls(X.T, y)
    union = r_squared(y, X.T @ coef)
    singles = []
    for row in X:
        c, _ = nnls(row[:, None], y)
        singles.append(r_squared(y, row * c[0]))
    best = int(np.nanargmax(singles))
    return union, coef, float(singles[best]), best


def breadth(profile, level=.5):
    """Fraction of grid points where a linear-scale profile is >= ``level`` x its maximum."""
    p = np.maximum(np.asarray(profile, float), 0)
    return float(np.mean(p >= level * p.max())) if p.max() > 0 else np.nan


def symbol_root(symbol):
    """Gene-family root, upper case: the symbol without its last member suffix.

    Slc22a6 -> SLC22A, Cyp4a12a -> CYP4A, Ces1d -> CES1, Gsta1 -> GSTA, Kap -> KAP. A heuristic
    for grouping family members by name only; it is not a homology call.
    """
    s = str(symbol).upper()
    for pattern in (r'^(.*?\d+[A-Z])\d+[A-Z]*$', r'^(.*\d)[A-Z]+$', r'^([A-Z]+)\d'):
        m = re.match(pattern, s)
        if m:
            return m.group(1)
    return s


# ---------------------------------------------------------------------------------------------
# Direct probe-to-paralog distance (no seed heuristics)
# ---------------------------------------------------------------------------------------------

_COMPLEMENT = str.maketrans('ACGTN', 'TGCAN')


def reverse_complement(seq):
    return str(seq).upper().translate(_COMPLEMENT)[::-1]


def transcript_sequences(fasta_paths, gene_ids):
    """Gene id (unversioned) -> list of transcript sequences, for the requested genes only."""
    wanted = set(gene_ids)
    out, current, chunks = {}, None, []
    for path in fasta_paths:
        opener = gzip.open if str(path).endswith('.gz') else open
        with opener(path, 'rt') as handle:
            for line in handle:
                if line[0] == '>':
                    if current is not None:
                        out.setdefault(current, []).append(''.join(chunks))
                    m = re.search(r'\sgene:(\S+)', line)
                    gene = m.group(1).split('.')[0] if m else None
                    current, chunks = (gene if gene in wanted else None), []
                elif current is not None:
                    chunks.append(line.strip().upper())
        if current is not None:
            out.setdefault(current, []).append(''.join(chunks))
            current, chunks = None, []
    return out


def min_window_mismatches(probe, sequences, *, reverse=True):
    """Smallest Hamming distance between a probe and any same-length window of ``sequences``.

    With ``reverse`` (the panel's probes are antisense to the transcript), the probe's reverse
    complement is compared with the sense transcript. Gapped matches are not considered. Returns
    NaN when no sequence is at least as long as the probe.
    """
    query = reverse_complement(probe) if reverse else str(probe).upper()
    q = np.frombuffer(query.encode(), np.uint8)
    best = np.inf
    for seq in sequences:
        s = np.frombuffer(seq.encode(), np.uint8)
        if len(s) < len(q):
            continue
        windows = np.lib.stride_tricks.sliding_window_view(s, len(q))
        best = min(best, int((windows != q).sum(axis=1).min()))
    return float(best) if np.isfinite(best) else np.nan
