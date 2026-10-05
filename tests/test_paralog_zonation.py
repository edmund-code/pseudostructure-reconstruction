import gzip

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from pseudospace import paralog_zonation as pz


def _hcop():
    rows = [
        # 1:1 with strong support
        ('ENSG1', 'A', 'ENSMUSG1', 'A', 'HGNC,Ensembl,OMA,Panther'),
        # one human gene, two mouse co-orthologs
        ('ENSG2', 'B', 'ENSMUSG2', 'B1', 'HGNC,Ensembl,OMA'),
        ('ENSG2', 'B', 'ENSMUSG3', 'B2', 'Panther,OrthoDB,OMA'),
        # a low-support link that must not join the family at support >= 3
        ('ENSG2', 'B', 'ENSMUSG9', 'Z', 'OrthoMCL'),
        # many:many
        ('ENSG4', 'C1', 'ENSMUSG4', 'C', 'HGNC,Ensembl,OMA'),
        ('ENSG5', 'C2', 'ENSMUSG4', 'C', 'HGNC,Ensembl,OMA'),
        ('ENSG5', 'C2', 'ENSMUSG5', 'Cx', 'HGNC,Ensembl,OMA'),
        # symbol-only key
        ('-', 'D', 'ENSMUSG6', 'D', 'HGNC,Ensembl,OMA'),
    ]
    return pd.DataFrame(rows, columns=['human_ensembl_gene', 'human_symbol', 'mouse_ensembl_gene', 'mouse_symbol',
                                       'support'])


def test_support_count_and_orthogroup_classes():
    assert pz.support_count('A,B, B') == 2 and pz.support_count('-') == 0
    groups = pz.hcop_orthogroups(_hcop(), min_support=3).set_index('key')
    assert groups.loc['H:ENSG1', 'ortholog_class'] == '1:1'
    assert groups.loc['H:ENSG2', 'ortholog_class'] == '1:many'
    assert groups.loc['M:ENSMUSG2', 'component'] == groups.loc['M:ENSMUSG3', 'component']
    assert groups.loc['M:ENSMUSG9', 'ortholog_class'] == 'no_consensus'
    assert groups.loc['M:ENSMUSG4', 'ortholog_class'] == 'many:many'
    assert groups.loc['H:sym:D', 'ortholog_class'] == '1:1'
    # at support >= 1 the low-support link joins the family
    loose = pz.hcop_orthogroups(_hcop(), min_support=1).set_index('key')
    assert loose.loc['H:ENSG2', 'n_mouse'] == 3


def test_classify_panel_matches_by_id_then_symbol():
    groups = pz.hcop_orthogroups(_hcop(), min_support=3)
    panel = pd.DataFrame({'gene_id': ['ENSMUSG2', 'ENSMUSGXX', 'ENSMUSG77']}, index=['B1', 'D', 'Kap'])
    out = pz.classify_panel(panel, 'mouse', groups)
    assert out.loc['B1', 'ortholog_class'] == '1:many' and out.loc['B1', 'matched_by'] == 'ensembl'
    assert out.loc['D', 'matched_by'] == 'symbol'
    assert out.loc['Kap', 'ortholog_class'] == 'no_ortholog' and out.loc['Kap', 'component'] == -1


def test_weighted_nested_f_matches_lstsq_dense_and_sparse():
    rng = np.random.default_rng(0)
    n = 200
    x = rng.uniform(size=n)
    w = rng.uniform(.5, 2, size=n)
    full = np.column_stack([np.ones(n), x, x ** 2])
    null = full[:, :1]
    Y = np.column_stack([3 + 2 * x + rng.normal(size=n), rng.poisson(1, size=n).astype(float)])
    coef, F, p, df1, df2 = pz.weighted_nested_f(Y, w, full, null)
    sw = np.sqrt(w)
    for g in range(2):
        b1, *_ = np.linalg.lstsq(full * sw[:, None], Y[:, g] * sw, rcond=None)
        r1 = np.sum(w * (Y[:, g] - full @ b1) ** 2)
        r0 = np.sum(w * (Y[:, g] - np.average(Y[:, g], weights=w)) ** 2)
        assert np.allclose(coef[:, g], b1)
        assert np.isclose(F[g], ((r0 - r1) / 2) / (r1 / (n - 3)))
    assert (df1, df2) == (2, n - 3)
    coef_s, F_s, p_s, *_ = pz.weighted_nested_f(sparse.csr_matrix(Y), w, full, null)
    assert np.allclose(coef_s, coef) and np.allclose(F_s, F) and p[0] < 1e-10


def test_separation_call_and_distances():
    grid = np.linspace(0, 1, 21)
    up, down = grid, 1 - grid
    profiles = np.stack([up, up + .01, down, down + .02])
    call, gap, within, between = pz.separation_call(profiles, ['m', 'm', 'h', 'h'], min_gap=.3)
    assert call and gap > .3 and within < between
    call2, *_ = pz.separation_call(profiles, ['a', 'b', 'a', 'b'], min_gap=.3)
    assert not call2
    assert np.isclose(pz.rms_shape_distance(up, up + 5), 0)


def test_matched_pair_draws_respect_strata_and_family():
    rng = np.random.default_rng(1)
    strata = np.array([0, 0, 1, 1, 0, 1, 0, 1])
    family = np.array([7, 7, 7, -1, -1, -1, 3, 3])
    eligible = np.ones(8, bool)
    a, b = pz.matched_pair_draws(np.array([[0, 2]]), strata, family, eligible, 500, rng)
    assert set(a.ravel()) <= {4, 6} and set(b.ravel()) <= {3, 5, 7}
    assert not np.any(a == b)
    assert pz.global_two_sided(10., np.arange(100.)) == pytest.approx(24 / 101)
    assert pz.global_two_sided(50., np.arange(100.)) == 1.


def test_nnls_union_and_breadth():
    grid = np.linspace(0, 1, 41)
    left = np.exp(-((grid - .2) / .1) ** 2)
    right = np.exp(-((grid - .8) / .1) ** 2)
    union, coef, best_single, best = pz.nnls_union(left + 2 * right, np.stack([left, right]))
    assert union == pytest.approx(1.0) and np.allclose(coef, [1, 2], atol=1e-6)
    assert best_single < .9 and best == 1
    assert pz.breadth(left + right) > pz.breadth(left)
    assert np.allclose(pz.unit_mean(np.stack([left, 3 * left])).mean(axis=1), 1)


def test_symbol_root():
    assert pz.symbol_root('Slc22a6') == 'SLC22A'
    assert pz.symbol_root('Cyp4a12a') == 'CYP4A'
    assert pz.symbol_root('Ces1d') == 'CES1'
    assert pz.symbol_root('Gsta1') == 'GSTA'
    assert pz.symbol_root('Kap') == 'KAP'


def _fasta(tmp_path):
    path = tmp_path / 'tx.fa.gz'
    headers = [
        '>T1.1 cdna chromosome:GRCm38:1:100:200:1 gene:G1.1 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:Own',
        '>T2.1 cdna chromosome:GRCm38:1:150:400:1 gene:G2.1 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:Readthrough',
        '>T3.1 cdna chromosome:GRCm38:5:100:200:1 gene:G3.1 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:Paralog',
        '>T4.1 cdna chromosome:GRCm38:7:100:200:-1 gene:G4.1 gene_biotype:processed_pseudogene transcript_biotype:processed_pseudogene gene_symbol:Pseudo',
        '>T5.1 ncrna scaffold:GRCm38:JH1:1:50:1 gene:G5.1 gene_biotype:lncRNA transcript_biotype:lncRNA',
        '>T6.1 cdna scaffold:GRCm38:JH2:1:50:1 gene:G6.1 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:Own',
    ]
    with gzip.open(path, 'wt') as handle:
        for h in headers:
            handle.write(h + '\nACGT\n')
    return path


def test_transcript_table_and_probe_offtarget(tmp_path):
    tx = pz.transcript_table([_fasta(tmp_path)])
    assert tx.loc['T1', 'gene_id'] == 'G1' and tx.loc['T5', 'symbol'] == '' and tx.loc['T4', 'strand'] == -1
    probe_gene = pd.Series({'p1': 'G1', 'p2': 'G1', 'p3': 'G1', 'p4': 'G3'})
    alignments = pd.DataFrame([
        ('p1', 'T1', False, 0), ('p1', 'T2', False, 0),    # read-through locus only -> unique
        ('p1', 'T6', False, 0),                            # same symbol on a scaffold -> ignored
        ('p2', 'T1', False, 0), ('p2', 'T3', False, 2),    # paralog at 2 mismatches -> cross
        ('p2', 'T3', True, 0),                             # wrong strand: ignored
        ('p3', 'T1', False, 0), ('p3', 'T4', False, 1), ('p3', 'T3', False, 5),  # pseudogene ignored -> near
        ('p4', 'T1', False, 4),                            # never hits its own gene -> unverified
    ], columns=['probe', 'transcript', 'reverse', 'mismatches'])
    status = pz.probe_offtarget(alignments, probe_gene, tx)
    assert status.loc['p1', 'status'] == 'unique'
    assert status.loc['p2', 'status'] == 'cross' and status.loc['p2', 'offtarget_genes'] == 'Paralog'
    assert status.loc['p3', 'status'] == 'near'
    assert status.loc['p4', 'status'] == 'unverified'
    with_pseudo = pz.probe_offtarget(alignments, probe_gene, tx, count_pseudogenes=True)
    assert with_pseudo.loc['p3', 'status'] == 'cross'
    genes = pz.gene_probe_status(status, probe_gene)
    assert genes.loc['G1', 'probe_status'] == 'ambiguous' and genes.loc['G1', 'n_probes'] == 3
    assert genes.loc['G3', 'probe_status'] == 'unverified'


def test_min_window_mismatches_and_sequences(tmp_path):
    path = tmp_path / 'tx.fa.gz'
    with gzip.open(path, 'wt') as handle:
        handle.write('>T1.1 cdna x gene:G1.1 y\nAACCGGTTAC\nGTACGT\n>T2.1 cdna x gene:G2.1 y\nTTTT\n')
    seqs = pz.transcript_sequences([path], ['G1'])
    assert seqs == {'G1': ['AACCGGTTACGTACGT']}
    probe_sense = 'GGTTACGT'
    assert pz.min_window_mismatches(probe_sense, seqs['G1'], reverse=False) == 0
    assert pz.min_window_mismatches(pz.reverse_complement(probe_sense), seqs['G1'], reverse=True) == 0
    assert pz.min_window_mismatches('GGTAACGT', seqs['G1'], reverse=False) == 1
    assert np.isnan(pz.min_window_mismatches('A' * 40, seqs['G1'], reverse=False))
