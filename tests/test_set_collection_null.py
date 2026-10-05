import numpy as np
import pandas as pd

from pseudospace.pathway_decoys import joint_matched_test
from pseudospace.pathway_remodeling import correlation_adjusted_rank_tests, residual_pathway_correlations
from pseudospace.set_collection_null import (JointReplicationZ, collection_null, label_permutation_p,
                                             stratified_permutation)
from pseudospace.specimen_null import matched_auc_null


def _example(seed=0):
    rng = np.random.default_rng(seed)
    genes = [f'g{i}' for i in range(200)]
    residuals = rng.normal(size=(12, 200))
    residuals[:, :15] += rng.normal(size=(12, 1)) * 1.5   # a correlated block
    statistic = pd.Series(rng.normal(size=200), index=genes)
    statistic.iloc[:15] += 1.
    strata = rng.integers(0, 4, 200)
    sets = {'a': genes[:20], 'b': genes[10:40], 'c': genes[100:130]}
    return statistic, strata, sets, residuals, genes


def test_identity_matches_notebook_37_joint_test():
    statistic, strata, sets, residuals, genes = _example()
    scorer = JointReplicationZ(statistic, strata, sets, residuals)
    corr = residual_pathway_correlations(residuals, np.repeat('donors', len(residuals)), genes, sets)
    vif = correlation_adjusted_rank_tests(statistic, sets, corr).set_index('pathway_id').variance_inflation
    expected = joint_matched_test(matched_auc_null(statistic.to_frame('a'), sets, strata), vif, keys=())
    np.testing.assert_allclose(scorer.z(), expected.set_index('pathway_id').loc[list(sets), 'z_joint'], rtol=1e-8)


def test_permutation_equals_rescoring_relabelled_sets():
    statistic, strata, sets, residuals, genes = _example(1)
    scorer = JointReplicationZ(statistic, strata, sets, residuals)
    perm = stratified_permutation(strata, np.random.default_rng(3))
    assert (strata[perm] == strata).all()
    moved = {p: [genes[perm[genes.index(g)]] for g in v] for p, v in sets.items()}
    assert len(set(moved['a']) & set(moved['b'])) == len(set(sets['a']) & set(sets['b']))
    np.testing.assert_allclose(scorer.z(perm), JointReplicationZ(statistic, strata, moved, residuals).z(), rtol=1e-8)


def test_collection_null_and_label_permutation():
    statistic, strata, sets, residuals, _ = _example(2)
    scorer = JointReplicationZ(statistic, strata, sets, residuals)
    observed, null = collection_null(scorer, {'first': ['a', 'b'], 'other': ['c']}, strata, 50, seed=1)
    assert null.shape == (50, 2) and np.isfinite(null.to_numpy()).all()
    assert observed['first'] > null['first'].median()
    p, diff = label_permutation_p([3., 2.5, 2.], [0., .1, -.2])
    assert diff > 0 and p == 1 / 20
