import numpy as np
import pytest
from scipy import sparse

from pseudospace.replication_rules import arm_passes, csr_group_sums, gene_verdict, story_verdict


def test_csr_group_sums_matches_dense_and_handles_blocks(tmp_path):
    h5py = pytest.importorskip('h5py')
    rng = np.random.default_rng(0)
    dense = rng.poisson(.4, size=(23, 7)).astype(float)
    matrix = sparse.csr_matrix(dense)
    path = tmp_path / 'm.h5'
    with h5py.File(path, 'w') as h:
        h['data'], h['indices'], h['indptr'] = matrix.data, matrix.indices, matrix.indptr
    rows = np.array([20, 0, 3, 4, 5, 11, 22])
    codes = np.array([1, 0, 0, 1, 1, 2, 2])
    with h5py.File(path, 'r') as h:
        sums, totals, counts = csr_group_sums(h['data'], h['indices'], h['indptr'], rows, codes, 3, 7, block_rows=4)
    for g in range(3):
        expected = dense[rows[codes == g]].sum(axis=0)
        np.testing.assert_allclose(sums[g], expected)
        assert totals[g] == pytest.approx(expected.sum())
    np.testing.assert_array_equal(counts, [2, 3, 2])


def test_csr_group_sums_rejects_duplicates():
    matrix = sparse.csr_matrix(np.eye(3))
    with pytest.raises(ValueError):
        csr_group_sums(matrix.data, matrix.indices, matrix.indptr, [0, 0], [0, 0], 1, 3)


def test_arm_passes_signs_dead_zone_and_flat():
    assert arm_passes('+', .5) is True
    assert arm_passes('+', .05) is False          # inside the dead zone is not a sign
    assert arm_passes('-', -.3) is True
    assert arm_passes('', 1.0) is None
    assert arm_passes('+', float('nan')) is None
    assert arm_passes('0', .2, other=-1.0) is True
    assert arm_passes('0', .8, other=-1.0) is False
    assert arm_passes('0', .01, other=.05) is False  # both flat: no species difference
    assert arm_passes('0', .01, other=None) is None


def test_gene_verdict_rules():
    critical = {'ours', 'human_lake'}
    ok = dict(ours=True, human_lake=True, mouse_sn_male=True, mouse_sn_female=True)
    assert gene_verdict(ok, critical=critical, probe_balanced=True, detection_ok=True) == 'headline'
    assert gene_verdict(ok, critical=critical, probe_balanced=False, detection_ok=True) == 'supporting'
    one = dict(ok, mouse_sn_female=False)
    assert gene_verdict(one, critical=critical, probe_balanced=True, detection_ok=True) == 'supporting'
    two = dict(one, mouse_sn_male=False)
    assert gene_verdict(two, critical=critical, probe_balanced=True, detection_ok=True) == 'not replicated'
    crit = dict(ok, human_lake=False)
    assert gene_verdict(crit, critical=critical, probe_balanced=True, detection_ok=True) == 'not replicated'
    untested = dict(ok, human_lake=None)
    assert gene_verdict(untested, critical=critical, probe_balanced=True, detection_ok=True) == 'supporting'
    assert gene_verdict({'a': None}, critical=critical, probe_balanced=True, detection_ok=True) == 'not testable'


def test_story_verdict():
    verdicts = {'A': 'headline', 'B': 'supporting', 'C': 'not replicated', 'D': 'headline'}
    balanced = {'A': True, 'B': True, 'C': True, 'D': False}
    assert story_verdict(verdicts, ['A'], balanced)[0] == 'lead'                 # 2 of 3 balanced pass
    assert story_verdict(verdicts, ['B'], balanced)[0] == 'supporting'           # lead not headline
    assert story_verdict(verdicts, ['A'], balanced, can_lead=False)[0] == 'supporting'
    verdicts['B'] = 'not replicated'
    assert story_verdict(verdicts, ['A'], balanced)[0] == 'not replicated'
