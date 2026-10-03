import numpy as np
import pytest
pytest.importorskip("scipy")
pytest.importorskip("sklearn")

from pseudospace.atlas_initialization import initialize_atlas_positions
from pseudospace.repeated_atlas_baselines import fit_pc1_reference, project_pc1_reference


def _data(seed=4):
    rng = np.random.default_rng(seed)
    anatomy = np.repeat([0, 1, 2], 30)
    z = np.repeat([0.12, 0.5, 0.88], 30) + rng.normal(0, 0.035, 90)
    specimens = np.tile(np.repeat(['A', 'B'], 15), 3)
    Y = np.column_stack([z, z**2, rng.normal(size=(90, 3))])
    Y += np.where(specimens[:, None] == 'B', np.array([0.15, -0.1, 0.05, 0.2, -0.1]), 0)
    return Y, anatomy, specimens


def test_pc1_matches_existing_frozen_wrapper_logic():
    Y, anatomy, specimens = _data()
    z, meta = initialize_atlas_positions(Y, anatomy, specimens, method='pc1')
    ref = fit_pc1_reference(Y, anatomy)
    np.testing.assert_allclose(z, project_pc1_reference(ref, Y))
    assert meta['low'] == ref['low'] and meta['high'] == ref['high']


def test_endpoint_contrast_equal_specimen_weighting():
    Y, anatomy, specimens = _data()
    _, m1 = initialize_atlas_positions(Y, anatomy, specimens, method='endpoint_contrast')
    # Duplicate all rows from specimen B; fixed observations should not change
    # the endpoint direction when specimens, rather than rows, are averaged.
    keep = specimens == 'B'
    Y2 = np.vstack([Y, Y[keep]])
    a2 = np.r_[anatomy, anatomy[keep]]
    s2 = np.r_[specimens, np.repeat('B', keep.sum())]
    _, m2 = initialize_atlas_positions(Y2, a2, s2, method='endpoint_contrast')
    np.testing.assert_allclose(m1['direction'], m2['direction'], atol=1e-12)


def test_soft_segment_seed_reproducible():
    Y, anatomy, specimens = _data()
    z1, _ = initialize_atlas_positions(Y, anatomy, specimens, method='soft_segment', seed=29)
    z2, _ = initialize_atlas_positions(Y, anatomy, specimens, method='soft_segment', seed=29)
    np.testing.assert_allclose(z1, z2)
    assert np.all((z1 >= 0) & (z1 <= 1))


@pytest.mark.parametrize('method', ['pc1', 'endpoint_contrast'])
def test_endpoint_methods_reject_missing_endpoint_labels(method):
    Y, anatomy, specimens = _data()
    with pytest.raises(ValueError, match='endpoint labels'):
        initialize_atlas_positions(Y, np.ones_like(anatomy), specimens, method=method)


def test_soft_segment_requires_all_three_labels():
    Y, anatomy, specimens = _data()
    keep = anatomy != 1
    with pytest.raises(ValueError, match='all labels'):
        initialize_atlas_positions(Y[keep], anatomy[keep], specimens[keep], method='soft_segment')


def test_endpoint_contrast_rejects_degenerate_direction_and_bad_method():
    rng = np.random.default_rng(8)
    anatomy = np.repeat([0, 2], 12)
    block = rng.normal(size=(12, 3))
    Y = np.vstack([block, block.copy()])
    with pytest.raises(ValueError, match='degenerate'):
        initialize_atlas_positions(Y, anatomy, np.repeat(['A', 'B'], 12), method='endpoint_contrast')
    with pytest.raises(ValueError, match='method'):
        initialize_atlas_positions(Y, anatomy, np.repeat(['A', 'B'], 12), method='dpt')
