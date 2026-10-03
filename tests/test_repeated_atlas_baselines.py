"""Synthetic checks for frozen PC1 and soft segment-mixture baselines."""
import numpy as np
import pytest

pytest.importorskip('sklearn')

from pseudospace.repeated_atlas_baselines import (
    fit_pc1_reference, fit_segment_mixture, predict_segment_mixture,
    project_pc1_reference,
)


def test_pc1_reference_is_training_oriented_and_query_changes_do_not_refit_it():
    rng = np.random.default_rng(2)
    z = np.linspace(0, 1, 90)
    anatomy = np.repeat([0, 1, 2], 30)
    Y = np.column_stack([z, 2*z, -z]) + rng.normal(0, .01, (90, 3))
    reference = fit_pc1_reference(Y, anatomy)
    coordinate = project_pc1_reference(reference, Y)
    assert np.median(coordinate[anatomy == 2]) > np.median(coordinate[anatomy == 0])
    extreme_query = Y[::10] + 1e5
    projected = project_pc1_reference(reference, extreme_query)
    assert np.all((projected >= 0) & (projected <= 1))
    changed_later_pcs = extreme_query.copy()
    changed_later_pcs[:, 1:] = -1e8
    np.testing.assert_array_equal(project_pc1_reference(reference, changed_later_pcs), projected)
    np.testing.assert_allclose(project_pc1_reference(reference, Y), coordinate)
    assert np.isfinite(reference['low']) and np.isfinite(reference['high'])


def test_pc1_orientation_needs_only_s1_s3_endpoints():
    train = np.array([[0., 8.], [.1, -4.], [.8, 10.], [1., -9.]])
    reference = fit_pc1_reference(train, np.array([0, 0, 2, 2]))
    assert reference['sign'] == 1.
    coord = project_pc1_reference(reference, train)
    assert np.median(coord[2:]) > np.median(coord[:2])


def _mixture_data():
    feature_centers = np.array([[-4., 0.], [0., 4.], [4., 0.]])
    class_response = np.array([[10., 0.], [0., 20.], [30., 20.]])
    sample_shift = {'A': np.array([2., -2.]), 'B': np.array([-2., 2.])}
    class_counts = {'A': [70, 20, 10], 'B': [10, 30, 60]}
    Ys, responses, labels, samples = [], [], [], []
    for sample in ('A', 'B'):
        for k, n in enumerate(class_counts[sample]):
            Ys.append(feature_centers[k] + np.zeros((n, 2)))
            responses.append(class_response[k] + sample_shift[sample] + np.zeros((n, 2)))
            labels.extend([k] * n)
            samples.extend([sample] * n)
    return (np.vstack(Ys), np.vstack(responses), np.asarray(labels),
            np.asarray(samples), class_response)


def test_soft_segment_mixture_predicts_withheld_means_and_balances_specimen_means():
    Y, expression, anatomy, specimens, expected = _mixture_data()
    model = fit_segment_mixture(Y, expression, anatomy, specimens, C=1.)
    # Equal average of specimen-specific means cancels opposite sample shifts.
    np.testing.assert_allclose(model['gene_means'], expected)
    query = np.array([[-4., 0.], [0., 4.], [4., 0.]])
    out = predict_segment_mixture(model, query)
    np.testing.assert_allclose(out['probabilities'].sum(axis=1), 1.)
    assert np.all(np.diag(out['probabilities']) > .9)
    np.testing.assert_allclose(out['prediction'], expected, atol=.5)


def test_query_projection_has_no_labels_or_responses_and_does_not_mutate_model():
    Y, expression, anatomy, specimens, _ = _mixture_data()
    model = fit_segment_mixture(Y, expression, anatomy, specimens)
    before = model['gene_means'].copy()
    result = predict_segment_mixture(model, Y[:4])
    assert result['prediction'].shape == (4, expression.shape[1])
    np.testing.assert_array_equal(model['gene_means'], before)


def test_segment_mixture_requires_all_training_classes_and_valid_features():
    Y, expression, anatomy, specimens, _ = _mixture_data()
    with pytest.raises(ValueError, match='all labels'):
        fit_segment_mixture(Y, expression, np.where(anatomy == 2, 1, anatomy), specimens)
    model = fit_segment_mixture(Y, expression, anatomy, specimens)
    with pytest.raises(ValueError, match='feature count'):
        predict_segment_mixture(model, np.zeros((3, 4)))
