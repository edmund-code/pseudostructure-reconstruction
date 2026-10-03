"""Synthetic checks for within-segment null and covariate oracle tools."""
import numpy as np
import pytest

from pseudospace.within_segment_validation import (
    fit_segment_covariate_oracle, permute_positions_within_segments,
    predict_segment_covariate_oracle,
)


def test_position_permutation_preserves_each_specimen_segment_multiset():
    z = np.linspace(.03, .97, 36)
    specimens = np.repeat(['A', 'B', 'C'], 12)
    anatomy = np.tile(np.repeat([0, 1, 2], 4), 3)
    original = z.copy()
    shuffled = permute_positions_within_segments(z, specimens, anatomy, seed=21)
    shuffled_again = permute_positions_within_segments(z, specimens, anatomy, seed=21)
    np.testing.assert_array_equal(shuffled, shuffled_again)
    np.testing.assert_array_equal(z, original)
    for specimen in np.unique(specimens):
        for segment in np.unique(anatomy):
            use = (specimens == specimen) & (anatomy == segment)
            np.testing.assert_array_equal(np.sort(shuffled[use]), np.sort(z[use]))


def _oracle_data():
    rows = []
    for specimen, counts in [('A', [20, 12, 8]), ('B', [8, 15, 20])]:
        for segment, n in enumerate(counts):
            x = np.linspace(-1, 1, n)
            cov = np.column_stack([x, x**2 + .1*x])
            intercept = np.array([2., -1.]) + segment*np.array([3., 2.])
            beta = np.array([[1.5, -.75], [.4, 1.1]])
            expr = intercept + cov @ beta
            rows.append((specimen, segment, cov, expr))
    specimens, anatomy, covariates, expression = [], [], [], []
    for specimen, segment, cov, expr in rows:
        specimens.extend([specimen]*len(cov))
        anatomy.extend([segment]*len(cov))
        covariates.append(cov)
        expression.append(expr)
    return np.vstack(expression), np.vstack(covariates), np.asarray(anatomy), np.asarray(specimens)


def test_segment_oracle_recovers_covariate_shift_and_is_specimen_balanced():
    expression, covariates, anatomy, specimens = _oracle_data()
    model = fit_segment_covariate_oracle(expression, covariates, anatomy,
                                         specimens, ridge=1e-7)
    query_cov = np.tile(np.array([[-.5, .3], [.6, .4]]), (3, 1))
    query_labels = np.repeat([0, 1, 2], 2)
    expected = []
    for segment in range(3):
        expected.extend((np.array([2., -1.]) + segment*np.array([3., 2.])
                         + query_cov[2*segment:2*segment+2]
                         @ np.array([[1.5, -.75], [.4, 1.1]])).tolist())
    prediction = predict_segment_covariate_oracle(model, query_cov, query_labels)
    np.testing.assert_allclose(prediction, np.asarray(expected), atol=2e-5)

    # A genuine specimen shift makes this sensitive to accidental pooling weights.
    shifted = expression + np.where(specimens == 'A', 4., -2.)[:, None]
    model = fit_segment_covariate_oracle(shifted, covariates, anatomy, specimens, ridge=.01)
    prediction = predict_segment_covariate_oracle(model, query_cov, query_labels)
    for segment in range(3):
        means = [shifted[(anatomy == segment) & (specimens == s)].mean(axis=0) for s in ['A','B']]
        np.testing.assert_allclose(model['segments'][segment]['intercept'], np.mean(means, axis=0))
    # Repeating all rows from one original specimen must not increase its mass.
    duplicate = np.flatnonzero(specimens == 'A')
    repeated_model = fit_segment_covariate_oracle(
        np.vstack([shifted, shifted[duplicate]]),
        np.vstack([covariates, covariates[duplicate]]),
        np.r_[anatomy, anatomy[duplicate]],
        np.r_[specimens, specimens[duplicate]], ridge=.01)
    repeated_prediction = predict_segment_covariate_oracle(
        repeated_model, query_cov, query_labels)
    np.testing.assert_allclose(repeated_prediction, prediction, atol=1e-10)


def test_permutation_and_oracle_reject_invalid_inputs_or_missing_segments():
    with pytest.raises(ValueError, match='finite'):
        permute_positions_within_segments([0., np.nan], ['a', 'a'], [0, 0], seed=1)
    with pytest.raises(ValueError, match='specimens'):
        permute_positions_within_segments([.1, .2], ['a'], [0, 0], seed=1)
    with pytest.raises(ValueError, match='nonmissing'):
        permute_positions_within_segments([.1, .2], ['a', None], [0, 0], seed=1)

    with pytest.raises(ValueError, match='anatomy'):
        permute_positions_within_segments([.1,.2], ['a','a'], [0,np.inf], seed=1)

    expression = np.arange(12, dtype=float).reshape(6, 2)
    covariates = np.arange(6, dtype=float).reshape(6, 1)
    labels = np.array([0, 0, 1, 1, 1, 1])
    specimens = np.array(['a', 'b', 'a', 'a', 'b', 'b'])
    with pytest.raises(ValueError, match='lack segment 2'):
        fit_segment_covariate_oracle(expression, covariates, labels, specimens)
    expression, covariates, anatomy, specimens = _oracle_data()
    model = fit_segment_covariate_oracle(expression, covariates, anatomy, specimens)
    with pytest.raises(ValueError, match='feature count'):
        predict_segment_covariate_oracle(model, np.zeros((2, 4)), [0, 1])


def _quadratic_data(seed=62, n_per_group=100):
    rng = np.random.default_rng(seed)
    expression, covariates, anatomy, specimens = [], [], [], []
    for specimen in ('A', 'B'):
        for segment in (0, 1, 2):
            cov = rng.uniform(-1, 1, size=(n_per_group, 2))
            x, y = cov.T
            intercept = np.array([2., -1.]) + segment*np.array([2., 3.])
            response = np.column_stack([
                intercept[0] + 1.2*x - .7*y + 1.5*x*x + .8*x*y - .4*y*y,
                intercept[1] - .3*x + 1.1*y - .5*x*x + .6*x*y + 1.3*y*y,
            ])
            response += rng.normal(0, .005, response.shape)
            expression.append(response)
            covariates.append(cov)
            anatomy.extend([segment]*n_per_group)
            specimens.extend([specimen]*n_per_group)
    return np.vstack(expression), np.vstack(covariates), np.asarray(anatomy), np.asarray(specimens)


def test_degree_two_recovers_quadratic_expression_at_shifted_queries():
    expression, covariates, anatomy, specimens = _quadratic_data()
    model = fit_segment_covariate_oracle(expression, covariates, anatomy,
                                         specimens, degree=2)
    query_cov = np.tile(np.array([[.62, -.37], [-.54, .46]]), (3, 1))
    query_labels = np.repeat([0, 1, 2], 2)
    x, y = query_cov.T
    expected = []
    for i, segment in enumerate(query_labels):
        base = np.array([2., -1.]) + segment*np.array([2., 3.])
        expected.append([base[0]+1.2*x[i]-.7*y[i]
                         +1.5*x[i]**2+.8*x[i]*y[i]
                         -.4*y[i]**2,
                         base[1]-.3*x[i]+1.1*y[i]
                         -.5*x[i]**2+.6*x[i]*y[i]
                         +1.3*y[i]**2])
    prediction = predict_segment_covariate_oracle(model, query_cov, query_labels)
    np.testing.assert_allclose(prediction, expected, atol=.04)


def test_degree_one_default_is_unchanged_and_degree_two_is_specimen_balanced():
    expression, covariates, anatomy, specimens = _oracle_data()
    default = fit_segment_covariate_oracle(expression, covariates, anatomy, specimens)
    explicit = fit_segment_covariate_oracle(expression, covariates, anatomy, specimens, degree=1)
    query = np.array([[-.3, .5], [.4, .2], [0., .1]])
    np.testing.assert_array_equal(
        predict_segment_covariate_oracle(default, query, [0, 1, 2]),
        predict_segment_covariate_oracle(explicit, query, [0, 1, 2]))

    shifted = expression.copy()
    shifted[specimens == 'A'] += [4., -3.]
    model = fit_segment_covariate_oracle(shifted, covariates, anatomy,
                                         specimens, degree=2)
    duplicate = np.flatnonzero(specimens == 'A')
    repeated = fit_segment_covariate_oracle(
        np.vstack([shifted, shifted[duplicate]]),
        np.vstack([covariates, covariates[duplicate]]),
        np.r_[anatomy, anatomy[duplicate]],
        np.r_[specimens, specimens[duplicate]], degree=2)
    np.testing.assert_allclose(
        predict_segment_covariate_oracle(model, query, [0, 1, 2]),
        predict_segment_covariate_oracle(repeated, query, [0, 1, 2]), atol=1e-10)


@pytest.mark.parametrize('degree', [0, 3, 1.5, True])
def test_oracle_rejects_invalid_polynomial_degrees(degree):
    expression, covariates, anatomy, specimens = _oracle_data()
    with pytest.raises(ValueError, match='degree'):
        fit_segment_covariate_oracle(expression, covariates, anatomy,
                                     specimens, degree=degree)
