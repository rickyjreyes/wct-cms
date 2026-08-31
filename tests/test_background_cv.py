import numpy as np

from cms_wct.background_cv import (
    BackgroundCandidate,
    analysis_eligibility,
    blocked_fold_ids,
    cross_validate_candidates,
    fit_candidate_train_mask,
    poisson_deviance,
)


def test_blocked_folds_do_not_bridge_excluded_gap():
    centers = np.arange(1.0, 31.0)
    eligible = analysis_eligibility(centers, excluded_windows=[(11.0, 15.0)])
    folds = blocked_fold_ids(eligible, n_folds=3, block_size=4)
    assert np.all(folds[10:15] == -1)
    assert np.all(folds[eligible] >= 0)


def test_training_fit_does_not_use_validation_counts():
    centers = np.geomspace(2.0, 120.0, 80)
    x = np.log(centers)
    counts = np.exp(5.0 - 0.4 * x + 0.03 * x**2)
    train = np.ones(len(counts), dtype=bool)
    train[25:40] = False
    candidate = BackgroundCandidate("cheb_d5", "chebyshev", degree=5)

    model_a = fit_candidate_train_mask(centers, counts, candidate, train)
    changed = counts.copy()
    changed[25:40] *= 50.0
    model_b = fit_candidate_train_mask(centers, changed, candidate, train)
    np.testing.assert_allclose(model_a, model_b, rtol=1e-12, atol=1e-12)


def test_poisson_deviance_is_zero_for_exact_prediction():
    y = np.array([0.0, 1.0, 4.0, 10.0])
    assert abs(poisson_deviance(y, y)) < 1e-12


def test_cross_validation_returns_all_candidates():
    centers = np.geomspace(2.0, 120.0, 100)
    x = np.log(centers)
    counts = np.exp(6.0 - 0.35 * x + 0.02 * x**2)
    candidates = [
        BackgroundCandidate("cheb_d5", "chebyshev", degree=5),
        BackgroundCandidate("bernstein_d5", "bernstein", degree=5),
        BackgroundCandidate("spline_s2", "spline", spline_smoothing_factor=2.0),
    ]
    rows = cross_validate_candidates(
        centers,
        counts,
        candidates,
        excluded_windows=[(8.5, 11.5)],
        n_folds=5,
        block_size=8,
    )
    assert [r["name"] for r in rows] == [c.name for c in candidates]
    assert all(r["n_validation_bins"] > 0 for r in rows)
    assert all(np.isfinite(r["total_poisson_deviance"]) for r in rows)
