import numpy as np

from cms_wct.signature import (
    permutation_null,
    scan_omegas,
    weighted_linear_sinusoid,
)


def test_recovers_injected_log_frequency():
    rng = np.random.default_rng(42)
    x = np.linspace(0.0, 4.0, 400)
    omega_true = 11.0
    y = 1.2 * np.cos(omega_true * x - 0.4) + rng.normal(0.0, 0.25, len(x))
    omegas = np.linspace(8.0, 14.0, 1201)
    _, best = scan_omegas(x, y, omegas)
    assert abs(best.omega - omega_true) < 0.05
    assert best.delta_chi2 > 100.0


def test_fixed_frequency_is_stronger_than_wrong_frequency():
    x = np.linspace(0.0, 3.0, 300)
    y = np.cos(7.0 * x)
    right = weighted_linear_sinusoid(x, y, 7.0)
    wrong = weighted_linear_sinusoid(x, y, 4.0)
    assert right.delta_chi2 > wrong.delta_chi2


def test_vectorized_scan_matches_scalar_least_squares():
    rng = np.random.default_rng(7)
    x = np.linspace(-1.0, 2.0, 180)
    y = rng.normal(size=len(x))
    omegas = np.linspace(0.5, 20.0, 80)

    scores, _ = scan_omegas(x, y, omegas)
    scalar_scores = np.array([
        weighted_linear_sinusoid(x, y, float(omega)).delta_chi2
        for omega in omegas
    ])
    np.testing.assert_allclose(scores, scalar_scores, rtol=1e-9, atol=1e-9)


def test_permutation_null_returns_reproducible_arrays():
    rng = np.random.default_rng(123)
    x = np.linspace(0.0, 3.0, 120)
    y = rng.normal(size=len(x))
    omegas = np.linspace(1.0, 15.0, 60)
    scores, best = scan_omegas(x, y, omegas)
    frozen = weighted_linear_sinusoid(x, y, 7.0)

    result1 = permutation_null(
        x, y, omegas, best.delta_chi2, 7.0, frozen.delta_chi2, 12, 99
    )
    result2 = permutation_null(
        x, y, omegas, best.delta_chi2, 7.0, frozen.delta_chi2, 12, 99
    )

    p_global, p_frozen, null_max, null_frozen = result1
    assert 0.0 < p_global <= 1.0
    assert p_frozen is not None and 0.0 < p_frozen <= 1.0
    assert null_max.shape == (12,)
    assert null_frozen is not None and null_frozen.shape == (12,)
    np.testing.assert_array_equal(result1[2], result2[2])
    np.testing.assert_array_equal(result1[3], result2[3])
