import numpy as np

from cms_wct.signature import scan_omegas, weighted_linear_sinusoid


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
