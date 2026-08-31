import numpy as np

from cms_wct.coherence_fast import paired_permutation_coherence_null_fast
from cms_wct.locked import (
    circular_mean,
    common_waveform_coherence,
    fit_phase_locked_waveform,
    paired_permutation_coherence_null,
    phase_locked_permutation_null,
)


def test_phase_locked_recovers_positive_amplitude():
    x = np.linspace(0.0, 4.0, 400)
    omega = 7.0
    phase = -0.2
    y = 0.9 * np.cos(omega * x - phase)
    result = fit_phase_locked_waveform(x, y, omega, phase)
    assert result.signed_amplitude > 0.89
    assert result.delta_chi2 > 100.0
    assert result.local_p_one_sided_chi2_1dof < 1e-10


def test_phase_locked_rejects_wrong_sign():
    x = np.linspace(0.0, 4.0, 400)
    omega = 7.0
    phase = -0.2
    y = -0.9 * np.cos(omega * x - phase)
    result = fit_phase_locked_waveform(x, y, omega, phase, positive_only=True)
    assert result.signed_amplitude < 0.0
    assert result.positive_amplitude == 0.0
    assert result.delta_chi2 == 0.0
    assert result.local_p_one_sided_chi2_1dof == 1.0


def test_phase_locked_permutation_is_reproducible():
    rng = np.random.default_rng(1)
    x = np.linspace(0.0, 3.0, 120)
    y = 0.6 * np.cos(5.0 * x + 0.3) + rng.normal(0.0, 0.4, len(x))
    obs = fit_phase_locked_waveform(x, y, 5.0, -0.3)
    a = phase_locked_permutation_null(x, y, 5.0, -0.3, obs.delta_chi2, 20, 99)
    b = phase_locked_permutation_null(x, y, 5.0, -0.3, obs.delta_chi2, 20, 99)
    assert a[0] == b[0]
    np.testing.assert_array_equal(a[1], b[1])


def test_circular_mean_close_phases():
    mean = circular_mean([-0.30599105090079975, -0.15679231968563603])
    assert abs(mean - (-0.2313916852932179)) < 1e-12


def test_common_waveform_rewards_coherence():
    x = np.linspace(0.0, 4.0, 300)
    omega = 7.0
    y1 = np.cos(omega * x - 0.2)
    y2 = 0.95 * np.cos(omega * x - 0.22)
    coherent = common_waveform_coherence(x, y1, x, y2, omega)
    anti = common_waveform_coherence(x, y1, x, -y2, omega)
    assert coherent["coherence_score"] > anti["coherence_score"]
    assert abs(coherent["phase_difference"]) < 0.05


def test_paired_coherence_null_shapes():
    rng = np.random.default_rng(42)
    x = np.linspace(0.0, 2.0, 80)
    y1 = rng.normal(size=len(x))
    y2 = rng.normal(size=len(x))
    observed = common_waveform_coherence(x, y1, x, y2, 4.0)
    p, null = paired_permutation_coherence_null(
        x, y1, x, y2, 4.0, observed["coherence_score"], 12, 123
    )
    assert p is not None and 0.0 < p <= 1.0
    assert null.shape == (12,)


def test_fast_paired_coherence_null_matches_scalar_seeded_draws():
    rng = np.random.default_rng(314)
    x1 = np.linspace(0.0, 2.0, 60)
    x2 = np.linspace(0.1, 2.2, 70)
    y1 = rng.normal(size=len(x1))
    y2 = rng.normal(size=len(x2))
    observed = common_waveform_coherence(x1, y1, x2, y2, 5.0)

    p_slow, n_slow = paired_permutation_coherence_null(
        x1, y1, x2, y2, 5.0, observed["coherence_score"], 23, 777
    )
    p_fast, n_fast = paired_permutation_coherence_null_fast(
        x1, y1, x2, y2, 5.0, observed["coherence_score"], 23, 777, batch_size=7
    )

    assert p_fast == p_slow
    np.testing.assert_allclose(n_fast, n_slow, rtol=1e-12, atol=1e-12)
