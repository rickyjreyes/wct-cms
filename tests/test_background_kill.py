import numpy as np

from cms_wct.background_cv import BackgroundCandidate
from cms_wct.background_kill import (
    PairKillConfig,
    deterministic_absorption_matrix,
    evaluate_candidate,
    injected_mean,
    spurious_signal_null_pair,
)


def synthetic_pair():
    centers = np.geomspace(2.0, 120.0, 90)
    x = np.log(centers)
    smooth_a = np.exp(6.3 - 0.28 * x + 0.018 * x**2)
    smooth_b = np.exp(6.4 - 0.30 * x + 0.020 * x**2)
    return centers, smooth_a, smooth_b


def test_injected_mean_uses_residual_amplitude_units():
    centers, smooth_a, _ = synthetic_pair()
    config = PairKillConfig(omega=7.0, phase=-0.2, min_model_count=5.0)
    amplitude = 0.8
    injected = injected_mean(centers, smooth_a, config, amplitude)
    expected_residual = (injected - smooth_a) / np.sqrt(smooth_a)
    target = amplitude * np.cos(config.omega * np.log(centers) - config.phase)
    np.testing.assert_allclose(expected_residual, target, rtol=1e-12, atol=1e-12)


def test_candidate_recovers_positive_locked_injection():
    centers, smooth_a, _ = synthetic_pair()
    config = PairKillConfig(omega=7.0, phase=-0.2, min_model_count=5.0)
    candidate = BackgroundCandidate("cheb_d3", "chebyshev", degree=3)
    injected = injected_mean(centers, smooth_a, config, amplitude=0.8)
    result, _ = evaluate_candidate(centers, injected, candidate, config)
    assert result["positive_sign"]
    assert result["locked_signed_amplitude"] > 0.2
    assert result["locked_delta_chi2"] > 0.0


def test_spurious_signal_null_runs_full_pair_pipeline():
    centers, smooth_a, smooth_b = synthetic_pair()
    config = PairKillConfig(
        omega=7.0,
        phase=-0.2,
        n_folds=3,
        block_size=10,
        min_model_count=5.0,
    )
    candidates = [
        BackgroundCandidate("cheb_d2", "chebyshev", degree=2),
        BackgroundCandidate("cheb_d3", "chebyshev", degree=3),
    ]
    result = spurious_signal_null_pair(
        centers,
        smooth_a,
        centers,
        smooth_b,
        candidates,
        config,
        n_trials=3,
        seed=1234,
        reselect_each_trial=True,
    )
    assert result["n_trials"] == 3
    assert len(result["trials"]) == 3
    assert 1.0 / 4.0 <= result["primary_add_one_p"] <= 1.0
    assert sum(result["selected_background_counts"].values()) == 3
    assert np.isfinite(result["observed"]["primary_pair_score_min_locked_delta_chi2"])


def test_absorption_matrix_covers_every_candidate_and_amplitude():
    centers, smooth_a, smooth_b = synthetic_pair()
    config = PairKillConfig(
        omega=7.0,
        phase=-0.2,
        n_folds=3,
        block_size=10,
        min_model_count=5.0,
    )
    candidates = [
        BackgroundCandidate("cheb_d2", "chebyshev", degree=2),
        BackgroundCandidate("cheb_d3", "chebyshev", degree=3),
    ]
    result = deterministic_absorption_matrix(
        centers,
        smooth_a,
        centers,
        smooth_b,
        candidates,
        config,
        injection_amplitudes=[0.4, 0.8],
    )
    assert len(result["rows"]) == 4
    assert {row["analysis_candidate"] for row in result["rows"]} == {
        "cheb_d2",
        "cheb_d3",
    }
    assert {row["injected_amplitude"] for row in result["rows"]} == {0.4, 0.8}
    assert all(np.isfinite(row["min_retention"]) for row in result["rows"])
