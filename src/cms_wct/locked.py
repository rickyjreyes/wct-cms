from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Iterable

import numpy as np
from scipy.stats import chi2


@dataclass
class PhaseLockedResult:
    omega: float
    phase: float
    c0: float
    signed_amplitude: float
    positive_amplitude: float
    delta_chi2: float
    local_p_one_sided_chi2_1dof: float

    def to_dict(self) -> dict:
        return asdict(self)


def _as_1d(name: str, values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def phase_locked_basis(x: np.ndarray, omega: float, phase: float) -> np.ndarray:
    x = _as_1d("x", x)
    return np.cos(float(omega) * x - float(phase))


def fit_phase_locked_waveform(
    x: np.ndarray,
    y: np.ndarray,
    omega: float,
    phase: float,
    *,
    positive_only: bool = True,
) -> PhaseLockedResult:
    """Fit c + A cos(omega*x - phase) with omega/phase fixed in advance.

    If positive_only=True, the preregistered alternative is A>0. Negative
    unconstrained amplitudes therefore receive zero improvement and p=1.
    For a positive fitted amplitude, the asymptotic one-sided fixed-waveform
    reference is the 50:50 boundary mixture, p = 0.5*chi2.sf(delta, 1).
    Empirical permutation/bootstrap calibration is preferred for claims.
    """
    x = _as_1d("x", x)
    y = _as_1d("y", y)
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    if len(x) < 3:
        raise ValueError("at least three samples are required")

    g = phase_locked_basis(x, omega, phase)
    x0 = np.ones((len(x), 1), dtype=float)
    x1 = np.column_stack((np.ones_like(x), g))

    beta0, *_ = np.linalg.lstsq(x0, y, rcond=None)
    beta1, *_ = np.linalg.lstsq(x1, y, rcond=None)
    r0 = y - x0 @ beta0
    r1 = y - x1 @ beta1

    c0 = float(beta1[0])
    amp = float(beta1[1])
    raw_delta = max(0.0, float(r0 @ r0 - r1 @ r1))

    if positive_only and amp <= 0.0:
        delta = 0.0
        p_local = 1.0
        positive_amp = 0.0
    else:
        delta = raw_delta
        positive_amp = max(0.0, amp) if positive_only else amp
        p_local = float(0.5 * chi2.sf(delta, df=1)) if positive_only else float(chi2.sf(delta, df=1))

    return PhaseLockedResult(
        omega=float(omega),
        phase=float(phase),
        c0=c0,
        signed_amplitude=amp,
        positive_amplitude=float(positive_amp),
        delta_chi2=delta,
        local_p_one_sided_chi2_1dof=p_local,
    )


def phase_locked_permutation_null(
    x: np.ndarray,
    y: np.ndarray,
    omega: float,
    phase: float,
    observed_delta: float,
    n_perm: int,
    seed: int,
    *,
    positive_only: bool = True,
) -> tuple[float | None, np.ndarray]:
    if n_perm <= 0:
        return None, np.array([], dtype=float)
    x = _as_1d("x", x)
    y = _as_1d("y", y)
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")

    rng = np.random.default_rng(seed)
    scores = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        yp = rng.permutation(y)
        scores[i] = fit_phase_locked_waveform(
            x, yp, omega, phase, positive_only=positive_only
        ).delta_chi2

    p = (np.count_nonzero(scores >= float(observed_delta)) + 1.0) / (n_perm + 1.0)
    return float(p), scores


def circular_mean(phases: Iterable[float]) -> float:
    vals = np.asarray(list(phases), dtype=float)
    if vals.size == 0 or not np.all(np.isfinite(vals)):
        raise ValueError("phases must contain finite values")
    z = np.mean(np.exp(1j * vals))
    if abs(z) < 1e-15:
        raise ValueError("circular mean undefined for cancelling phases")
    return float(np.angle(z))


def _centered_harmonic_components(x: np.ndarray, y: np.ndarray, omega: float):
    x = _as_1d("x", x)
    y = _as_1d("y", y)
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    yc = y - np.mean(y)
    c = np.cos(float(omega) * x)
    s = np.sin(float(omega) * x)
    c -= np.mean(c)
    s -= np.mean(s)
    X = np.column_stack((c, s))
    gram = X.T @ X
    h = X.T @ yc
    pinv = np.linalg.pinv(gram, rcond=1e-12)
    beta = pinv @ h
    score = max(0.0, float(h @ pinv @ h))
    return gram, h, beta, score


def common_waveform_coherence(
    x_a: np.ndarray,
    y_a: np.ndarray,
    x_b: np.ndarray,
    y_b: np.ndarray,
    omega: float,
) -> dict:
    """Retrospective H/G common-waveform diagnostic at fixed omega.

    The common-component likelihood-ratio score rewards aligned complex
    coefficients automatically. The heterogeneity score is the extra fit
    gained by allowing separate coefficients. The frozen coherence statistic
    is T = common_score - heterogeneity = 2*common_score - separate_score.
    It must be calibrated end-to-end with paired null spectra.
    """
    ga, ha, ba, sa = _centered_harmonic_components(x_a, y_a, omega)
    gb, hb, bb, sb = _centered_harmonic_components(x_b, y_b, omega)
    g = ga + gb
    h = ha + hb
    pinv = np.linalg.pinv(g, rcond=1e-12)
    beta_common = pinv @ h
    common_score = max(0.0, float(h @ pinv @ h))
    separate_score = float(sa + sb)
    heterogeneity = max(0.0, separate_score - common_score)
    coherence_score = float(common_score - heterogeneity)

    def amp_phase(beta):
        return float(math.hypot(beta[0], beta[1])), float(math.atan2(beta[1], beta[0]))

    amp_a, phase_a = amp_phase(ba)
    amp_b, phase_b = amp_phase(bb)
    amp_c, phase_c = amp_phase(beta_common)
    phase_difference = float(np.angle(np.exp(1j * (phase_a - phase_b))))

    return {
        "omega": float(omega),
        "common_score": common_score,
        "separate_score": separate_score,
        "heterogeneity": heterogeneity,
        "coherence_score": coherence_score,
        "amplitude_a": amp_a,
        "phase_a": phase_a,
        "amplitude_b": amp_b,
        "phase_b": phase_b,
        "common_amplitude": amp_c,
        "common_phase": phase_c,
        "phase_difference": phase_difference,
    }


def paired_permutation_coherence_null(
    x_a: np.ndarray,
    y_a: np.ndarray,
    x_b: np.ndarray,
    y_b: np.ndarray,
    omega: float,
    observed_coherence_score: float,
    n_perm: int,
    seed: int,
) -> tuple[float | None, np.ndarray]:
    if n_perm <= 0:
        return None, np.array([], dtype=float)
    rng = np.random.default_rng(seed)
    scores = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        pa = rng.permutation(y_a)
        pb = rng.permutation(y_b)
        scores[i] = common_waveform_coherence(x_a, pa, x_b, pb, omega)["coherence_score"]
    p = (np.count_nonzero(scores >= float(observed_coherence_score)) + 1.0) / (n_perm + 1.0)
    return float(p), scores
