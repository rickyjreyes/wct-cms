from __future__ import annotations

import math
from typing import Optional

import numpy as np
from scipy.stats import chi2

from .models import ScanResult


def weighted_linear_sinusoid(
    x: np.ndarray,
    y: np.ndarray,
    omega: float,
    weights: Optional[np.ndarray] = None,
) -> ScanResult:
    if weights is None:
        weights = np.ones_like(y, dtype=float)

    sw = np.sqrt(weights)
    x0 = np.ones((len(x), 1), dtype=float)
    x1 = np.column_stack(
        [np.ones_like(x), np.cos(omega * x), np.sin(omega * x)]
    )

    x0w = x0 * sw[:, None]
    x1w = x1 * sw[:, None]
    yw = y * sw

    beta0, *_ = np.linalg.lstsq(x0w, yw, rcond=None)
    beta1, *_ = np.linalg.lstsq(x1w, yw, rcond=None)

    r0 = yw - x0w @ beta0
    r1 = yw - x1w @ beta1
    delta = max(0.0, float(r0 @ r0 - r1 @ r1))

    c0, a, b = map(float, beta1)
    return ScanResult(
        omega=float(omega),
        amplitude=math.hypot(a, b),
        phase=math.atan2(b, a),
        delta_chi2=delta,
        local_p_chi2_2dof=float(chi2.sf(delta, df=2)),
        c0=c0,
        cos_coeff=a,
        sin_coeff=b,
    )


def scan_omegas(
    x: np.ndarray,
    y: np.ndarray,
    omegas: np.ndarray,
) -> tuple[np.ndarray, ScanResult]:
    scores = np.empty_like(omegas, dtype=float)
    best = None
    for i, omega in enumerate(omegas):
        result = weighted_linear_sinusoid(x, y, float(omega))
        scores[i] = result.delta_chi2
        if best is None or result.delta_chi2 > best.delta_chi2:
            best = result
    assert best is not None
    return scores, best


def permutation_null(
    x: np.ndarray,
    residuals: np.ndarray,
    omegas: np.ndarray,
    observed_best_score: float,
    frozen_omega: Optional[float],
    observed_frozen_score: Optional[float],
    n_perm: int,
    seed: int,
) -> tuple[Optional[float], Optional[float], np.ndarray, Optional[np.ndarray]]:
    if n_perm <= 0:
        return None, None, np.array([]), None

    rng = np.random.default_rng(seed)
    max_scores = np.empty(n_perm, dtype=float)
    frozen_scores = np.empty(n_perm, dtype=float) if frozen_omega is not None else None

    for i in range(n_perm):
        shuffled = rng.permutation(residuals)
        scores, _ = scan_omegas(x, shuffled, omegas)
        max_scores[i] = float(np.max(scores))
        if frozen_scores is not None and frozen_omega is not None:
            frozen_scores[i] = weighted_linear_sinusoid(x, shuffled, frozen_omega).delta_chi2

    p_global = (np.count_nonzero(max_scores >= observed_best_score) + 1.0) / (n_perm + 1.0)

    p_frozen = None
    if frozen_scores is not None and observed_frozen_score is not None:
        p_frozen = (
            np.count_nonzero(frozen_scores >= observed_frozen_score) + 1.0
        ) / (n_perm + 1.0)

    return float(p_global), None if p_frozen is None else float(p_frozen), max_scores, frozen_scores
