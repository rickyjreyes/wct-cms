from __future__ import annotations

import math
from typing import Optional

import numpy as np
from scipy.stats import chi2

from .models import ScanResult


def _as_1d_float(name: str, values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional array")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def weighted_linear_sinusoid(
    x: np.ndarray,
    y: np.ndarray,
    omega: float,
    weights: Optional[np.ndarray] = None,
) -> ScanResult:
    x = _as_1d_float("x", x)
    y = _as_1d_float("y", y)
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    if len(x) < 3:
        raise ValueError("at least three samples are required")

    if weights is None:
        weights = np.ones_like(y, dtype=float)
    else:
        weights = _as_1d_float("weights", weights)
        if len(weights) != len(y):
            raise ValueError("weights must have the same length as y")
        if np.any(weights < 0.0):
            raise ValueError("weights must be non-negative")

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


def _prepare_unweighted_basis(
    x: np.ndarray,
    omegas: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Precompute centered sinusoid bases and their 2x2 pseudoinverses.

    Including an intercept in a linear sinusoid fit is equivalent to centering
    the cosine/sine columns and the response, then solving only for the two
    oscillatory coefficients.  This lets an entire omega scan use matrix
    products instead of hundreds or thousands of repeated least-squares calls.
    """
    x = _as_1d_float("x", x)
    omegas = _as_1d_float("omegas", omegas)
    if len(omegas) == 0:
        raise ValueError("omegas must not be empty")

    phase = omegas[:, None] * x[None, :]
    c = np.cos(phase)
    s = np.sin(phase)
    c -= np.mean(c, axis=1, keepdims=True)
    s -= np.mean(s, axis=1, keepdims=True)

    g11 = np.einsum("ij,ij->i", c, c)
    g22 = np.einsum("ij,ij->i", s, s)
    g12 = np.einsum("ij,ij->i", c, s)

    gram = np.empty((len(omegas), 2, 2), dtype=float)
    gram[:, 0, 0] = g11
    gram[:, 0, 1] = g12
    gram[:, 1, 0] = g12
    gram[:, 1, 1] = g22
    gram_pinv = np.linalg.pinv(gram, rcond=1e-12)
    return c, s, gram_pinv


def _scores_from_basis(
    y: np.ndarray,
    c: np.ndarray,
    s: np.ndarray,
    gram_pinv: np.ndarray,
) -> np.ndarray:
    y = _as_1d_float("y", y)
    yc = y - np.mean(y)
    h = np.column_stack((c @ yc, s @ yc))
    scores = np.einsum("ki,kij,kj->k", h, gram_pinv, h)
    return np.maximum(scores, 0.0)


def _scores_from_basis_batch(
    ys: np.ndarray,
    c: np.ndarray,
    s: np.ndarray,
    gram_pinv: np.ndarray,
) -> np.ndarray:
    ys = np.asarray(ys, dtype=float)
    if ys.ndim != 2:
        raise ValueError("ys must be a two-dimensional array")
    if ys.shape[1] != c.shape[1]:
        raise ValueError("permutation sample length does not match x")

    yc = ys - np.mean(ys, axis=1, keepdims=True)
    h1 = yc @ c.T
    h2 = yc @ s.T
    h = np.stack((h1, h2), axis=-1)
    scores = np.einsum("bki,kij,bkj->bk", h, gram_pinv, h)
    return np.maximum(scores, 0.0)


def scan_omegas(
    x: np.ndarray,
    y: np.ndarray,
    omegas: np.ndarray,
) -> tuple[np.ndarray, ScanResult]:
    x = _as_1d_float("x", x)
    y = _as_1d_float("y", y)
    omegas = _as_1d_float("omegas", omegas)
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")

    c, s, gram_pinv = _prepare_unweighted_basis(x, omegas)
    scores = _scores_from_basis(y, c, s, gram_pinv)
    best_index = int(np.argmax(scores))

    # Use the original least-squares implementation once at the winning
    # frequency so amplitudes/phases retain the same definition as the fixed
    # frequency test while the expensive full scan stays vectorized.
    best = weighted_linear_sinusoid(x, y, float(omegas[best_index]))
    scores[best_index] = best.delta_chi2
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

    x = _as_1d_float("x", x)
    residuals = _as_1d_float("residuals", residuals)
    omegas = _as_1d_float("omegas", omegas)
    if len(x) != len(residuals):
        raise ValueError("x and residuals must have the same length")

    rng = np.random.default_rng(seed)
    max_scores = np.empty(n_perm, dtype=float)
    frozen_scores = np.empty(n_perm, dtype=float) if frozen_omega is not None else None

    c, s, gram_pinv = _prepare_unweighted_basis(x, omegas)
    frozen_basis = None
    if frozen_omega is not None:
        frozen_basis = _prepare_unweighted_basis(
            x, np.asarray([float(frozen_omega)], dtype=float)
        )

    # Batch the permutations so each batch uses BLAS matrix products instead
    # of running one least-squares solve for every omega of every permutation.
    batch_size = 128
    for start in range(0, n_perm, batch_size):
        stop = min(start + batch_size, n_perm)
        batch = np.stack(
            [rng.permutation(residuals) for _ in range(stop - start)], axis=0
        )
        batch_scores = _scores_from_basis_batch(batch, c, s, gram_pinv)
        max_scores[start:stop] = np.max(batch_scores, axis=1)

        if frozen_scores is not None and frozen_basis is not None:
            fc, fs, fpinv = frozen_basis
            frozen_scores[start:stop] = _scores_from_basis_batch(
                batch, fc, fs, fpinv
            )[:, 0]

    p_global = (np.count_nonzero(max_scores >= observed_best_score) + 1.0) / (n_perm + 1.0)

    p_frozen = None
    if frozen_scores is not None and observed_frozen_score is not None:
        p_frozen = (
            np.count_nonzero(frozen_scores >= observed_frozen_score) + 1.0
        ) / (n_perm + 1.0)

    return float(p_global), None if p_frozen is None else float(p_frozen), max_scores, frozen_scores
