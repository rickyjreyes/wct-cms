from __future__ import annotations

import numpy as np


def _as_1d(name: str, values) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def _basis(x: np.ndarray, omega: float):
    x = _as_1d("x", x)
    c = np.cos(float(omega) * x)
    s = np.sin(float(omega) * x)
    c -= np.mean(c)
    s -= np.mean(s)
    X = np.column_stack((c, s))
    gram = X.T @ X
    pinv = np.linalg.pinv(gram, rcond=1e-12)
    return X, gram, pinv


def paired_permutation_coherence_null_fast(
    x_a,
    y_a,
    x_b,
    y_b,
    omega: float,
    observed_coherence_score: float,
    n_perm: int,
    seed: int,
    *,
    batch_size: int = 512,
):
    """Vectorized paired residual-permutation null for the frozen H/G statistic.

    This is algebraically the same statistic as
    ``locked.paired_permutation_coherence_null`` but precomputes the harmonic
    Gram matrices and evaluates permutations in BLAS-friendly batches.  No
    analysis definition changes: H and G are still independently permuted in
    every paired trial, and

        T_coh = 2 * common_score - separate_score.
    """
    if n_perm <= 0:
        return None, np.array([], dtype=float)
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")

    y_a = _as_1d("y_a", y_a)
    y_b = _as_1d("y_b", y_b)
    Xa, ga, pa = _basis(x_a, omega)
    Xb, gb, pb = _basis(x_b, omega)
    if Xa.shape[0] != y_a.size or Xb.shape[0] != y_b.size:
        raise ValueError("x/y lengths differ")

    # The residual vectors supplied by the runner are centered already, but
    # centering here makes the function safe for direct use as well. A
    # permutation preserves the mean, so no per-trial recentering is needed.
    ya = y_a - np.mean(y_a)
    yb = y_b - np.mean(y_b)

    common_pinv = np.linalg.pinv(ga + gb, rcond=1e-12)
    rng = np.random.default_rng(seed)
    scores = np.empty(n_perm, dtype=float)

    for start in range(0, n_perm, batch_size):
        stop = min(start + batch_size, n_perm)
        n = stop - start
        ba = np.stack([rng.permutation(ya) for _ in range(n)], axis=0)
        bb = np.stack([rng.permutation(yb) for _ in range(n)], axis=0)

        ha = ba @ Xa
        hb = bb @ Xb
        h = ha + hb

        score_a = np.einsum("bi,ij,bj->b", ha, pa, ha)
        score_b = np.einsum("bi,ij,bj->b", hb, pb, hb)
        common = np.einsum("bi,ij,bj->b", h, common_pinv, h)

        separate = np.maximum(score_a, 0.0) + np.maximum(score_b, 0.0)
        common = np.maximum(common, 0.0)
        scores[start:stop] = 2.0 * common - separate

    p = (np.count_nonzero(scores >= float(observed_coherence_score)) + 1.0) / (n_perm + 1.0)
    return float(p), scores
