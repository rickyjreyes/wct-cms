from __future__ import annotations

from typing import Optional

import numpy as np

from .background import robust_smooth_background
from .signature import (
    _prepare_unweighted_basis,
    _scores_from_basis,
)


def parametric_background_null(
    centers: np.ndarray,
    null_model: np.ndarray,
    analysis_mask: np.ndarray,
    x: np.ndarray,
    omegas: np.ndarray,
    observed_best_score: float,
    frozen_omega: Optional[float],
    observed_frozen_score: Optional[float],
    degree: int,
    iterations: int,
    clip_sigma: float,
    excluded_windows,
    n_bootstrap: int,
    seed: int,
) -> tuple[Optional[float], Optional[float], np.ndarray, Optional[np.ndarray]]:
    """Parametric Poisson null with a complete background refit per trial.

    Each pseudo-spectrum is generated from the observed smooth-background model,
    then passed back through the same robust background fitter before its
    residual spectrum is scanned.  This propagates Poisson counting noise and
    the effect of re-estimating the continuum into the null distribution.

    The observed analysis mask is kept fixed.  Bin edges, resonance exclusions,
    and the minimum-model-count decision are analysis choices and should not
    fluctuate opportunistically from pseudo-experiment to pseudo-experiment.
    """
    if n_bootstrap <= 0:
        return None, None, np.array([]), None

    centers = np.asarray(centers, dtype=float)
    null_model = np.asarray(null_model, dtype=float)
    analysis_mask = np.asarray(analysis_mask, dtype=bool)
    x = np.asarray(x, dtype=float)
    omegas = np.asarray(omegas, dtype=float)

    if centers.ndim != 1 or null_model.ndim != 1 or analysis_mask.ndim != 1:
        raise ValueError("centers, null_model, and analysis_mask must be 1D")
    if not (len(centers) == len(null_model) == len(analysis_mask)):
        raise ValueError("centers, null_model, and analysis_mask lengths differ")
    if len(x) != int(np.count_nonzero(analysis_mask)):
        raise ValueError("x length must equal the number of analysis bins")
    if np.any(~np.isfinite(null_model)) or np.any(null_model < 0.0):
        raise ValueError("null_model must be finite and non-negative")

    rng = np.random.default_rng(seed)
    c, s, gram_pinv = _prepare_unweighted_basis(x, omegas)

    frozen_basis = None
    if frozen_omega is not None:
        frozen_basis = _prepare_unweighted_basis(
            x, np.asarray([float(frozen_omega)], dtype=float)
        )

    max_scores = np.empty(n_bootstrap, dtype=float)
    frozen_scores = (
        np.empty(n_bootstrap, dtype=float) if frozen_omega is not None else None
    )

    for i in range(n_bootstrap):
        pseudo_counts = rng.poisson(null_model).astype(float)
        pseudo_model, _ = robust_smooth_background(
            centers,
            pseudo_counts,
            degree=degree,
            iterations=iterations,
            clip_sigma=clip_sigma,
            excluded_windows=excluded_windows,
        )
        pseudo_residuals = (
            pseudo_counts - pseudo_model
        ) / np.sqrt(np.maximum(pseudo_model, 1.0))
        y = pseudo_residuals[analysis_mask]

        scores = _scores_from_basis(y, c, s, gram_pinv)
        max_scores[i] = float(np.max(scores))

        if frozen_scores is not None and frozen_basis is not None:
            fc, fs, fpinv = frozen_basis
            frozen_scores[i] = float(
                _scores_from_basis(y, fc, fs, fpinv)[0]
            )

    p_global = (
        np.count_nonzero(max_scores >= observed_best_score) + 1.0
    ) / (n_bootstrap + 1.0)

    p_frozen = None
    if frozen_scores is not None and observed_frozen_score is not None:
        p_frozen = (
            np.count_nonzero(frozen_scores >= observed_frozen_score) + 1.0
        ) / (n_bootstrap + 1.0)

    return (
        float(p_global),
        None if p_frozen is None else float(p_frozen),
        max_scores,
        frozen_scores,
    )
