from __future__ import annotations

from dataclasses import dataclass
from math import comb

import numpy as np
from numpy.polynomial import Chebyshev
from scipy.interpolate import UnivariateSpline


@dataclass(frozen=True)
class BackgroundCandidate:
    name: str
    family: str
    degree: int = 7
    spline_smoothing_factor: float = 1.0

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "family": self.family,
            "degree": int(self.degree),
            "spline_smoothing_factor": float(self.spline_smoothing_factor),
        }


def frozen_background_candidates() -> list[BackgroundCandidate]:
    """Continuum candidates frozen before the H/G CV result is inspected."""
    out = [BackgroundCandidate(f"cheb_d{d}", "chebyshev", degree=d) for d in range(5, 13)]
    out += [BackgroundCandidate(f"bernstein_d{d}", "bernstein", degree=d) for d in (5, 7, 9, 12)]
    out += [
        BackgroundCandidate(f"spline_s{s:g}", "spline", degree=7, spline_smoothing_factor=s)
        for s in (0.5, 1.0, 2.0)
    ]
    return out


def analysis_eligibility(centers: np.ndarray, excluded_windows=()) -> np.ndarray:
    centers = np.asarray(centers, dtype=float)
    eligible = np.isfinite(centers) & (centers > 0.0)
    for lo, hi in excluded_windows:
        eligible &= ~((centers >= float(lo)) & (centers <= float(hi)))
    return eligible


def blocked_fold_ids(eligible: np.ndarray, *, n_folds: int = 5, block_size: int = 16) -> np.ndarray:
    """Assign contiguous eligible stretches to cyclic blocked CV folds.

    Blocks never bridge an excluded gap. A block is a contiguous run of up to
    ``block_size`` eligible histogram bins. Block labels cycle 0..n_folds-1.
    """
    eligible = np.asarray(eligible, dtype=bool)
    if eligible.ndim != 1:
        raise ValueError("eligible must be one-dimensional")
    if n_folds < 2:
        raise ValueError("n_folds must be at least 2")
    if block_size < 1:
        raise ValueError("block_size must be positive")

    fold_ids = np.full(len(eligible), -1, dtype=int)
    idx = np.flatnonzero(eligible)
    if idx.size == 0:
        raise ValueError("no eligible bins")

    block_number = 0
    start = 0
    while start < idx.size:
        run_end = start + 1
        while run_end < idx.size and idx[run_end] == idx[run_end - 1] + 1:
            run_end += 1
        run = idx[start:run_end]
        for offset in range(0, len(run), block_size):
            block = run[offset : offset + block_size]
            fold_ids[block] = block_number % n_folds
            block_number += 1
        start = run_end

    return fold_ids


def poisson_deviance(observed: np.ndarray, expected: np.ndarray) -> float:
    y = np.asarray(observed, dtype=float)
    mu = np.asarray(expected, dtype=float)
    if y.shape != mu.shape:
        raise ValueError("observed and expected shapes differ")
    if np.any(y < 0.0) or np.any(~np.isfinite(y)):
        raise ValueError("observed counts must be finite and non-negative")
    mu = np.clip(mu, 1e-12, np.inf)
    term = np.where(y > 0.0, y * np.log(y / mu), 0.0)
    return float(2.0 * np.sum(mu - y + term))


def _bernstein_design(x: np.ndarray, degree: int) -> np.ndarray:
    xmin = float(np.min(x))
    xmax = float(np.max(x))
    if not xmax > xmin:
        raise ValueError("background coordinate has zero span")
    t = (x - xmin) / (xmax - xmin)
    return np.column_stack(
        [comb(degree, j) * np.power(t, j) * np.power(1.0 - t, degree - j) for j in range(degree + 1)]
    )


def _fit_linear_log_model(
    design: np.ndarray,
    counts: np.ndarray,
    train_mask: np.ndarray,
    *,
    minimum: int,
    iterations: int,
    clip_sigma: float,
) -> np.ndarray:
    y = np.log(counts + 0.5)
    fit_mask = train_mask & (counts > 0.0)
    if np.count_nonzero(fit_mask) < minimum:
        raise RuntimeError("Too few training bins for background fit")

    for _ in range(iterations):
        w = np.sqrt(np.maximum(counts[fit_mask], 1.0))
        beta, *_ = np.linalg.lstsq(design[fit_mask] * w[:, None], y[fit_mask] * w, rcond=None)
        model = np.exp(np.clip(design @ beta, -30.0, 50.0))
        residual = (counts - model) / np.sqrt(np.maximum(model, 1.0))
        new_mask = fit_mask & (np.abs(residual) <= clip_sigma)
        if np.array_equal(new_mask, fit_mask):
            break
        fit_mask = new_mask
        if np.count_nonzero(fit_mask) < minimum:
            raise RuntimeError("Robust training fit clipped too many bins")

    w = np.sqrt(np.maximum(counts[fit_mask], 1.0))
    beta, *_ = np.linalg.lstsq(design[fit_mask] * w[:, None], y[fit_mask] * w, rcond=None)
    return np.exp(np.clip(design @ beta, -30.0, 50.0))


def fit_candidate_train_mask(
    centers: np.ndarray,
    counts: np.ndarray,
    candidate: BackgroundCandidate,
    train_mask: np.ndarray,
    *,
    iterations: int = 6,
    clip_sigma: float = 3.5,
) -> np.ndarray:
    """Fit one continuum candidate using only training bins."""
    centers = np.asarray(centers, dtype=float)
    counts = np.asarray(counts, dtype=float)
    train_mask = np.asarray(train_mask, dtype=bool)
    if not (centers.ndim == counts.ndim == train_mask.ndim == 1):
        raise ValueError("centers, counts, train_mask must be one-dimensional")
    if not (len(centers) == len(counts) == len(train_mask)):
        raise ValueError("centers, counts, train_mask lengths differ")

    x = np.log(centers)
    family = candidate.family.lower()

    if family == "chebyshev":
        degree = int(candidate.degree)
        # Chebyshev.vander is not exposed through the convenience class, so
        # fit the coefficients on the training coordinate and evaluate on all bins.
        mask = train_mask & (counts > 0.0)
        if np.count_nonzero(mask) < degree + 3:
            raise RuntimeError("Too few training bins for Chebyshev fit")
        fit_mask = mask.copy()
        y = np.log(counts + 0.5)
        domain = [float(np.min(x)), float(np.max(x))]
        for _ in range(iterations):
            w = np.sqrt(np.maximum(counts[fit_mask], 1.0))
            poly = Chebyshev.fit(x[fit_mask], y[fit_mask], deg=degree, w=w, domain=domain)
            model = np.exp(np.clip(poly(x), -30.0, 50.0))
            residual = (counts - model) / np.sqrt(np.maximum(model, 1.0))
            new_mask = fit_mask & (np.abs(residual) <= clip_sigma)
            if np.array_equal(new_mask, fit_mask):
                break
            fit_mask = new_mask
            if np.count_nonzero(fit_mask) < degree + 3:
                raise RuntimeError("Robust Chebyshev training fit clipped too many bins")
        w = np.sqrt(np.maximum(counts[fit_mask], 1.0))
        poly = Chebyshev.fit(x[fit_mask], y[fit_mask], deg=degree, w=w, domain=domain)
        return np.exp(np.clip(poly(x), -30.0, 50.0))

    if family == "bernstein":
        degree = int(candidate.degree)
        design = _bernstein_design(x, degree)
        return _fit_linear_log_model(
            design,
            counts,
            train_mask,
            minimum=degree + 3,
            iterations=iterations,
            clip_sigma=clip_sigma,
        )

    if family == "spline":
        fit_mask = train_mask & (counts > 0.0)
        if np.count_nonzero(fit_mask) < 8:
            raise RuntimeError("Too few training bins for spline fit")
        y = np.log(counts + 0.5)

        def fit(mask: np.ndarray) -> np.ndarray:
            w = np.sqrt(np.maximum(counts[mask], 1.0))
            smoothing = float(candidate.spline_smoothing_factor) * float(np.count_nonzero(mask))
            spline = UnivariateSpline(x[mask], y[mask], w=w, k=3, s=smoothing, ext=0)
            return np.exp(np.clip(spline(x), -30.0, 50.0))

        for _ in range(iterations):
            model = fit(fit_mask)
            residual = (counts - model) / np.sqrt(np.maximum(model, 1.0))
            new_mask = fit_mask & (np.abs(residual) <= clip_sigma)
            if np.array_equal(new_mask, fit_mask):
                break
            fit_mask = new_mask
            if np.count_nonzero(fit_mask) < 8:
                raise RuntimeError("Robust spline training fit clipped too many bins")
        return fit(fit_mask)

    raise ValueError(f"Unknown background family: {candidate.family}")


def cross_validate_candidates(
    centers: np.ndarray,
    counts: np.ndarray,
    candidates: list[BackgroundCandidate],
    *,
    excluded_windows=(),
    n_folds: int = 5,
    block_size: int = 16,
    iterations: int = 6,
    clip_sigma: float = 3.5,
) -> list[dict]:
    """Score continuum models without using any WCT frequency/amplitude."""
    centers = np.asarray(centers, dtype=float)
    counts = np.asarray(counts, dtype=float)
    eligible = analysis_eligibility(centers, excluded_windows)
    folds = blocked_fold_ids(eligible, n_folds=n_folds, block_size=block_size)

    results = []
    for candidate in candidates:
        fold_deviances = []
        fold_bins = []
        failed = None
        try:
            for fold in range(n_folds):
                validation = eligible & (folds == fold)
                training = eligible & (folds != fold)
                if np.count_nonzero(validation) == 0:
                    continue
                model = fit_candidate_train_mask(
                    centers,
                    counts,
                    candidate,
                    training,
                    iterations=iterations,
                    clip_sigma=clip_sigma,
                )
                fold_deviances.append(poisson_deviance(counts[validation], model[validation]))
                fold_bins.append(int(np.count_nonzero(validation)))
        except Exception as exc:  # retained in the result rather than silently dropped
            failed = f"{type(exc).__name__}: {exc}"

        total_bins = int(sum(fold_bins))
        total_deviance = float(sum(fold_deviances)) if failed is None else float("inf")
        results.append(
            {
                **candidate.to_dict(),
                "failed": failed,
                "n_validation_bins": total_bins,
                "total_poisson_deviance": total_deviance,
                "deviance_per_bin": total_deviance / total_bins if total_bins > 0 else float("inf"),
                "fold_deviances": fold_deviances,
                "fold_bins": fold_bins,
            }
        )
    return results
