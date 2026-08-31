from __future__ import annotations

from math import comb

import numpy as np
from scipy.interpolate import UnivariateSpline

from .background import robust_smooth_background


def _initial_mask(centers, counts, excluded_windows):
    mask = np.asarray(counts, dtype=float) > 0
    for lo, hi in excluded_windows:
        mask &= ~((centers >= lo) & (centers <= hi))
    return mask


def _bernstein_design(x: np.ndarray, degree: int) -> np.ndarray:
    xmin = float(np.min(x))
    xmax = float(np.max(x))
    if not xmax > xmin:
        raise ValueError("background coordinate has zero span")
    t = (x - xmin) / (xmax - xmin)
    cols = [comb(degree, j) * np.power(t, j) * np.power(1.0 - t, degree - j) for j in range(degree + 1)]
    return np.column_stack(cols)


def _robust_linear_log_background(
    centers,
    counts,
    design,
    minimum,
    iterations,
    clip_sigma,
    excluded_windows,
):
    counts = np.asarray(counts, dtype=float)
    y = np.log(counts + 0.5)
    fit_mask = _initial_mask(centers, counts, excluded_windows)
    if np.count_nonzero(fit_mask) < minimum:
        raise RuntimeError("Too few bins for background fit")

    for _ in range(iterations):
        w = np.sqrt(np.maximum(counts[fit_mask], 1.0))
        Xw = design[fit_mask] * w[:, None]
        yw = y[fit_mask] * w
        beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        model = np.exp(np.clip(design @ beta, -30, 50))
        residual = (counts - model) / np.sqrt(np.maximum(model, 1.0))
        new_mask = fit_mask & (np.abs(residual) <= clip_sigma)
        if np.array_equal(new_mask, fit_mask):
            break
        fit_mask = new_mask
        if np.count_nonzero(fit_mask) < minimum:
            raise RuntimeError("Robust fit clipped too many bins")

    w = np.sqrt(np.maximum(counts[fit_mask], 1.0))
    beta, *_ = np.linalg.lstsq(design[fit_mask] * w[:, None], y[fit_mask] * w, rcond=None)
    model = np.exp(np.clip(design @ beta, -30, 50))
    return model, fit_mask


def robust_bernstein_background(
    centers,
    counts,
    degree=7,
    iterations=6,
    clip_sigma=3.5,
    excluded_windows=(),
):
    x = np.log(np.asarray(centers, dtype=float))
    design = _bernstein_design(x, int(degree))
    return _robust_linear_log_background(
        np.asarray(centers, dtype=float),
        counts,
        design,
        int(degree) + 3,
        iterations,
        clip_sigma,
        excluded_windows,
    )


def robust_spline_background(
    centers,
    counts,
    iterations=6,
    clip_sigma=3.5,
    excluded_windows=(),
    smoothing_factor=1.0,
):
    centers = np.asarray(centers, dtype=float)
    counts = np.asarray(counts, dtype=float)
    x = np.log(centers)
    y = np.log(counts + 0.5)
    fit_mask = _initial_mask(centers, counts, excluded_windows)
    if np.count_nonzero(fit_mask) < 8:
        raise RuntimeError("Too few bins for spline background fit")

    def fit(mask):
        w = np.sqrt(np.maximum(counts[mask], 1.0))
        s = float(smoothing_factor) * float(np.count_nonzero(mask))
        spline = UnivariateSpline(x[mask], y[mask], w=w, k=3, s=s, ext=0)
        return np.exp(np.clip(spline(x), -30, 50))

    for _ in range(iterations):
        model = fit(fit_mask)
        residual = (counts - model) / np.sqrt(np.maximum(model, 1.0))
        new_mask = fit_mask & (np.abs(residual) <= clip_sigma)
        if np.array_equal(new_mask, fit_mask):
            break
        fit_mask = new_mask
        if np.count_nonzero(fit_mask) < 8:
            raise RuntimeError("Robust spline fit clipped too many bins")

    return fit(fit_mask), fit_mask


def fit_background_family(
    centers,
    counts,
    *,
    family="chebyshev",
    degree=7,
    iterations=6,
    clip_sigma=3.5,
    excluded_windows=(),
    spline_smoothing_factor=1.0,
):
    family = str(family).lower()
    if family == "chebyshev":
        return robust_smooth_background(
            centers,
            counts,
            degree=degree,
            iterations=iterations,
            clip_sigma=clip_sigma,
            excluded_windows=excluded_windows,
        )
    if family == "bernstein":
        return robust_bernstein_background(
            centers,
            counts,
            degree=degree,
            iterations=iterations,
            clip_sigma=clip_sigma,
            excluded_windows=excluded_windows,
        )
    if family == "spline":
        return robust_spline_background(
            centers,
            counts,
            iterations=iterations,
            clip_sigma=clip_sigma,
            excluded_windows=excluded_windows,
            smoothing_factor=spline_smoothing_factor,
        )
    raise ValueError(f"Unknown background family: {family}")
