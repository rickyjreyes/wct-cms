from __future__ import annotations

import numpy as np
from numpy.polynomial import Chebyshev


def make_histogram(masses, mass_min, mass_max, bins, log_bins):
    selected = masses[(masses >= mass_min) & (masses <= mass_max)]
    if log_bins:
        if mass_min <= 0:
            raise ValueError("mass_min must be > 0 with logarithmic bins")
        edges = np.geomspace(mass_min, mass_max, bins + 1)
        centers = np.sqrt(edges[:-1] * edges[1:])
    else:
        edges = np.linspace(mass_min, mass_max, bins + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])
    counts, _ = np.histogram(selected, bins=edges)
    return counts.astype(float), edges, centers


def robust_smooth_background(
    centers,
    counts,
    degree=7,
    iterations=6,
    clip_sigma=3.5,
    excluded_windows=(),
):
    x = np.log(centers)
    y = np.log(counts + 0.5)
    fit_mask = counts > 0
    for lo, hi in excluded_windows:
        fit_mask &= ~((centers >= lo) & (centers <= hi))

    minimum = degree + 3
    if np.count_nonzero(fit_mask) < minimum:
        raise RuntimeError("Too few bins for background fit")

    for _ in range(iterations):
        w = np.sqrt(np.maximum(counts[fit_mask], 1.0))
        poly = Chebyshev.fit(
            x[fit_mask], y[fit_mask], deg=degree, w=w, domain=[x.min(), x.max()]
        )
        model = np.exp(np.clip(poly(x), -30, 50))
        residual = (counts - model) / np.sqrt(np.maximum(model, 1.0))
        new_mask = fit_mask & (np.abs(residual) <= clip_sigma)
        if np.array_equal(new_mask, fit_mask):
            break
        fit_mask = new_mask
        if np.count_nonzero(fit_mask) < minimum:
            raise RuntimeError("Robust fit clipped too many bins")

    w = np.sqrt(np.maximum(counts[fit_mask], 1.0))
    poly = Chebyshev.fit(
        x[fit_mask], y[fit_mask], deg=degree, w=w, domain=[x.min(), x.max()]
    )
    model = np.exp(np.clip(poly(x), -30, 50))
    return model, fit_mask
