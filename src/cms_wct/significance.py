from __future__ import annotations

import math

from scipy.stats import beta, norm


def one_sided_gaussian_p(z: float) -> float:
    """Return the one-sided Gaussian tail probability for a Z value."""
    z = float(z)
    if not math.isfinite(z):
        raise ValueError("z must be finite")
    return float(norm.sf(z))


def add_one_monte_carlo_p(exceedances: int, trials: int) -> float:
    """Finite-sample Monte Carlo p-value using the conservative add-one rule."""
    exceedances = int(exceedances)
    trials = int(trials)
    if trials < 0:
        raise ValueError("trials must be non-negative")
    if exceedances < 0 or exceedances > trials:
        raise ValueError("exceedances must satisfy 0 <= exceedances <= trials")
    return float((exceedances + 1.0) / (trials + 1.0))


def minimum_trials_for_add_one_resolution(target_p: float) -> int:
    """Smallest N for which 1/(N+1) is at or below target_p."""
    target_p = float(target_p)
    if not (0.0 < target_p < 1.0):
        raise ValueError("target_p must be between zero and one")
    return int(math.ceil(1.0 / target_p - 1.0))


def clopper_pearson_upper_bound(
    exceedances: int,
    trials: int,
    confidence: float = 0.95,
) -> float:
    """Exact one-sided upper confidence bound for a binomial exceedance rate."""
    exceedances = int(exceedances)
    trials = int(trials)
    confidence = float(confidence)
    if trials <= 0:
        raise ValueError("trials must be positive")
    if exceedances < 0 or exceedances > trials:
        raise ValueError("exceedances must satisfy 0 <= exceedances <= trials")
    if not (0.0 < confidence < 1.0):
        raise ValueError("confidence must be between zero and one")
    if exceedances == trials:
        return 1.0
    return float(beta.ppf(confidence, exceedances + 1, trials - exceedances))


def minimum_zero_exceedance_trials(
    target_p: float,
    confidence: float = 0.95,
) -> int:
    """Trials needed so zero exceedances imply an exact upper bound <= target_p."""
    target_p = float(target_p)
    confidence = float(confidence)
    if not (0.0 < target_p < 1.0):
        raise ValueError("target_p must be between zero and one")
    if not (0.0 < confidence < 1.0):
        raise ValueError("confidence must be between zero and one")
    alpha = 1.0 - confidence
    return int(math.ceil(math.log(alpha) / math.log1p(-target_p)))


def five_sigma_requirements() -> dict:
    """Canonical direct-Monte-Carlo requirements for a one-sided 5-sigma tail."""
    p5 = one_sided_gaussian_p(5.0)
    return {
        "z": 5.0,
        "one_sided_p": p5,
        "minimum_trials_add_one_resolution": minimum_trials_for_add_one_resolution(p5),
        "minimum_zero_exceedance_trials_95pct": minimum_zero_exceedance_trials(p5, 0.95),
        "minimum_zero_exceedance_trials_99pct": minimum_zero_exceedance_trials(p5, 0.99),
    }
