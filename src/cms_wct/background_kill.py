from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Iterable

import numpy as np

from .background_cv import (
    BackgroundCandidate,
    analysis_eligibility,
    cross_validate_candidates,
    frozen_background_candidates,
)
from .background_families import fit_background_family
from .locked import fit_phase_locked_waveform
from .signature import weighted_linear_sinusoid


@dataclass(frozen=True)
class PairKillConfig:
    omega: float
    phase: float
    excluded_windows: tuple[tuple[float, float], ...] = ()
    n_folds: int = 5
    block_size: int = 16
    iterations: int = 6
    clip_sigma: float = 3.5
    min_model_count: float = 5.0

    def to_dict(self) -> dict:
        return {
            "omega": float(self.omega),
            "phase": float(self.phase),
            "excluded_windows": [[float(lo), float(hi)] for lo, hi in self.excluded_windows],
            "n_folds": int(self.n_folds),
            "block_size": int(self.block_size),
            "iterations": int(self.iterations),
            "clip_sigma": float(self.clip_sigma),
            "min_model_count": float(self.min_model_count),
            "primary_pair_statistic": "min(locked_delta_chi2_a, locked_delta_chi2_b)",
            "secondary_pair_statistic": "locked_delta_chi2_a + locked_delta_chi2_b",
        }


def _as_spectrum(centers: np.ndarray, counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    centers = np.asarray(centers, dtype=float)
    counts = np.asarray(counts, dtype=float)
    if centers.ndim != 1 or counts.ndim != 1:
        raise ValueError("centers and counts must be one-dimensional")
    if len(centers) != len(counts):
        raise ValueError("centers and counts lengths differ")
    if np.any(~np.isfinite(centers)) or np.any(centers <= 0.0):
        raise ValueError("centers must be finite and positive")
    if np.any(~np.isfinite(counts)) or np.any(counts < 0.0):
        raise ValueError("counts must be finite and non-negative")
    return centers, counts


def _fit_model(
    centers: np.ndarray,
    counts: np.ndarray,
    candidate: BackgroundCandidate,
    config: PairKillConfig,
) -> tuple[np.ndarray, np.ndarray]:
    model, fit_mask = fit_background_family(
        centers,
        counts,
        family=candidate.family,
        degree=int(candidate.degree),
        iterations=int(config.iterations),
        clip_sigma=float(config.clip_sigma),
        excluded_windows=config.excluded_windows,
        spline_smoothing_factor=float(candidate.spline_smoothing_factor),
    )
    model = np.asarray(model, dtype=float)
    if model.shape != counts.shape or np.any(~np.isfinite(model)) or np.any(model <= 0.0):
        raise RuntimeError(f"invalid model returned by {candidate.name}")
    return model, np.asarray(fit_mask, dtype=bool)


def evaluate_candidate(
    centers: np.ndarray,
    counts: np.ndarray,
    candidate: BackgroundCandidate,
    config: PairKillConfig,
) -> tuple[dict, np.ndarray]:
    """Detrend one spectrum and evaluate the frozen frequency and phase.

    The background fit never receives the WCT frequency or phase. They are used
    only after detrending, when the fixed-frequency/free-phase and fully locked
    statistics are evaluated.
    """
    centers, counts = _as_spectrum(centers, counts)
    model, fit_mask = _fit_model(centers, counts, candidate, config)
    residuals = (counts - model) / np.sqrt(np.maximum(model, 1.0))

    analysis_mask = analysis_eligibility(centers, config.excluded_windows)
    analysis_mask &= np.isfinite(residuals)
    analysis_mask &= model >= float(config.min_model_count)
    if np.count_nonzero(analysis_mask) < 8:
        raise RuntimeError("too few bins remain for frozen-waveform evaluation")

    x = np.log(centers[analysis_mask])
    y = residuals[analysis_mask]
    free = weighted_linear_sinusoid(x, y, float(config.omega))
    locked = fit_phase_locked_waveform(
        x,
        y,
        float(config.omega),
        float(config.phase),
        positive_only=True,
    )

    return {
        "candidate": candidate.to_dict(),
        "analysis_bins": int(np.count_nonzero(analysis_mask)),
        "background_fit_bins": int(np.count_nonzero(fit_mask)),
        "residual_rms": float(np.sqrt(np.mean(np.square(y)))),
        "residual_max_abs": float(np.max(np.abs(y))),
        "fixed_omega_amplitude": float(free.amplitude),
        "fixed_omega_phase": float(free.phase),
        "fixed_omega_delta_chi2": float(free.delta_chi2),
        "phase_offset_from_prediction": float(
            np.angle(np.exp(1j * (float(free.phase) - float(config.phase))))
        ),
        "locked_signed_amplitude": float(locked.signed_amplitude),
        "locked_delta_chi2": float(locked.delta_chi2),
        "locked_local_p_one_sided": float(locked.local_p_one_sided_chi2_1dof),
        "positive_sign": bool(locked.signed_amplitude > 0.0),
    }, model


def pair_cv_scores(
    centers_a: np.ndarray,
    counts_a: np.ndarray,
    centers_b: np.ndarray,
    counts_b: np.ndarray,
    candidates: Iterable[BackgroundCandidate] | None,
    config: PairKillConfig,
) -> list[dict]:
    """Rank backgrounds by combined WCT-blind held-out Poisson deviance."""
    centers_a, counts_a = _as_spectrum(centers_a, counts_a)
    centers_b, counts_b = _as_spectrum(centers_b, counts_b)
    if candidates is None:
        candidates = frozen_background_candidates()
    candidates = list(candidates)
    if not candidates:
        raise ValueError("at least one background candidate is required")

    cv_a = cross_validate_candidates(
        centers_a,
        counts_a,
        candidates,
        excluded_windows=config.excluded_windows,
        n_folds=int(config.n_folds),
        block_size=int(config.block_size),
        iterations=int(config.iterations),
        clip_sigma=float(config.clip_sigma),
    )
    cv_b = cross_validate_candidates(
        centers_b,
        counts_b,
        candidates,
        excluded_windows=config.excluded_windows,
        n_folds=int(config.n_folds),
        block_size=int(config.block_size),
        iterations=int(config.iterations),
        clip_sigma=float(config.clip_sigma),
    )
    by_a = {row["name"]: row for row in cv_a}
    by_b = {row["name"]: row for row in cv_b}

    rows = []
    for candidate in candidates:
        a = by_a[candidate.name]
        b = by_b[candidate.name]
        n_bins = int(a["n_validation_bins"] + b["n_validation_bins"])
        total = float(a["total_poisson_deviance"] + b["total_poisson_deviance"])
        failed = None
        if a["failed"] is not None or b["failed"] is not None:
            failed = {"a": a["failed"], "b": b["failed"]}
            total = float("inf")
        rows.append({
            **candidate.to_dict(),
            "a_failed": a["failed"],
            "b_failed": b["failed"],
            "failed": failed,
            "a_total_poisson_deviance": float(a["total_poisson_deviance"]),
            "b_total_poisson_deviance": float(b["total_poisson_deviance"]),
            "a_deviance_per_bin": float(a["deviance_per_bin"]),
            "b_deviance_per_bin": float(b["deviance_per_bin"]),
            "combined_validation_bins": n_bins,
            "combined_total_poisson_deviance": total,
            "combined_deviance_per_bin": total / n_bins if n_bins > 0 else float("inf"),
        })
    return rows


def select_pair_background(
    centers_a: np.ndarray,
    counts_a: np.ndarray,
    centers_b: np.ndarray,
    counts_b: np.ndarray,
    candidates: Iterable[BackgroundCandidate] | None,
    config: PairKillConfig,
) -> tuple[BackgroundCandidate, list[dict], list[str]]:
    candidates = frozen_background_candidates() if candidates is None else list(candidates)
    rows = pair_cv_scores(
        centers_a, counts_a, centers_b, counts_b, candidates, config
    )
    finite = [row for row in rows if np.isfinite(row["combined_total_poisson_deviance"])]
    if not finite:
        raise RuntimeError("all background candidates failed pair cross-validation")
    winner_row = min(finite, key=lambda row: row["combined_total_poisson_deviance"])
    winner = next(candidate for candidate in candidates if candidate.name == winner_row["name"])
    best = float(winner_row["combined_total_poisson_deviance"])
    tol = max(1e-9, 1e-10 * abs(best))
    ties = [
        row["name"]
        for row in finite
        if abs(float(row["combined_total_poisson_deviance"]) - best) <= tol
    ]
    return winner, rows, ties


def evaluate_pair_pipeline(
    centers_a: np.ndarray,
    counts_a: np.ndarray,
    centers_b: np.ndarray,
    counts_b: np.ndarray,
    candidates: Iterable[BackgroundCandidate] | None,
    config: PairKillConfig,
    *,
    reselect_background: bool = True,
    fixed_candidate: BackgroundCandidate | None = None,
) -> tuple[dict, np.ndarray, np.ndarray]:
    candidates = frozen_background_candidates() if candidates is None else list(candidates)
    if reselect_background:
        winner, cv_rows, ties = select_pair_background(
            centers_a, counts_a, centers_b, counts_b, candidates, config
        )
    else:
        if fixed_candidate is None:
            raise ValueError("fixed_candidate is required when reselect_background=False")
        winner = fixed_candidate
        cv_rows = []
        ties = [winner.name]

    eval_a, model_a = evaluate_candidate(centers_a, counts_a, winner, config)
    eval_b, model_b = evaluate_candidate(centers_b, counts_b, winner, config)
    primary = float(min(eval_a["locked_delta_chi2"], eval_b["locked_delta_chi2"]))
    secondary = float(eval_a["locked_delta_chi2"] + eval_b["locked_delta_chi2"])

    return {
        "selected_background": winner.to_dict(),
        "winner_ties_within_tolerance": ties,
        "cv_scores": cv_rows,
        "spectrum_a": eval_a,
        "spectrum_b": eval_b,
        "primary_pair_score_min_locked_delta_chi2": primary,
        "secondary_pair_score_sum_locked_delta_chi2": secondary,
        "positive_sign_in_both": bool(eval_a["positive_sign"] and eval_b["positive_sign"]),
    }, model_a, model_b


def _add_one_p(scores: np.ndarray, observed: float) -> tuple[int, float]:
    scores = np.asarray(scores, dtype=float)
    exceedances = int(np.count_nonzero(scores >= float(observed)))
    return exceedances, float((exceedances + 1.0) / (len(scores) + 1.0))


def spurious_signal_null_pair(
    centers_a: np.ndarray,
    counts_a: np.ndarray,
    centers_b: np.ndarray,
    counts_b: np.ndarray,
    candidates: Iterable[BackgroundCandidate] | None,
    config: PairKillConfig,
    *,
    n_trials: int = 1000,
    seed: int = 20260901,
    reselect_each_trial: bool = True,
) -> dict:
    """End-to-end detrending null.

    Smooth pseudo-spectra are generated from the WCT-blind CV winner. Every
    pseudoexperiment then repeats background selection (by default), detrending,
    and the frozen waveform test. This directly asks whether the analysis
    pipeline can manufacture the locked waveform under a smooth null.
    """
    if n_trials < 1:
        raise ValueError("n_trials must be positive")
    candidates = frozen_background_candidates() if candidates is None else list(candidates)
    observed, model_a, model_b = evaluate_pair_pipeline(
        centers_a, counts_a, centers_b, counts_b, candidates, config
    )
    generating_candidate = next(
        c for c in candidates if c.name == observed["selected_background"]["name"]
    )

    rng = np.random.default_rng(seed)
    primary_scores = np.empty(n_trials, dtype=float)
    secondary_scores = np.empty(n_trials, dtype=float)
    amp_a = np.empty(n_trials, dtype=float)
    amp_b = np.empty(n_trials, dtype=float)
    selected_names: list[str] = []

    for i in range(n_trials):
        pseudo_a = rng.poisson(model_a).astype(float)
        pseudo_b = rng.poisson(model_b).astype(float)
        result, _, _ = evaluate_pair_pipeline(
            centers_a,
            pseudo_a,
            centers_b,
            pseudo_b,
            candidates,
            config,
            reselect_background=bool(reselect_each_trial),
            fixed_candidate=generating_candidate,
        )
        primary_scores[i] = result["primary_pair_score_min_locked_delta_chi2"]
        secondary_scores[i] = result["secondary_pair_score_sum_locked_delta_chi2"]
        amp_a[i] = result["spectrum_a"]["locked_signed_amplitude"]
        amp_b[i] = result["spectrum_b"]["locked_signed_amplitude"]
        selected_names.append(result["selected_background"]["name"])

    primary_exceed, primary_p = _add_one_p(
        primary_scores, observed["primary_pair_score_min_locked_delta_chi2"]
    )
    secondary_exceed, secondary_p = _add_one_p(
        secondary_scores, observed["secondary_pair_score_sum_locked_delta_chi2"]
    )
    return {
        "classification": "retrospective_end_to_end_spurious_signal_null",
        "config": config.to_dict(),
        "n_trials": int(n_trials),
        "seed": int(seed),
        "reselect_background_each_trial": bool(reselect_each_trial),
        "generating_background": generating_candidate.to_dict(),
        "observed": observed,
        "primary_exceedances": primary_exceed,
        "primary_add_one_p": primary_p,
        "secondary_exceedances": secondary_exceed,
        "secondary_add_one_p": secondary_p,
        "primary_score_quantiles": {
            "q50": float(np.quantile(primary_scores, 0.50)),
            "q90": float(np.quantile(primary_scores, 0.90)),
            "q95": float(np.quantile(primary_scores, 0.95)),
            "q99": float(np.quantile(primary_scores, 0.99)),
        },
        "selected_background_counts": dict(Counter(selected_names)),
        "trials": [
            {
                "trial": i,
                "selected_background": selected_names[i],
                "primary_pair_score": float(primary_scores[i]),
                "secondary_pair_score": float(secondary_scores[i]),
                "locked_amplitude_a": float(amp_a[i]),
                "locked_amplitude_b": float(amp_b[i]),
            }
            for i in range(n_trials)
        ],
    }


def injected_mean(
    centers: np.ndarray,
    smooth_model: np.ndarray,
    config: PairKillConfig,
    amplitude: float,
) -> np.ndarray:
    """Inject a waveform whose amplitude is expressed in Pearson-residual units."""
    centers = np.asarray(centers, dtype=float)
    smooth_model = np.asarray(smooth_model, dtype=float)
    if centers.shape != smooth_model.shape:
        raise ValueError("centers and smooth_model shapes differ")
    mean = smooth_model.copy()
    mask = analysis_eligibility(centers, config.excluded_windows)
    mask &= smooth_model >= float(config.min_model_count)
    basis = np.cos(float(config.omega) * np.log(centers[mask]) - float(config.phase))
    mean[mask] += float(amplitude) * np.sqrt(np.maximum(smooth_model[mask], 1.0)) * basis
    return np.maximum(mean, 1e-9)


def deterministic_absorption_matrix(
    centers_a: np.ndarray,
    counts_a: np.ndarray,
    centers_b: np.ndarray,
    counts_b: np.ndarray,
    candidates: Iterable[BackgroundCandidate] | None,
    config: PairKillConfig,
    *,
    injection_amplitudes: Iterable[float] = (0.25, 0.5, 0.75, 1.0),
) -> dict:
    """Measure noiseless signal absorption by every prespecified background."""
    candidates = frozen_background_candidates() if candidates is None else list(candidates)
    observed, base_a, base_b = evaluate_pair_pipeline(
        centers_a, counts_a, centers_b, counts_b, candidates, config
    )
    rows = []
    for amplitude in injection_amplitudes:
        amplitude = float(amplitude)
        if amplitude <= 0.0:
            raise ValueError("injection amplitudes must be positive")
        mean_a = injected_mean(centers_a, base_a, config, amplitude)
        mean_b = injected_mean(centers_b, base_b, config, amplitude)
        for candidate in candidates:
            eva, _ = evaluate_candidate(centers_a, mean_a, candidate, config)
            evb, _ = evaluate_candidate(centers_b, mean_b, candidate, config)
            rows.append({
                "injected_amplitude": amplitude,
                "analysis_candidate": candidate.name,
                "recovered_amplitude_a": float(eva["locked_signed_amplitude"]),
                "recovered_amplitude_b": float(evb["locked_signed_amplitude"]),
                "retention_a": float(eva["locked_signed_amplitude"] / amplitude),
                "retention_b": float(evb["locked_signed_amplitude"] / amplitude),
                "min_retention": float(
                    min(eva["locked_signed_amplitude"], evb["locked_signed_amplitude"]) / amplitude
                ),
                "primary_pair_score": float(
                    min(eva["locked_delta_chi2"], evb["locked_delta_chi2"])
                ),
            })
    return {
        "classification": "deterministic_signal_absorption_matrix",
        "generating_background": observed["selected_background"],
        "rows": rows,
    }


def injection_recovery_pair(
    centers_a: np.ndarray,
    counts_a: np.ndarray,
    centers_b: np.ndarray,
    counts_b: np.ndarray,
    candidates: Iterable[BackgroundCandidate] | None,
    config: PairKillConfig,
    *,
    injection_amplitudes: Iterable[float] = (0.25, 0.5, 0.75, 1.0),
    n_trials_per_amplitude: int = 200,
    seed: int = 20260902,
    null_primary_threshold: float | None = None,
    reselect_each_trial: bool = True,
) -> dict:
    """Poisson injection/recovery with the full background-selection pipeline."""
    if n_trials_per_amplitude < 1:
        raise ValueError("n_trials_per_amplitude must be positive")
    candidates = frozen_background_candidates() if candidates is None else list(candidates)
    observed, base_a, base_b = evaluate_pair_pipeline(
        centers_a, counts_a, centers_b, counts_b, candidates, config
    )
    generating_candidate = next(
        c for c in candidates if c.name == observed["selected_background"]["name"]
    )
    rng = np.random.default_rng(seed)
    rows = []
    summaries = []

    for amplitude in injection_amplitudes:
        amplitude = float(amplitude)
        if amplitude <= 0.0:
            raise ValueError("injection amplitudes must be positive")
        mean_a = injected_mean(centers_a, base_a, config, amplitude)
        mean_b = injected_mean(centers_b, base_b, config, amplitude)
        recovered_a = []
        recovered_b = []
        pair_scores = []
        selected_names = []

        for trial in range(n_trials_per_amplitude):
            pseudo_a = rng.poisson(mean_a).astype(float)
            pseudo_b = rng.poisson(mean_b).astype(float)
            result, _, _ = evaluate_pair_pipeline(
                centers_a,
                pseudo_a,
                centers_b,
                pseudo_b,
                candidates,
                config,
                reselect_background=bool(reselect_each_trial),
                fixed_candidate=generating_candidate,
            )
            ra = float(result["spectrum_a"]["locked_signed_amplitude"])
            rb = float(result["spectrum_b"]["locked_signed_amplitude"])
            score = float(result["primary_pair_score_min_locked_delta_chi2"])
            selected = result["selected_background"]["name"]
            recovered_a.append(ra)
            recovered_b.append(rb)
            pair_scores.append(score)
            selected_names.append(selected)
            rows.append({
                "injected_amplitude": amplitude,
                "trial": trial,
                "selected_background": selected,
                "recovered_amplitude_a": ra,
                "recovered_amplitude_b": rb,
                "primary_pair_score": score,
                "positive_sign_in_both": bool(result["positive_sign_in_both"]),
            })

        ra_arr = np.asarray(recovered_a)
        rb_arr = np.asarray(recovered_b)
        score_arr = np.asarray(pair_scores)
        min_amp = np.minimum(ra_arr, rb_arr)
        summary = {
            "injected_amplitude": amplitude,
            "median_recovered_amplitude_a": float(np.median(ra_arr)),
            "median_recovered_amplitude_b": float(np.median(rb_arr)),
            "median_min_recovered_amplitude": float(np.median(min_amp)),
            "median_min_retention": float(np.median(min_amp) / amplitude),
            "positive_sign_in_both_fraction": float(np.mean((ra_arr > 0.0) & (rb_arr > 0.0))),
            "median_primary_pair_score": float(np.median(score_arr)),
            "selected_background_counts": dict(Counter(selected_names)),
        }
        if null_primary_threshold is not None:
            summary["power_above_null_threshold"] = float(
                np.mean(score_arr >= float(null_primary_threshold))
            )
            summary["null_primary_threshold"] = float(null_primary_threshold)
        summaries.append(summary)

    return {
        "classification": "end_to_end_signal_injection_recovery",
        "config": config.to_dict(),
        "seed": int(seed),
        "n_trials_per_amplitude": int(n_trials_per_amplitude),
        "reselect_background_each_trial": bool(reselect_each_trial),
        "generating_background": observed["selected_background"],
        "summaries": summaries,
        "trials": rows,
    }
