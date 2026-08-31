from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from cms_wct.background import make_histogram
from cms_wct.background_cv import cross_validate_candidates, frozen_background_candidates
from cms_wct.background_families import fit_background_family
from cms_wct.cmsio import expand_inputs, extract_dimuon_masses, load_golden_json
from cms_wct.locked import fit_phase_locked_waveform
from cms_wct.signature import weighted_linear_sinusoid

OMEGA = 7.025825825825827
PHASE_PRED = -0.2313916852932179
MASKS = [(2.9, 3.3), (3.55, 3.85), (8.5, 11.5), (80.0, 100.0)]


def make_args(manifest: str, max_events: int):
    return SimpleNamespace(
        input=[manifest],
        tree="Events",
        muon_pt_min=4.0,
        muon_eta_max=2.4,
        tight_id=True,
        medium_id=False,
        max_events=max_events,
        step_size="100 MB",
        mass_min=2.0,
        mass_max=120.0,
        bins=350,
        log_bins=True,
        min_model_count=5.0,
        fit_degree=7,
        fit_iterations=6,
        clip_sigma=3.5,
        m0=1.0,
    )


def load_histogram(label: str, manifest: str, golden, max_events: int):
    args = make_args(manifest, max_events)
    inputs = expand_inputs(args.input)
    print(f"Extracting {label} baseline event selection")
    masses, counters = extract_dimuon_masses(inputs, args, golden)
    counts, edges, centers = make_histogram(masses, 2.0, 120.0, 350, True)
    return {
        "label": label,
        "args": args,
        "counts": counts,
        "edges": edges,
        "centers": centers,
        "counters": counters,
    }


def winner_waveform(hist: dict, candidate: dict):
    centers = hist["centers"]
    counts = hist["counts"]
    model, fit_mask = fit_background_family(
        centers,
        counts,
        family=candidate["family"],
        degree=int(candidate["degree"]),
        iterations=6,
        clip_sigma=3.5,
        excluded_windows=MASKS,
        spline_smoothing_factor=float(candidate["spline_smoothing_factor"]),
    )
    residuals = (counts - model) / np.sqrt(np.maximum(model, 1.0))
    analysis_mask = np.isfinite(residuals) & (model >= 5.0) & (centers > 0.0)
    for lo, hi in MASKS:
        analysis_mask &= ~((centers >= lo) & (centers <= hi))

    x = np.log(centers[analysis_mask])
    raw_y = residuals[analysis_mask]
    y = raw_y - np.mean(raw_y)
    free = weighted_linear_sinusoid(x, y, OMEGA)
    locked = fit_phase_locked_waveform(x, y, OMEGA, PHASE_PRED, positive_only=True)
    phase_offset = float(np.angle(np.exp(1j * (free.phase - PHASE_PRED))))
    return {
        "analysis_bins": int(np.count_nonzero(analysis_mask)),
        "background_fit_bins": int(np.count_nonzero(fit_mask)),
        "residual_rms": float(np.sqrt(np.mean(np.square(raw_y)))),
        "residual_max_abs": float(np.max(np.abs(raw_y))),
        "fixed_omega_amplitude": float(free.amplitude),
        "fixed_omega_phase": float(free.phase),
        "phase_offset_from_prediction": phase_offset,
        "fixed_omega_delta_chi2": float(free.delta_chi2),
        "locked_signed_amplitude": float(locked.signed_amplitude),
        "locked_delta_chi2": float(locked.delta_chi2),
        "locked_local_p_one_sided": float(locked.local_p_one_sided_chi2_1dof),
        "positive_sign": bool(locked.signed_amplitude > 0.0),
    }


def main():
    p = argparse.ArgumentParser(
        description="Select the H/G continuum by WCT-blind blocked Poisson cross-validation, then evaluate the frozen waveform"
    )
    p.add_argument("--h-manifest", default="data/files_replication.txt")
    p.add_argument("--g-manifest", default="data/files_run2016g.txt")
    p.add_argument(
        "--golden-json",
        default="data/cert/14221/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt",
    )
    p.add_argument("--max-events", type=int, default=100000)
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--block-size", type=int, default=16)
    p.add_argument("--output-dir", default="results/hg_background_cv")
    args = p.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    candidates = frozen_background_candidates()
    freeze = {
        "classification": "retrospective_wct_blind_background_selection",
        "selection_statistic": "combined held-out Poisson deviance",
        "n_folds": args.folds,
        "block_size": args.block_size,
        "candidate_order": [c.to_dict() for c in candidates],
        "mass_range_GeV": [2.0, 120.0],
        "bins": 350,
        "log_bins": True,
        "masks_GeV": MASKS,
        "muon_pt_min_GeV": 4.0,
        "muon_eta_max": 2.4,
        "muon_id": "tight",
        "max_events": args.max_events,
        "wct_frequency_not_used_for_selection": True,
        "omega_evaluated_only_after_selection": OMEGA,
        "phase_evaluated_only_after_selection": PHASE_PRED,
    }
    (outdir / "selection_freeze.json").write_text(json.dumps(freeze, indent=2, sort_keys=True) + "\n")

    golden = load_golden_json(args.golden_json)
    h = load_histogram("Run2016H_file2", args.h_manifest, golden, args.max_events)
    g = load_histogram("Run2016G_file1", args.g_manifest, golden, args.max_events)

    print("Cross-validating continuum candidates on H (WCT-blind)")
    cv_h = cross_validate_candidates(
        h["centers"], h["counts"], candidates,
        excluded_windows=MASKS, n_folds=args.folds, block_size=args.block_size,
    )
    print("Cross-validating continuum candidates on G (WCT-blind)")
    cv_g = cross_validate_candidates(
        g["centers"], g["counts"], candidates,
        excluded_windows=MASKS, n_folds=args.folds, block_size=args.block_size,
    )

    by_h = {r["name"]: r for r in cv_h}
    by_g = {r["name"]: r for r in cv_g}
    rows = []
    for candidate in candidates:
        rh = by_h[candidate.name]
        rg = by_g[candidate.name]
        bins_total = rh["n_validation_bins"] + rg["n_validation_bins"]
        dev_total = rh["total_poisson_deviance"] + rg["total_poisson_deviance"]
        rows.append({
            **candidate.to_dict(),
            "h_failed": rh["failed"],
            "g_failed": rg["failed"],
            "h_validation_bins": rh["n_validation_bins"],
            "g_validation_bins": rg["n_validation_bins"],
            "h_total_poisson_deviance": rh["total_poisson_deviance"],
            "g_total_poisson_deviance": rg["total_poisson_deviance"],
            "h_deviance_per_bin": rh["deviance_per_bin"],
            "g_deviance_per_bin": rg["deviance_per_bin"],
            "combined_validation_bins": bins_total,
            "combined_total_poisson_deviance": dev_total,
            "combined_deviance_per_bin": dev_total / bins_total if bins_total else float("inf"),
        })

    finite = [r for r in rows if np.isfinite(r["combined_total_poisson_deviance"])]
    if not finite:
        raise RuntimeError("All background candidates failed cross-validation")
    winner = min(finite, key=lambda r: r["combined_total_poisson_deviance"])
    best_dev = float(winner["combined_total_poisson_deviance"])
    tie_tol = max(1e-9, 1e-10 * abs(best_dev))
    ties = [r["name"] for r in finite if abs(r["combined_total_poisson_deviance"] - best_dev) <= tie_tol]

    score_path = outdir / "cv_scores.csv"
    with score_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    selected = {
        "name": winner["name"],
        "family": winner["family"],
        "degree": int(winner["degree"]),
        "spline_smoothing_factor": float(winner["spline_smoothing_factor"]),
    }

    # Only now, after WCT-blind selection is complete, evaluate the frozen waveform.
    h_wave = winner_waveform(h, selected)
    g_wave = winner_waveform(g, selected)

    ranking = sorted(finite, key=lambda r: r["combined_total_poisson_deviance"])
    summary = {
        "classification": "retrospective_wct_blind_background_selection",
        "selection_used_wct_frequency_or_phase": False,
        "selection_statistic": "combined held-out Poisson deviance",
        "n_folds": args.folds,
        "block_size": args.block_size,
        "winner": selected,
        "winner_ties_within_tolerance": ties,
        "winner_combined_total_poisson_deviance": best_dev,
        "winner_combined_deviance_per_bin": float(winner["combined_deviance_per_bin"]),
        "top5_cv": [
            {
                "name": r["name"],
                "family": r["family"],
                "degree": int(r["degree"]),
                "spline_smoothing_factor": float(r["spline_smoothing_factor"]),
                "combined_total_poisson_deviance": float(r["combined_total_poisson_deviance"]),
                "combined_deviance_per_bin": float(r["combined_deviance_per_bin"]),
                "h_deviance_per_bin": float(r["h_deviance_per_bin"]),
                "g_deviance_per_bin": float(r["g_deviance_per_bin"]),
            }
            for r in ranking[:5]
        ],
        "post_selection_frozen_waveform": {
            "omega": OMEGA,
            "phase_prediction": PHASE_PRED,
            "Run2016H_file2": h_wave,
            "Run2016G_file1": g_wave,
            "positive_sign_in_both": bool(h_wave["positive_sign"] and g_wave["positive_sign"]),
        },
        "warning": "H/G are already-observed data. Background choice is WCT-blind, but this is a retrospective model-selection diagnostic, not a new prospective replication.",
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    print("CV winner:", selected)
    print("Winner ties:", ties)
    print("Combined held-out deviance/bin:", summary["winner_combined_deviance_per_bin"])
    print("H locked A, DeltaChi2:", h_wave["locked_signed_amplitude"], h_wave["locked_delta_chi2"])
    print("G locked A, DeltaChi2:", g_wave["locked_signed_amplitude"], g_wave["locked_delta_chi2"])
    print("Positive sign in both:", summary["post_selection_frozen_waveform"]["positive_sign_in_both"])
    print("Summary:", outdir / "summary.json")


if __name__ == "__main__":
    main()
