from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from cms_wct.background import make_histogram
from cms_wct.background_families import fit_background_family
from cms_wct.cmsio import extract_dimuon_masses, load_golden_json
from cms_wct.locked import fit_phase_locked_waveform, phase_locked_permutation_null


def parse_windows(raw):
    out = []
    for item in raw:
        lo, hi = map(float, item.split(":", 1))
        if lo >= hi:
            raise ValueError(f"Invalid mask window: {item}")
        out.append((lo, hi))
    return out


def parser():
    p = argparse.ArgumentParser(description="Prospective CMS fixed-omega fixed-phase signed-amplitude test")
    p.add_argument("--input", nargs="+", required=True)
    p.add_argument("--tree", default="Events")
    p.add_argument("--golden-json", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--omega", type=float, required=True)
    p.add_argument("--phase", type=float, required=True)
    p.add_argument("--mass-min", type=float, default=2.0)
    p.add_argument("--mass-max", type=float, default=120.0)
    p.add_argument("--bins", type=int, default=350)
    p.add_argument("--log-bins", action="store_true")
    p.add_argument("--muon-pt-min", type=float, default=4.0)
    p.add_argument("--muon-eta-max", type=float, default=2.4)
    p.add_argument("--tight-id", action="store_true")
    p.add_argument("--medium-id", action="store_true")
    p.add_argument("--max-events", type=int, default=100000)
    p.add_argument("--step-size", default="100 MB")
    p.add_argument("--background-family", choices=["chebyshev", "bernstein", "spline"], default="chebyshev")
    p.add_argument("--fit-degree", type=int, default=7)
    p.add_argument("--fit-iterations", type=int, default=6)
    p.add_argument("--clip-sigma", type=float, default=3.5)
    p.add_argument("--spline-smoothing-factor", type=float, default=1.0)
    p.add_argument("--min-model-count", type=float, default=5.0)
    p.add_argument("--mask-window", action="append", default=[])
    p.add_argument("--m0", type=float, default=1.0)
    p.add_argument("--permutations", type=int, default=1000)
    p.add_argument("--parametric-bootstrap", type=int, default=500)
    p.add_argument("--seed", type=int, default=20260831)
    return p


def main():
    args = parser().parse_args()
    if args.tight_id and args.medium_id:
        raise SystemExit("Choose only one of --tight-id or --medium-id")
    if not args.tight_id and not args.medium_id:
        raise SystemExit("Explicitly select --tight-id or --medium-id")

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    windows = parse_windows(args.mask_window)
    golden = load_golden_json(args.golden_json)

    masses, counters = extract_dimuon_masses(args.input, args, golden)
    counts, edges, centers = make_histogram(masses, args.mass_min, args.mass_max, args.bins, args.log_bins)
    model, bg_mask = fit_background_family(
        centers,
        counts,
        family=args.background_family,
        degree=args.fit_degree,
        iterations=args.fit_iterations,
        clip_sigma=args.clip_sigma,
        excluded_windows=windows,
        spline_smoothing_factor=args.spline_smoothing_factor,
    )
    residuals = (counts - model) / np.sqrt(np.maximum(model, 1.0))
    mask = np.isfinite(residuals) & (model >= args.min_model_count) & (centers > 0)
    for lo, hi in windows:
        mask &= ~((centers >= lo) & (centers <= hi))
    if np.count_nonzero(mask) < 20:
        raise RuntimeError("Too few usable mass bins")

    x = np.log(centers[mask] / args.m0)
    raw_y = residuals[mask]
    y = raw_y - np.mean(raw_y)
    observed = fit_phase_locked_waveform(x, y, args.omega, args.phase, positive_only=True)
    p_perm, perm_scores = phase_locked_permutation_null(
        x, y, args.omega, args.phase, observed.delta_chi2, args.permutations, args.seed, positive_only=True
    )

    rng = np.random.default_rng(args.seed + 1)
    bootstrap_scores = np.empty(args.parametric_bootstrap, dtype=float)
    for i in range(args.parametric_bootstrap):
        pseudo_counts = rng.poisson(model).astype(float)
        pseudo_model, _ = fit_background_family(
            centers,
            pseudo_counts,
            family=args.background_family,
            degree=args.fit_degree,
            iterations=args.fit_iterations,
            clip_sigma=args.clip_sigma,
            excluded_windows=windows,
            spline_smoothing_factor=args.spline_smoothing_factor,
        )
        pseudo_residuals = (pseudo_counts - pseudo_model) / np.sqrt(np.maximum(pseudo_model, 1.0))
        py = pseudo_residuals[mask]
        py = py - np.mean(py)
        bootstrap_scores[i] = fit_phase_locked_waveform(
            x, py, args.omega, args.phase, positive_only=True
        ).delta_chi2

    p_boot = None
    if args.parametric_bootstrap > 0:
        p_boot = float((np.count_nonzero(bootstrap_scores >= observed.delta_chi2) + 1.0) / (args.parametric_bootstrap + 1.0))

    summary = {
        "test": "fixed_omega_fixed_phase_positive_amplitude",
        "omega": float(args.omega),
        "phase": float(args.phase),
        "positive_amplitude_required": True,
        "observed": observed.to_dict(),
        "permutation_p": p_perm,
        "parametric_bootstrap_p": p_boot,
        "permutations": int(args.permutations),
        "parametric_bootstrap": int(args.parametric_bootstrap),
        "background_family": args.background_family,
        "fit_degree": int(args.fit_degree),
        "spline_smoothing_factor": float(args.spline_smoothing_factor),
        "residual_rms": float(np.sqrt(np.mean(np.square(raw_y)))),
        "residual_max_abs": float(np.max(np.abs(raw_y))),
        "analysis_bins": int(np.count_nonzero(mask)),
        "events_read": int(counters["events_read"]),
        "events_after_json": int(counters["events_after_json"]),
        "selected_muons": int(counters["selected_muons"]),
        "opposite_sign_pairs": int(counters["opposite_sign_pairs"]),
        "masses_in_range": int(np.sum(counts)),
        "config": vars(args),
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    np.savetxt(outdir / "permutation_locked.csv", perm_scores, delimiter=",", header="delta_chi2", comments="")
    if args.parametric_bootstrap > 0:
        np.savetxt(outdir / "bootstrap_locked.csv", bootstrap_scores, delimiter=",", header="delta_chi2", comments="")
    np.savetxt(
        outdir / "spectrum.csv",
        np.column_stack([edges[:-1], edges[1:], centers, counts, model, residuals, mask.astype(int), bg_mask.astype(int)]),
        delimiter=",",
        header="mass_low,mass_high,mass_center,count,background,residual,analysis_mask,background_fit_mask",
        comments="",
    )

    print("Locked omega:", observed.omega)
    print("Locked phase:", observed.phase)
    print("Signed amplitude:", observed.signed_amplitude)
    print("Locked delta chi-square:", observed.delta_chi2)
    print("Local one-sided p (1 dof boundary):", observed.local_p_one_sided_chi2_1dof)
    print("Locked permutation p:", p_perm)
    print("Locked refit-bootstrap p:", p_boot)
    print("Residual RMS:", summary["residual_rms"])
    print("Maximum |residual|:", summary["residual_max_abs"])
    print("Summary:", outdir / "summary.json")


if __name__ == "__main__":
    main()
