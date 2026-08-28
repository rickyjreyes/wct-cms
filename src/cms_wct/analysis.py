from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from .background import make_histogram, robust_smooth_background
from .cmsio import expand_inputs, extract_dimuon_masses, load_golden_json
from .models import AnalysisSummary
from .plots import write_plots
from .signature import permutation_null, scan_omegas, weighted_linear_sinusoid


def parse_windows(raw):
    out = []
    for item in raw:
        lo, hi = map(float, item.split(":", 1))
        if lo >= hi:
            raise ValueError(f"Invalid mask window: {item}")
        out.append((lo, hi))
    return out


def run_analysis(args):
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    inputs = expand_inputs(args.input)
    golden = load_golden_json(args.golden_json)
    windows = parse_windows(args.mask_window)

    (outdir / "config.json").write_text(json.dumps(vars(args), indent=2, sort_keys=True) + "\n")
    (outdir / "inputs.txt").write_text("\n".join(inputs) + "\n")

    masses, counters = extract_dimuon_masses(inputs, args, golden)
    counts, edges, centers = make_histogram(
        masses, args.mass_min, args.mass_max, args.bins, args.log_bins
    )
    model, bg_mask = robust_smooth_background(
        centers, counts, args.fit_degree, args.fit_iterations, args.clip_sigma, windows
    )
    residuals = (counts - model) / np.sqrt(np.maximum(model, 1.0))
    mask = np.isfinite(residuals) & (model >= args.min_model_count) & (centers > 0)
    if np.count_nonzero(mask) < 20:
        raise RuntimeError("Too few usable mass bins")

    x = np.log(centers[mask] / args.m0)
    y = residuals[mask] - np.mean(residuals[mask])
    omegas = np.linspace(args.omega_min, args.omega_max, args.omega_steps)
    scores, best = scan_omegas(x, y, omegas)
    frozen = None if args.frozen_omega is None else weighted_linear_sinusoid(x, y, args.frozen_omega)

    p_global, p_frozen, null_max, null_frozen = permutation_null(
        x, y, omegas, best.delta_chi2, args.frozen_omega,
        None if frozen is None else frozen.delta_chi2,
        args.permutations, args.seed,
    )

    np.savetxt(outdir / "omega_scan.csv", np.column_stack([omegas, scores]), delimiter=",", header="omega,delta_chi2", comments="")
    np.savetxt(outdir / "spectrum.csv", np.column_stack([
        edges[:-1], edges[1:], centers, counts, model, residuals, mask.astype(int), bg_mask.astype(int)
    ]), delimiter=",", header="mass_low,mass_high,mass_center,count,background,residual,analysis_mask,background_fit_mask", comments="")
    if len(null_max):
        np.savetxt(outdir / "permutation_global_max.csv", null_max, delimiter=",", header="max_delta_chi2", comments="")
    if null_frozen is not None:
        np.savetxt(outdir / "permutation_frozen.csv", null_frozen, delimiter=",", header="frozen_delta_chi2", comments="")

    summary = AnalysisSummary(
        events_read=counters["events_read"],
        events_after_json=counters["events_after_json"],
        selected_muons=counters["selected_muons"],
        opposite_sign_pairs=counters["opposite_sign_pairs"],
        masses_in_range=int(np.sum(counts)),
        best_scan=asdict(best),
        frozen_scan=None if frozen is None else asdict(frozen),
        global_permutation_p=p_global,
        frozen_permutation_p=p_frozen,
        mass_min=args.mass_min,
        mass_max=args.mass_max,
        bins=args.bins,
        log_bins=args.log_bins,
        fit_degree=args.fit_degree,
        omega_min=args.omega_min,
        omega_max=args.omega_max,
        omega_steps=args.omega_steps,
        permutations=args.permutations,
        seed=args.seed,
    )
    (outdir / "summary.json").write_text(json.dumps(asdict(summary), indent=2, sort_keys=True) + "\n")

    write_plots(outdir, centers, counts, model, residuals, mask, omegas, scores, best, frozen, x, y, null_max, args.log_bins)
    return summary
