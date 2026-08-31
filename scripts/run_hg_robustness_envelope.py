from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from cms_wct.background import make_histogram
from cms_wct.background_families import fit_background_family
from cms_wct.cmsio import expand_inputs, extract_dimuon_masses, load_golden_json
from cms_wct.locked import fit_phase_locked_waveform
from cms_wct.signature import weighted_linear_sinusoid

OMEGA = 7.025825825825827
PHASE_PRED = -0.2313916853457179
BASE_MASKS = [(2.9, 3.3), (3.55, 3.85), (8.5, 11.5), (80.0, 100.0)]


def suite():
    out = [{"name": "baseline", "family": "chebyshev", "degree": 7}]
    out += [{"name": f"cheb_d{d}", "family": "chebyshev", "degree": d} for d in range(5, 13) if d != 7]
    out += [{"name": f"bernstein_d{d}", "family": "bernstein", "degree": d} for d in (5, 7, 9, 12)]
    out += [
        {"name": f"spline_s{s:g}", "family": "spline", "degree": 7, "spline": s}
        for s in (0.5, 1.0, 2.0)
    ]
    out += [
        {"name": "bins_300", "bins": 300},
        {"name": "bins_400", "bins": 400},
        {"name": "pt_5", "pt": 5.0},
        {"name": "pt_6", "pt": 6.0},
        {"name": "eta_2p1", "eta": 2.1},
        {"name": "medium_id", "id": "medium"},
        {"name": "mass_2_to_70", "mass_max": 70.0, "masks": BASE_MASKS[:3]},
        {"name": "mass_3_to_120", "mass_min": 3.0},
        {
            "name": "masks_narrow",
            "masks": [(2.95, 3.25), (3.60, 3.80), (8.8, 11.2), (82.0, 98.0)],
        },
        {
            "name": "masks_wide",
            "masks": [(2.8, 3.4), (3.45, 3.95), (8.0, 12.0), (78.0, 102.0)],
        },
    ]
    return out


def build_args(config, manifest, max_events):
    return SimpleNamespace(
        input=[manifest],
        tree="Events",
        golden_json=None,
        output_dir=None,
        muon_pt_min=float(config.get("pt", 4.0)),
        muon_eta_max=float(config.get("eta", 2.4)),
        tight_id=config.get("id", "tight") == "tight",
        medium_id=config.get("id", "tight") == "medium",
        max_events=max_events,
        step_size="100 MB",
        mass_min=float(config.get("mass_min", 2.0)),
        mass_max=float(config.get("mass_max", 120.0)),
        bins=int(config.get("bins", 350)),
        log_bins=True,
        min_model_count=5.0,
        fit_degree=int(config.get("degree", 7)),
        fit_iterations=6,
        clip_sigma=3.5,
        m0=1.0,
    )


def evaluate(label, manifest, golden_path, config, max_events):
    args = build_args(config, manifest, max_events)
    golden = load_golden_json(golden_path)
    inputs = expand_inputs(args.input)
    masses, counters = extract_dimuon_masses(inputs, args, golden)
    counts, _, centers = make_histogram(masses, args.mass_min, args.mass_max, args.bins, True)
    masks = config.get("masks", BASE_MASKS)
    family = config.get("family", "chebyshev")
    degree = int(config.get("degree", 7))
    spline = float(config.get("spline", 1.0))
    model, _ = fit_background_family(
        centers,
        counts,
        family=family,
        degree=degree,
        iterations=6,
        clip_sigma=3.5,
        excluded_windows=masks,
        spline_smoothing_factor=spline,
    )
    residuals = (counts - model) / np.sqrt(np.maximum(model, 1.0))
    mask = np.isfinite(residuals) & (model >= 5.0) & (centers > 0)
    for lo, hi in masks:
        mask &= ~((centers >= lo) & (centers <= hi))
    if np.count_nonzero(mask) < 20:
        raise RuntimeError(f"{label}/{config['name']}: too few analysis bins")
    x = np.log(centers[mask])
    raw_y = residuals[mask]
    y = raw_y - np.mean(raw_y)
    free_phase = weighted_linear_sinusoid(x, y, OMEGA)
    locked = fit_phase_locked_waveform(x, y, OMEGA, PHASE_PRED, positive_only=True)
    phase_offset = float(np.angle(np.exp(1j * (free_phase.phase - PHASE_PRED))))
    return {
        "period": label,
        "config": config["name"],
        "background_family": family,
        "degree": degree,
        "spline_smoothing_factor": spline,
        "bins": args.bins,
        "mass_min": args.mass_min,
        "mass_max": args.mass_max,
        "pt_min": args.muon_pt_min,
        "eta_max": args.muon_eta_max,
        "id": "tight" if args.tight_id else "medium",
        "analysis_bins": int(np.count_nonzero(mask)),
        "events_read": int(counters["events_read"]),
        "masses_in_range": int(np.sum(counts)),
        "residual_rms": float(np.sqrt(np.mean(np.square(raw_y)))),
        "residual_max_abs": float(np.max(np.abs(raw_y))),
        "fixed_omega_amplitude": float(free_phase.amplitude),
        "fixed_omega_phase": float(free_phase.phase),
        "phase_offset_from_prediction": phase_offset,
        "fixed_omega_delta_chi2": float(free_phase.delta_chi2),
        "locked_signed_amplitude": float(locked.signed_amplitude),
        "locked_delta_chi2": float(locked.delta_chi2),
        "locked_local_p_one_sided": float(locked.local_p_one_sided_chi2_1dof),
        "positive_sign": bool(locked.signed_amplitude > 0.0),
    }


def main():
    p = argparse.ArgumentParser(description="Prespecified H/G continuum and detector-analysis robustness envelope")
    p.add_argument("--h-manifest", default="data/files_replication.txt")
    p.add_argument("--g-manifest", default="data/files_run2016g.txt")
    p.add_argument("--golden-json", default="data/cert/14221/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt")
    p.add_argument("--max-events", type=int, default=100000)
    p.add_argument("--output-dir", default="results/hg_robustness_envelope")
    p.add_argument("--only", help="Comma-separated configuration names")
    args = p.parse_args()

    configs = suite()
    if args.only:
        wanted = {x.strip() for x in args.only.split(",") if x.strip()}
        configs = [c for c in configs if c["name"] in wanted]
        missing = wanted - {c["name"] for c in configs}
        if missing:
            raise SystemExit(f"Unknown robustness configs: {sorted(missing)}")

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "frozen_suite.json").write_text(json.dumps({"omega": OMEGA, "phase_prediction": PHASE_PRED, "configs": configs}, indent=2) + "\n")

    rows = []
    for config in configs:
        for label, manifest in (("Run2016H_file2", args.h_manifest), ("Run2016G_file1", args.g_manifest)):
            print(f"Running {label}: {config['name']}")
            rows.append(evaluate(label, manifest, args.golden_json, config, args.max_events))

    csv_path = outdir / "robustness_envelope.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    phases = np.asarray([abs(r["phase_offset_from_prediction"]) for r in rows], dtype=float)
    amps = np.asarray([r["locked_signed_amplitude"] for r in rows], dtype=float)
    deltas = np.asarray([r["locked_delta_chi2"] for r in rows], dtype=float)
    summary = {
        "omega": OMEGA,
        "phase_prediction": PHASE_PRED,
        "n_rows": len(rows),
        "n_configs": len(configs),
        "all_positive_sign": bool(np.all(amps > 0.0)),
        "min_locked_signed_amplitude": float(np.min(amps)),
        "min_locked_delta_chi2": float(np.min(deltas)),
        "max_abs_phase_offset_rad": float(np.max(phases)),
        "median_abs_phase_offset_rad": float(np.median(phases)),
        "note": "Retrospective robustness envelope; not new prospective significance.",
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2))
    print("Table:", csv_path)


if __name__ == "__main__":
    main()
