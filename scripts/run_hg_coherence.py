from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from cms_wct.background import make_histogram
from cms_wct.background_families import fit_background_family
from cms_wct.cmsio import expand_inputs, extract_dimuon_masses, load_golden_json
from cms_wct.locked import common_waveform_coherence, paired_permutation_coherence_null

OMEGA = 7.025825825825827
MASKS = [(2.9, 3.3), (3.55, 3.85), (8.5, 11.5), (80.0, 100.0)]


def make_args(manifest, max_events):
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


def residual_vector(manifest, golden_path, max_events):
    args = make_args(manifest, max_events)
    golden = load_golden_json(golden_path)
    inputs = expand_inputs(args.input)
    masses, counters = extract_dimuon_masses(inputs, args, golden)
    counts, _, centers = make_histogram(masses, 2.0, 120.0, 350, True)
    model, _ = fit_background_family(
        centers,
        counts,
        family="chebyshev",
        degree=7,
        iterations=6,
        clip_sigma=3.5,
        excluded_windows=MASKS,
    )
    residuals = (counts - model) / np.sqrt(np.maximum(model, 1.0))
    mask = np.isfinite(residuals) & (model >= 5.0) & (centers > 0)
    for lo, hi in MASKS:
        mask &= ~((centers >= lo) & (centers <= hi))
    x = np.log(centers[mask])
    raw_y = residuals[mask]
    y = raw_y - np.mean(raw_y)
    diagnostics = {
        "events_read": int(counters["events_read"]),
        "masses_in_range": int(np.sum(counts)),
        "analysis_bins": int(np.count_nonzero(mask)),
        "residual_rms": float(np.sqrt(np.mean(np.square(raw_y)))),
        "residual_max_abs": float(np.max(np.abs(raw_y))),
    }
    return x, y, diagnostics


def main():
    p = argparse.ArgumentParser(description="Retrospective Run2016H/G common-waveform coherence calibration")
    p.add_argument("--h-manifest", default="data/files_replication.txt")
    p.add_argument("--g-manifest", default="data/files_run2016g.txt")
    p.add_argument("--golden-json", default="data/cert/14221/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt")
    p.add_argument("--max-events", type=int, default=100000)
    p.add_argument("--permutations", type=int, default=10000)
    p.add_argument("--seed", type=int, default=20260901)
    p.add_argument("--output-dir", default="results/hg_coherence")
    args = p.parse_args()

    xh, yh, dh = residual_vector(args.h_manifest, args.golden_json, args.max_events)
    xg, yg, dg = residual_vector(args.g_manifest, args.golden_json, args.max_events)

    observed = common_waveform_coherence(xh, yh, xg, yg, OMEGA)
    p_empirical, null_scores = paired_permutation_coherence_null(
        xh,
        yh,
        xg,
        yg,
        OMEGA,
        observed["coherence_score"],
        args.permutations,
        args.seed,
    )

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    summary = {
        "classification": "retrospective_common_waveform_diagnostic",
        "omega": OMEGA,
        "statistic_definition": "T_coh = common_score - heterogeneity = 2*common_score - separate_score",
        "observed": observed,
        "empirical_p": p_empirical,
        "permutations": args.permutations,
        "seed": args.seed,
        "run2016h": dh,
        "run2016g": dg,
        "warning": "H and G were already observed before this coherence statistic was introduced; use as calibrated retrospective evidence, not prospective replication.",
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    np.savetxt(outdir / "coherence_null.csv", null_scores, delimiter=",", header="coherence_score", comments="")

    print("Common score:", observed["common_score"])
    print("Heterogeneity:", observed["heterogeneity"])
    print("Coherence score:", observed["coherence_score"])
    print("Phase difference (rad):", observed["phase_difference"])
    print("Empirical paired-null p:", p_empirical)
    print("Summary:", outdir / "summary.json")


if __name__ == "__main__":
    main()
