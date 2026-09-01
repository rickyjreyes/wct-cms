from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from cms_wct.background import make_histogram
from cms_wct.background_cv import frozen_background_candidates
from cms_wct.background_kill import (
    PairKillConfig,
    deterministic_absorption_matrix,
    injection_recovery_pair,
    spurious_signal_null_pair,
)
from cms_wct.cmsio import expand_inputs, extract_dimuon_masses, load_golden_json

OMEGA = 7.025825825825827
PHASE = -0.2313916852932179
MASKS = ((2.9, 3.3), (3.55, 3.85), (8.5, 11.5), (80.0, 100.0))


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
    print(f"Extracting {label}")
    masses, counters = extract_dimuon_masses(inputs, args, golden)
    counts, edges, centers = make_histogram(masses, 2.0, 120.0, 350, True)
    return {
        "label": label,
        "counts": np.asarray(counts, dtype=float),
        "edges": np.asarray(edges, dtype=float),
        "centers": np.asarray(centers, dtype=float),
        "counters": counters,
    }


def write_rows(path: Path, rows: list[dict]):
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def compact_without_trials(record: dict) -> dict:
    return {key: value for key, value in record.items() if key != "trials"}


def main():
    p = argparse.ArgumentParser(
        description=(
            "CMS Background Kill Test: WCT-blind predictive background selection, "
            "end-to-end spurious-signal null, and signal-injection/absorption calibration"
        )
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
    p.add_argument("--null-trials", type=int, default=200)
    p.add_argument("--injection-trials", type=int, default=100)
    p.add_argument(
        "--injection-amplitudes",
        type=float,
        nargs="+",
        default=[0.25, 0.5, 0.75, 1.0],
    )
    p.add_argument("--seed", type=int, default=20260901)
    p.add_argument(
        "--no-reselect",
        action="store_true",
        help=(
            "Do not rerun WCT-blind background CV inside pseudoexperiments. "
            "This is faster but is not the preferred end-to-end kill test."
        ),
    )
    p.add_argument("--output-dir", default="results/cms_background_kill")
    args = p.parse_args()

    if args.null_trials < 1:
        p.error("--null-trials must be positive")
    if args.injection_trials < 1:
        p.error("--injection-trials must be positive")
    if any(a <= 0.0 for a in args.injection_amplitudes):
        p.error("--injection-amplitudes must all be positive")

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    candidates = frozen_background_candidates()
    config = PairKillConfig(
        omega=OMEGA,
        phase=PHASE,
        excluded_windows=MASKS,
        n_folds=args.folds,
        block_size=args.block_size,
        iterations=6,
        clip_sigma=3.5,
        min_model_count=5.0,
    )
    reselect = not args.no_reselect

    freeze = {
        "classification": "retrospective_background_kill_test",
        "purpose": (
            "Test whether signal-independent flexible-background selection and detrending "
            "can manufacture or absorb the frozen CMS waveform."
        ),
        "frequency_or_phase_used_for_background_selection": False,
        "config": config.to_dict(),
        "candidate_order": [candidate.to_dict() for candidate in candidates],
        "h_manifest": args.h_manifest,
        "g_manifest": args.g_manifest,
        "golden_json": args.golden_json,
        "max_events": args.max_events,
        "null_trials": args.null_trials,
        "injection_trials_per_amplitude": args.injection_trials,
        "injection_amplitudes_residual_units": list(map(float, args.injection_amplitudes)),
        "reselect_background_each_pseudoexperiment": reselect,
        "seed_spurious_null": args.seed,
        "seed_injection": args.seed + 1,
        "primary_statistic": (
            "minimum of the two positive-only phase-locked Delta chi-square values; "
            "this requires both periods to contribute"
        ),
        "claim_lock": (
            "This is retrospective because H/G have already been inspected. "
            "Its purpose is model identifiability, not prospective replication."
        ),
    }
    (outdir / "selection_freeze.json").write_text(
        json.dumps(freeze, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    golden = load_golden_json(args.golden_json)
    h = load_histogram("Run2016H_file2", args.h_manifest, golden, args.max_events)
    g = load_histogram("Run2016G_file1", args.g_manifest, golden, args.max_events)

    print("Running end-to-end spurious-signal null")
    spurious = spurious_signal_null_pair(
        h["centers"],
        h["counts"],
        g["centers"],
        g["counts"],
        candidates,
        config,
        n_trials=args.null_trials,
        seed=args.seed,
        reselect_each_trial=reselect,
    )
    null_threshold = float(spurious["primary_score_quantiles"]["q95"])

    print("Running deterministic signal-absorption matrix")
    absorption = deterministic_absorption_matrix(
        h["centers"],
        h["counts"],
        g["centers"],
        g["counts"],
        candidates,
        config,
        injection_amplitudes=args.injection_amplitudes,
    )

    print("Running end-to-end signal injection/recovery")
    injection = injection_recovery_pair(
        h["centers"],
        h["counts"],
        g["centers"],
        g["counts"],
        candidates,
        config,
        injection_amplitudes=args.injection_amplitudes,
        n_trials_per_amplitude=args.injection_trials,
        seed=args.seed + 1,
        null_primary_threshold=null_threshold,
        reselect_each_trial=reselect,
    )

    observed = spurious["observed"]
    cv_rows = observed["cv_scores"]
    write_rows(outdir / "cv_scores.csv", cv_rows)
    write_rows(outdir / "spurious_null_trials.csv", spurious["trials"])
    write_rows(outdir / "absorption_matrix.csv", absorption["rows"])
    write_rows(outdir / "injection_trials.csv", injection["trials"])

    (outdir / "spurious_null_summary.json").write_text(
        json.dumps(compact_without_trials(spurious), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (outdir / "injection_summary.json").write_text(
        json.dumps(compact_without_trials(injection), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    worst_absorption = {}
    for amplitude in args.injection_amplitudes:
        rows = [
            row for row in absorption["rows"]
            if np.isclose(row["injected_amplitude"], float(amplitude))
        ]
        worst = min(rows, key=lambda row: row["min_retention"])
        worst_absorption[str(float(amplitude))] = {
            "candidate": worst["analysis_candidate"],
            "min_retention": float(worst["min_retention"]),
        }

    summary = {
        "classification": "retrospective_cms_background_kill_test",
        "config": config.to_dict(),
        "selected_background": observed["selected_background"],
        "winner_ties_within_tolerance": observed["winner_ties_within_tolerance"],
        "observed_primary_pair_score": observed[
            "primary_pair_score_min_locked_delta_chi2"
        ],
        "observed_secondary_pair_score": observed[
            "secondary_pair_score_sum_locked_delta_chi2"
        ],
        "observed_positive_sign_in_both": observed["positive_sign_in_both"],
        "observed_spectrum_a": observed["spectrum_a"],
        "observed_spectrum_b": observed["spectrum_b"],
        "spurious_signal_null": {
            "n_trials": spurious["n_trials"],
            "primary_exceedances": spurious["primary_exceedances"],
            "primary_add_one_p": spurious["primary_add_one_p"],
            "secondary_exceedances": spurious["secondary_exceedances"],
            "secondary_add_one_p": spurious["secondary_add_one_p"],
            "primary_score_quantiles": spurious["primary_score_quantiles"],
            "selected_background_counts": spurious["selected_background_counts"],
        },
        "signal_injection_recovery": injection["summaries"],
        "worst_deterministic_absorption_by_injected_amplitude": worst_absorption,
        "interpretation_rules": [
            (
                "If smooth pseudo-spectra frequently produce a primary score at least "
                "as large as observed, the detrending/background-selection pipeline can "
                "manufacture the frozen waveform and the baseline evidence is weakened."
            ),
            (
                "If predictively competitive backgrounds retain only a small fraction of "
                "injected amplitude, the pipeline lacks power to distinguish a real waveform "
                "from continuum flexibility at that amplitude."
            ),
            (
                "A low spurious-null p-value is not sufficient by itself: acceptable "
                "injection recovery is required so a flexible continuum cannot simply erase "
                "the alternative."
            ),
            (
                "This H/G analysis is retrospective and cannot become a new prospective "
                "replication regardless of the result."
            ),
        ],
    }
    (outdir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    print("Selected background:", observed["selected_background"]["name"])
    print(
        "Observed primary min locked DeltaChi2:",
        observed["primary_pair_score_min_locked_delta_chi2"],
    )
    print(
        "Spurious-null exceedances:",
        spurious["primary_exceedances"],
        "/",
        spurious["n_trials"],
        "add-one p =",
        spurious["primary_add_one_p"],
    )
    print("Null 95% primary-score threshold:", null_threshold)
    for row in injection["summaries"]:
        print(
            "Injection",
            row["injected_amplitude"],
            "median min retention =",
            row["median_min_retention"],
            "power =",
            row.get("power_above_null_threshold"),
        )
    print("Summary:", outdir / "summary.json")


if __name__ == "__main__":
    main()
