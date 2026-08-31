from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


def as_bool(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def main():
    p = argparse.ArgumentParser(description="Summarize the existing H/G robustness envelope without rerunning ROOT extraction")
    p.add_argument("--input", default="results/hg_robustness_envelope/robustness_envelope.csv")
    p.add_argument("--output", default="results/hg_robustness_envelope/detailed_summary.json")
    args = p.parse_args()

    with Path(args.input).open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit("No robustness rows found")

    for r in rows:
        for key in (
            "locked_signed_amplitude",
            "locked_delta_chi2",
            "phase_offset_from_prediction",
            "residual_rms",
            "residual_max_abs",
            "fixed_omega_delta_chi2",
        ):
            r[key] = float(r[key])
        r["positive_sign"] = as_bool(r["positive_sign"])

    nonpositive = [r for r in rows if not r["positive_sign"]]
    phase_gt_05 = [r for r in rows if abs(r["phase_offset_from_prediction"]) > 0.5]
    phase_gt_10 = [r for r in rows if abs(r["phase_offset_from_prediction"]) > 1.0]
    worst_phase = sorted(rows, key=lambda r: abs(r["phase_offset_from_prediction"]), reverse=True)[:10]
    weakest = sorted(rows, key=lambda r: r["locked_delta_chi2"])[:10]

    by_period = defaultdict(lambda: {"rows": 0, "positive": 0, "nonpositive": 0, "phase_gt_0p5": 0, "phase_gt_1p0": 0})
    for r in rows:
        d = by_period[r["period"]]
        d["rows"] += 1
        d["positive"] += int(r["positive_sign"])
        d["nonpositive"] += int(not r["positive_sign"])
        d["phase_gt_0p5"] += int(abs(r["phase_offset_from_prediction"]) > 0.5)
        d["phase_gt_1p0"] += int(abs(r["phase_offset_from_prediction"]) > 1.0)

    config_periods = defaultdict(list)
    for r in rows:
        config_periods[r["config"]].append(r)
    both_positive = []
    any_nonpositive = []
    for config, cr in sorted(config_periods.items()):
        if all(r["positive_sign"] for r in cr):
            both_positive.append(config)
        else:
            any_nonpositive.append(config)

    def compact(r):
        return {
            "period": r["period"],
            "config": r["config"],
            "background_family": r["background_family"],
            "locked_signed_amplitude": r["locked_signed_amplitude"],
            "locked_delta_chi2": r["locked_delta_chi2"],
            "phase_offset_from_prediction": r["phase_offset_from_prediction"],
            "fixed_omega_delta_chi2": r["fixed_omega_delta_chi2"],
            "residual_rms": r["residual_rms"],
            "residual_max_abs": r["residual_max_abs"],
        }

    summary = {
        "n_rows": len(rows),
        "n_configs": len(config_periods),
        "n_nonpositive_rows": len(nonpositive),
        "nonpositive_rows": [compact(r) for r in nonpositive],
        "n_phase_offset_gt_0p5_rad": len(phase_gt_05),
        "n_phase_offset_gt_1p0_rad": len(phase_gt_10),
        "worst_phase_rows": [compact(r) for r in worst_phase],
        "weakest_locked_rows": [compact(r) for r in weakest],
        "period_summary": dict(by_period),
        "n_configs_positive_in_both_periods": len(both_positive),
        "configs_positive_in_both_periods": both_positive,
        "configs_with_any_nonpositive_period": any_nonpositive,
        "background_family_counts": dict(Counter(r["background_family"] for r in rows)),
        "interpretation": "Retrospective robustness diagnostic only. Nonpositive rows and large phase drifts are failures of uniform waveform robustness and must remain visible in reporting.",
    }

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("Rows:", summary["n_rows"])
    print("Configs:", summary["n_configs"])
    print("Nonpositive rows:", summary["n_nonpositive_rows"])
    for r in summary["nonpositive_rows"]:
        print("  FAIL SIGN:", r["period"], r["config"], "A=", r["locked_signed_amplitude"], "phase_offset=", r["phase_offset_from_prediction"])
    print("|phase offset| > 0.5 rad:", summary["n_phase_offset_gt_0p5_rad"])
    print("|phase offset| > 1.0 rad:", summary["n_phase_offset_gt_1p0_rad"])
    print("Configs positive in both periods:", summary["n_configs_positive_in_both_periods"], "/", summary["n_configs"])
    print("Detailed summary:", out)


if __name__ == "__main__":
    main()
