#!/usr/bin/env python3
"""Issue #9: read-only comparison of preserved CMS reproduction checkpoints.

This reconstructs the FIRST DIFFERING **ARCHIVED** stage, not necessarily
the first computational operation that diverged. It never edits the frozen
historical outputs or treats numerical agreement as background-robust evidence.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HIST_G2 = "docs/CMS_RUN2016G_FILE2_PHASE_LOCK_RESULT_2026-08-31.json"
AUDIT_G2 = "audit/2026-09-06/baseline/G2_locked/summary.json"
AUDIT_PAIR = "audit/2026-09-06/baseline/pair_observed.json"
HIST_PAIR_SCORE = 108.49776  # Published historical summary, not an invented JSON precision.


def stage_equal(a, b, *, atol=1e-12):
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return math.isclose(a, b, rel_tol=0, abs_tol=atol)
    return a == b


def compare(historic: dict, audit: dict, pair: dict) -> dict:
    stages = [
        ("configuration", ["bins", "log_bins", "mass_min", "mass_max",
                           "max_events", "fit_degree", "fit_iterations",
                           "clip_sigma", "mask_window", "tight_id", "omega",
                           "phase", "permutations", "parametric_bootstrap",
                           "seed"]),
        ("event_counts", ["events_read", "events_after_json", "selected_muons",
                          "opposite_sign_pairs", "masses_in_range", "analysis_bins"]),
        ("residual_summary", ["residual_rms", "residual_max_abs"]),
    ]
    lines = []
    for stage, keys in stages:
        for key in keys:
            left = historic["config"].get(key) if stage == "configuration" else historic.get(key)
            right = audit["config"].get(key) if stage == "configuration" else audit.get(key)
            lines.append(dict(stage=stage, name=key, historical=left, audit=right,
                              equal=stage_equal(left, right),
                              delta=right-left if isinstance(left, (int,float)) and
                                   isinstance(right, (int,float)) and
                                   not isinstance(left,bool) and not isinstance(right,bool) else None))
    for key in ("signed_amplitude", "delta_chi2", "c0"):
        l, r = historic["observed"].get(key), audit["observed"].get(key)
        lines.append(dict(stage="locked_result", name=key, historical=l, audit=r,
                          equal=stage_equal(l,r),
                          delta=r-l if isinstance(l,(float,int)) and isinstance(r,(float,int)) else None))
    new = pair["primary_pair_score_min_locked_delta_chi2"]
    lines.append(dict(stage="pair_pipeline", name="selected_pair_Q",
                      historical=HIST_PAIR_SCORE, audit=new,
                      equal=stage_equal(HIST_PAIR_SCORE,new), delta=new-HIST_PAIR_SCORE))
    first = next((v for v in lines if not v["equal"]), None)
    # Input locations differ by construction; step_size is a real configuration
    # difference, not proof of different selected events. It is kept separate
    # from comparable scientific parameters and must be tested experimentally.
    return {
        "schema": "CMS_NUMERIC_DISCREPANCY_ARCHIVE_COMPARISON_V1",
        "classification": "BOUNDED_OUTPUT_DISCREPANCY_CAUSE_UNRESOLVED",
        "first_differing_archived_stage": first["stage"] if first else None,
        "first_differing_archived_field": first["name"] if first else None,
        "matched_selected_counts": all(row["equal"] for row in lines
            if row["stage"]=="event_counts"),
        "pair_historical_precision": "Approximately 108.49776 from historical summary",
        "pair_absolute_difference": abs(new-HIST_PAIR_SCORE),
        "pair_relative_difference_pct": 100*abs(new-HIST_PAIR_SCORE)/HIST_PAIR_SCORE,
        "observed_pair_winner": pair["selected_background"]["name"],
        "historical_pair_winner": "spline_s2",
        "historical_step_size": historic["config"].get("step_size"),
        "audit_step_size": audit["config"].get("step_size"),
        "step_size_is_a_cause_proven": False,
        "full_original_root_checksums_verified_by_this_script": False,
        "first_computational_divergence_proven": False,
        "physical_significance_inferred": False,
        "comparison": lines,
        "remaining_work": (
            "Re-run checksum-verified original ROOT files with both step sizes "
            "and recorded dependency environments; compare event keys, float64 "
            "mass vectors, histograms/bin assignments, fit masks, continuum "
            "predictions, residual arrays, blocked CV folds, and final scores "
            "at the first differing intermediate byte/stage. Historical input "
            "arrays necessary for that proof are not present in this comparison."
        )
    }


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",type=Path,default=ROOT)
    p.add_argument("--out",type=Path,default=None)
    a=p.parse_args(argv)
    h=json.loads((a.root/HIST_G2).read_text())
    x=json.loads((a.root/AUDIT_G2).read_text())
    pair=json.loads((a.root/AUDIT_PAIR).read_text())
    result=compare(h,x,pair)
    out=a.out or a.root/"results"/"cms_numeric_archive_comparison.json"
    if out.resolve() in {(a.root/path).resolve() for path in (HIST_G2,AUDIT_G2,AUDIT_PAIR)}:
        p.error("Refusing to overwrite archived source evidence")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(f"{result['classification']}: first differing archived field "
          f"{result['first_differing_archived_stage']}/"
          f"{result['first_differing_archived_field']}; report {out}")
    return 0


if __name__=="__main__":
    raise SystemExit(main())
