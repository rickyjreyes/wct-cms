from __future__ import annotations

import argparse

from .analysis import run_analysis


def parser():
    p = argparse.ArgumentParser(description="CMS NanoAOD blind WCT cross-validation")
    p.add_argument("--input", nargs="+", required=True)
    p.add_argument("--tree", default="Events")
    p.add_argument("--golden-json")
    p.add_argument("--output-dir", default="results/cms_dimuon")
    p.add_argument("--muon-pt-min", type=float, default=4.0)
    p.add_argument("--muon-eta-max", type=float, default=2.4)
    p.add_argument("--tight-id", action="store_true")
    p.add_argument("--medium-id", action="store_true")
    p.add_argument("--max-events", type=int)
    p.add_argument("--step-size", default="100 MB")
    p.add_argument("--mass-min", type=float, default=2.0)
    p.add_argument("--mass-max", type=float, default=120.0)
    p.add_argument("--bins", type=int, default=350)
    p.add_argument("--log-bins", action="store_true")
    p.add_argument("--min-model-count", type=float, default=5.0)
    p.add_argument("--fit-degree", type=int, default=7)
    p.add_argument("--fit-iterations", type=int, default=6)
    p.add_argument("--clip-sigma", type=float, default=3.5)
    p.add_argument("--mask-window", action="append", default=[])
    p.add_argument("--omega-min", type=float, default=0.5)
    p.add_argument("--omega-max", type=float, default=80.0)
    p.add_argument("--omega-steps", type=int, default=3000)
    p.add_argument("--frozen-omega", type=float)
    p.add_argument("--m0", type=float, default=1.0)
    p.add_argument("--permutations", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20260827)
    return p


def main():
    args = parser().parse_args()
    if args.tight_id and args.medium_id:
        raise SystemExit("Choose only one of --tight-id or --medium-id")
    summary = run_analysis(args)
    print("Best exploratory omega:", summary.best_scan["omega"])
    print("Cycles across analyzed log-mass span:", summary.best_cycles_across_span)
    print("Residual RMS:", summary.residual_rms)
    print("Maximum |residual|:", summary.residual_max_abs)
    print("Global permutation p:", summary.global_permutation_p)

    if summary.best_scan_at_boundary:
        print(
            "WARNING: best exploratory omega is on a scan boundary; "
            "do not interpret it as a resolved spectral peak."
        )
    if summary.best_cycles_across_span < 1.0:
        print(
            "WARNING: best exploratory omega spans less than one complete cycle; "
            "it is degenerate with broad background structure."
        )
    elif summary.best_cycles_across_span < 2.0:
        print(
            "CAUTION: best exploratory omega spans fewer than two complete cycles; "
            "background degeneracy remains substantial."
        )

    if summary.frozen_scan is None:
        print("No frozen omega supplied: result is exploratory, not blind replication.")
    else:
        print("Frozen omega:", summary.frozen_scan["omega"])
        print("Frozen permutation p:", summary.frozen_permutation_p)


if __name__ == "__main__":
    main()
