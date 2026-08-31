from __future__ import annotations

import argparse

from cms_wct.significance import (
    add_one_monte_carlo_p,
    clopper_pearson_upper_bound,
    five_sigma_requirements,
    one_sided_gaussian_p,
)


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Plan or evaluate direct Monte Carlo resolution for a one-sided Gaussian significance target"
    )
    p.add_argument("--z", type=float, default=5.0, help="One-sided Gaussian target Z (default: 5)")
    p.add_argument("--trials", type=int, help="Completed null pseudoexperiments")
    p.add_argument("--exceedances", type=int, default=0, help="Null trials at or above the observed statistic")
    p.add_argument("--confidence", type=float, default=0.95, help="One-sided exact binomial confidence level")
    return p


def main() -> None:
    args = parser().parse_args()
    target_p = one_sided_gaussian_p(args.z)

    if args.z == 5.0:
        req = five_sigma_requirements()
        print("One-sided 5 sigma p threshold:", req["one_sided_p"])
        print("Minimum trials for add-one numerical resolution:", req["minimum_trials_add_one_resolution"])
        print("Zero-exceedance trials for 95% upper bound <= 5 sigma p:", req["minimum_zero_exceedance_trials_95pct"])
        print("Zero-exceedance trials for 99% upper bound <= 5 sigma p:", req["minimum_zero_exceedance_trials_99pct"])
    else:
        print(f"One-sided {args.z:g} sigma p threshold:", target_p)

    if args.trials is None:
        return
    if args.trials <= 0:
        raise SystemExit("--trials must be positive")
    if args.exceedances < 0 or args.exceedances > args.trials:
        raise SystemExit("--exceedances must satisfy 0 <= exceedances <= trials")

    p_add_one = add_one_monte_carlo_p(args.exceedances, args.trials)
    p_upper = clopper_pearson_upper_bound(args.exceedances, args.trials, args.confidence)
    print("Observed null trials:", args.trials)
    print("Observed exceedances:", args.exceedances)
    print("Add-one Monte Carlo p:", p_add_one)
    print(f"Exact one-sided {100.0 * args.confidence:g}% upper bound on null exceedance probability:", p_upper)
    print("Target p:", target_p)
    print("Add-one numerical target reached:", p_add_one <= target_p)
    print("Confidence-bound target established:", p_upper <= target_p)


if __name__ == "__main__":
    main()
