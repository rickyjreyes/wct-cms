#!/usr/bin/env python3
"""
cms_wct_cross_validation.py

Blind/cross-validation pipeline for searching a pre-specified oscillatory or
log-periodic residual signature in CMS NanoAOD dimuon data.

What it does
------------
1. Reads one or more CMS NanoAOD ROOT files with uproot.
2. Optionally applies a CMS golden/validated run-luminosity JSON.
3. Selects muons and builds opposite-sign dimuon pairs.
4. Histograms dimuon invariant mass.
5. Fits a smooth background in log(counts) using a Chebyshev polynomial.
6. Computes Pearson-like residuals:
       r = (N - model) / sqrt(model)
7. Searches residuals for:
       r(log m) ~ a cos(omega log(m/m0)) + b sin(omega log(m/m0))
   using an explicit weighted linear regression at every omega.
8. Reports:
   - best-fit omega
   - amplitude and phase
   - local improvement Delta-chi^2
   - global permutation p-value from the maximum scan statistic
   - a pre-specified/frozen omega test if --frozen-omega is supplied
9. Writes CSV/JSON/PNG outputs for reproducibility.

This is deliberately a discovery-resistant workflow: if you already have a
candidate WCT frequency from GWTC/LHC/JUNO/photodiode analyses, pass that value
with --frozen-omega and treat the frozen test as the primary result.

Dependencies
------------
    pip install numpy scipy matplotlib uproot awkward

Example
-------
    python cms_wct_cross_validation.py \
        --input files.txt \
        --output-dir cms_dimuon_wct \
        --mass-min 2.0 \
        --mass-max 120.0 \
        --bins 350 \
        --log-bins \
        --muon-pt-min 4.0 \
        --muon-eta-max 2.4 \
        --tight-id \
        --fit-degree 7 \
        --omega-min 0.5 \
        --omega-max 80 \
        --omega-steps 3000 \
        --frozen-omega 12.345 \
        --permutations 2000 \
        --seed 20260827

Input may be:
  * one ROOT path/URL:
        --input file.root
  * several ROOT paths/URLs:
        --input a.root b.root root://...
  * a text manifest:
        --input files.txt
    with one path/URL per line.

For remote CMS files, URLs such as root://eospublic.cern.ch//eos/opendata/...
may be used if your uproot/XRootD environment supports them.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Optional

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import uproot
from numpy.polynomial import Chebyshev
from scipy.stats import chi2


# ----------------------------
# Configuration/result classes
# ----------------------------

@dataclass
class ScanResult:
    omega: float
    amplitude: float
    phase: float
    delta_chi2: float
    local_p_chi2_2dof: float
    c0: float
    cos_coeff: float
    sin_coeff: float


@dataclass
class AnalysisSummary:
    events_read: int
    events_after_json: int
    selected_muons: int
    opposite_sign_pairs: int
    masses_in_range: int
    best_scan: dict
    frozen_scan: Optional[dict]
    global_permutation_p: Optional[float]
    frozen_permutation_p: Optional[float]
    mass_min: float
    mass_max: float
    bins: int
    log_bins: bool
    fit_degree: int
    omega_min: float
    omega_max: float
    omega_steps: int
    permutations: int
    seed: int


# ----------------------------
# CLI
# ----------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="CMS NanoAOD dimuon blind residual-signature cross-validation"
    )

    p.add_argument(
        "--input", nargs="+", required=True,
        help="ROOT file(s), XRootD URL(s), or a single text manifest."
    )
    p.add_argument(
        "--tree", default="Events",
        help="NanoAOD tree name (default: Events)."
    )
    p.add_argument(
        "--golden-json", default=None,
        help="Optional CMS validated run/luminosity JSON."
    )
    p.add_argument(
        "--output-dir", default="cms_wct_results",
        help="Output directory."
    )

    # Event/object selection
    p.add_argument("--muon-pt-min", type=float, default=4.0)
    p.add_argument("--muon-eta-max", type=float, default=2.4)
    p.add_argument(
        "--tight-id", action="store_true",
        help="Require Muon_tightId if branch exists."
    )
    p.add_argument(
        "--medium-id", action="store_true",
        help="Require Muon_mediumId if branch exists."
    )
    p.add_argument(
        "--max-events", type=int, default=None,
        help="Optional global event cap for quick tests."
    )
    p.add_argument(
        "--step-size", default="100 MB",
        help="uproot iterate step size, e.g. '100 MB'."
    )

    # Histogram
    p.add_argument("--mass-min", type=float, default=2.0)
    p.add_argument("--mass-max", type=float, default=120.0)
    p.add_argument("--bins", type=int, default=350)
    p.add_argument(
        "--log-bins", action="store_true",
        help="Use logarithmically spaced invariant-mass bins."
    )
    p.add_argument(
        "--min-model-count", type=float, default=5.0,
        help="Only use bins whose fitted expected count exceeds this threshold."
    )

    # Smooth background
    p.add_argument(
        "--fit-degree", type=int, default=7,
        help="Chebyshev degree for log(count) smooth background."
    )
    p.add_argument(
        "--fit-iterations", type=int, default=6,
        help="Robust clipping iterations for background fit."
    )
    p.add_argument(
        "--clip-sigma", type=float, default=3.5,
        help="Clip large residual bins while learning smooth background."
    )
    p.add_argument(
        "--mask-window", action="append", default=[],
        help=(
            "Mass window to exclude from background fit as 'low:high'. "
            "May be repeated, e.g. --mask-window 8.5:11.5"
        )
    )

    # Signature scan
    p.add_argument("--omega-min", type=float, default=0.5)
    p.add_argument("--omega-max", type=float, default=80.0)
    p.add_argument("--omega-steps", type=int, default=3000)
    p.add_argument(
        "--frozen-omega", type=float, default=None,
        help="Pre-specified angular frequency in log(m/m0); primary blind test."
    )
    p.add_argument(
        "--m0", type=float, default=1.0,
        help="Reference mass in GeV inside log(m/m0)."
    )

    # Null tests
    p.add_argument(
        "--permutations", type=int, default=1000,
        help="Residual permutations for empirical global/frozen p-values."
    )
    p.add_argument("--seed", type=int, default=20260827)

    return p.parse_args()


# ----------------------------
# Input helpers
# ----------------------------

def expand_inputs(raw_inputs: list[str]) -> list[str]:
    """Expand a one-file .txt/.list manifest, otherwise return inputs verbatim."""
    if len(raw_inputs) == 1:
        p = Path(raw_inputs[0])
        if p.exists() and p.suffix.lower() in {".txt", ".list", ".manifest"}:
            out = []
            for line in p.read_text(encoding="utf-8").splitlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                out.append(line)
            if not out:
                raise ValueError(f"Manifest is empty: {p}")
            return out
    return raw_inputs


def load_golden_json(path: Optional[str]) -> Optional[dict[int, list[tuple[int, int]]]]:
    """
    CMS certification JSON format:
      {"278820": [[1, 50], [53, 80]], ...}
    """
    if path is None:
        return None

    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    parsed: dict[int, list[tuple[int, int]]] = {}
    for run_str, ranges in raw.items():
        parsed[int(run_str)] = [(int(a), int(b)) for a, b in ranges]
    return parsed


def golden_mask(
    runs: np.ndarray,
    lumis: np.ndarray,
    golden: Optional[dict[int, list[tuple[int, int]]]]
) -> np.ndarray:
    if golden is None:
        return np.ones(len(runs), dtype=bool)

    result = np.zeros(len(runs), dtype=bool)
    # Group by run to avoid a Python loop over every event.
    for run in np.unique(runs):
        ranges = golden.get(int(run))
        if not ranges:
            continue
        idx = (runs == run)
        l = lumis[idx]
        ok = np.zeros(len(l), dtype=bool)
        for lo, hi in ranges:
            ok |= (l >= lo) & (l <= hi)
        result[idx] = ok
    return result


def needed_branches(args: argparse.Namespace) -> list[str]:
    branches = [
        "run",
        "luminosityBlock",
        "Muon_pt",
        "Muon_eta",
        "Muon_phi",
        "Muon_mass",
        "Muon_charge",
    ]
    if args.tight_id:
        branches.append("Muon_tightId")
    if args.medium_id:
        branches.append("Muon_mediumId")
    return branches


# ----------------------------
# Physics reconstruction
# ----------------------------

def invariant_mass_from_pairs(pairs: ak.Array) -> ak.Array:
    """
    Compute invariant mass from NanoAOD pt, eta, phi, mass values.

    E  = sqrt(px^2 + py^2 + pz^2 + m^2)
    M2 = (E1+E2)^2 - |p1+p2|^2
    """
    a = pairs["a"]
    b = pairs["b"]

    def components(mu):
        pt = mu["pt"]
        eta = mu["eta"]
        phi = mu["phi"]
        mass = mu["mass"]

        px = pt * np.cos(phi)
        py = pt * np.sin(phi)
        pz = pt * np.sinh(eta)
        p2 = px * px + py * py + pz * pz
        e = np.sqrt(np.maximum(p2 + mass * mass, 0.0))
        return px, py, pz, e

    px1, py1, pz1, e1 = components(a)
    px2, py2, pz2, e2 = components(b)

    e = e1 + e2
    px = px1 + px2
    py = py1 + py2
    pz = pz1 + pz2

    m2 = e * e - px * px - py * py - pz * pz
    return np.sqrt(np.maximum(m2, 0.0))


def extract_dimuon_masses(
    inputs: list[str],
    args: argparse.Namespace,
    golden: Optional[dict[int, list[tuple[int, int]]]],
) -> tuple[np.ndarray, dict[str, int]]:
    """
    Stream NanoAOD. We form all unique opposite-sign selected muon pairs.

    This avoids choosing a resonance-specific "best pair", which is preferable
    for a broad blind residual-spectrum test.
    """
    branches = needed_branches(args)

    counters = {
        "events_read": 0,
        "events_after_json": 0,
        "selected_muons": 0,
        "opposite_sign_pairs": 0,
    }
    mass_chunks: list[np.ndarray] = []

    # uproot accepts ["file.root:Events", ...]
    sources = [f"{src}:{args.tree}" for src in inputs]

    for arrays in uproot.iterate(
        sources,
        expressions=branches,
        step_size=args.step_size,
        library="ak",
        how=dict,
    ):
        n_chunk = len(arrays["run"])

        if args.max_events is not None:
            remaining = args.max_events - counters["events_read"]
            if remaining <= 0:
                break
            if n_chunk > remaining:
                arrays = {k: v[:remaining] for k, v in arrays.items()}
                n_chunk = remaining

        counters["events_read"] += n_chunk

        runs = ak.to_numpy(arrays["run"])
        lumis = ak.to_numpy(arrays["luminosityBlock"])
        gmask_np = golden_mask(runs, lumis, golden)
        counters["events_after_json"] += int(np.count_nonzero(gmask_np))

        gmask = ak.Array(gmask_np)

        pt = arrays["Muon_pt"][gmask]
        eta = arrays["Muon_eta"][gmask]
        phi = arrays["Muon_phi"][gmask]
        mass = arrays["Muon_mass"][gmask]
        charge = arrays["Muon_charge"][gmask]

        keep = (pt >= args.muon_pt_min) & (abs(eta) <= args.muon_eta_max)

        if args.tight_id:
            if "Muon_tightId" not in arrays:
                raise KeyError("Requested --tight-id but Muon_tightId is absent.")
            keep = keep & arrays["Muon_tightId"][gmask]

        if args.medium_id:
            if "Muon_mediumId" not in arrays:
                raise KeyError("Requested --medium-id but Muon_mediumId is absent.")
            keep = keep & arrays["Muon_mediumId"][gmask]

        mu = ak.zip({
            "pt": pt[keep],
            "eta": eta[keep],
            "phi": phi[keep],
            "mass": mass[keep],
            "charge": charge[keep],
        })

        counters["selected_muons"] += int(ak.sum(ak.num(mu, axis=1)))

        pairs = ak.combinations(mu, 2, fields=["a", "b"])
        os_mask = (pairs["a"]["charge"] * pairs["b"]["charge"]) < 0
        os_pairs = pairs[os_mask]

        counters["opposite_sign_pairs"] += int(ak.sum(ak.num(os_pairs, axis=1)))

        m = invariant_mass_from_pairs(os_pairs)
        flat = ak.to_numpy(ak.flatten(m, axis=None))
        flat = flat[np.isfinite(flat)]

        if flat.size:
            mass_chunks.append(flat.astype(np.float64, copy=False))

        print(
            f"\rEvents: {counters['events_read']:,} | "
            f"OS pairs: {counters['opposite_sign_pairs']:,}",
            end="",
            flush=True,
        )

    print()

    if not mass_chunks:
        raise RuntimeError("No opposite-sign dimuon masses survived selection.")

    return np.concatenate(mass_chunks), counters


# ----------------------------
# Histogram/background model
# ----------------------------

def parse_mask_windows(raw: Iterable[str]) -> list[tuple[float, float]]:
    out = []
    for item in raw:
        try:
            lo_s, hi_s = item.split(":", 1)
            lo, hi = float(lo_s), float(hi_s)
        except Exception as exc:
            raise ValueError(
                f"Invalid --mask-window '{item}'. Expected low:high."
            ) from exc
        if not (lo < hi):
            raise ValueError(f"Invalid mask window {item}: low must be < high.")
        out.append((lo, hi))
    return out


def make_histogram(
    masses: np.ndarray, args: argparse.Namespace
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    in_range = masses[
        (masses >= args.mass_min) &
        (masses <= args.mass_max)
    ]

    if args.log_bins:
        if args.mass_min <= 0:
            raise ValueError("--mass-min must be > 0 with --log-bins.")
        edges = np.geomspace(args.mass_min, args.mass_max, args.bins + 1)
        centers = np.sqrt(edges[:-1] * edges[1:])
    else:
        edges = np.linspace(args.mass_min, args.mass_max, args.bins + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])

    counts, _ = np.histogram(in_range, bins=edges)
    return counts.astype(float), edges, centers


def robust_smooth_background(
    centers: np.ndarray,
    counts: np.ndarray,
    degree: int,
    iterations: int,
    clip_sigma: float,
    excluded_windows: list[tuple[float, float]],
) -> tuple[np.ndarray, np.ndarray]:
    """
    Fit log(count + 0.5) against log(m) with a Chebyshev polynomial.

    Iteratively exclude bins with large Pearson residuals so narrow physical
    resonances do not force the smooth baseline to follow them.

    Returns:
        expected counts
        final mask of bins used to learn background
    """
    if np.any(centers <= 0):
        raise ValueError("Mass centers must be positive for log-space fit.")

    x = np.log(centers)
    y = np.log(counts + 0.5)

    fit_mask = counts > 0
    for lo, hi in excluded_windows:
        fit_mask &= ~((centers >= lo) & (centers <= hi))

    min_points = degree + 3
    if np.count_nonzero(fit_mask) < min_points:
        raise RuntimeError("Too few nonzero/unmasked bins for background fit.")

    model = np.maximum(counts, 0.5)

    for _ in range(iterations):
        if np.count_nonzero(fit_mask) < min_points:
            raise RuntimeError("Robust fit clipped too many bins.")

        # Weight log(count) approximately by sqrt(N).
        w = np.sqrt(np.maximum(counts[fit_mask], 1.0))
        poly = Chebyshev.fit(
            x[fit_mask],
            y[fit_mask],
            deg=degree,
            w=w,
            domain=[x.min(), x.max()],
        )
        log_model = poly(x)
        # Protect against numerical overflow.
        log_model = np.clip(log_model, -30, 50)
        model = np.exp(log_model)

        resid = (counts - model) / np.sqrt(np.maximum(model, 1.0))
        new_mask = fit_mask & (np.abs(resid) <= clip_sigma)

        # Do not allow clipped bins to re-enter in later iterations.
        if np.array_equal(new_mask, fit_mask):
            break
        fit_mask = new_mask

    # Final refit after clipping.
    w = np.sqrt(np.maximum(counts[fit_mask], 1.0))
    poly = Chebyshev.fit(
        x[fit_mask],
        y[fit_mask],
        deg=degree,
        w=w,
        domain=[x.min(), x.max()],
    )
    model = np.exp(np.clip(poly(x), -30, 50))

    return model, fit_mask


# ----------------------------
# Signature statistic
# ----------------------------

def weighted_linear_sinusoid(
    x: np.ndarray,
    y: np.ndarray,
    omega: float,
    weights: Optional[np.ndarray] = None,
) -> ScanResult:
    """
    Compare:
        H0: y = c
        H1: y = c + a cos(omega*x) + b sin(omega*x)

    Delta chi^2 = chi2(H0) - chi2(H1)

    With fixed omega, H1 adds 2 linear parameters. The chi-square p-value is a
    useful asymptotic local reference; permutation p-values are also computed.
    """
    if weights is None:
        weights = np.ones_like(y)

    sw = np.sqrt(weights)

    X0 = np.ones((len(x), 1))
    X1 = np.column_stack([
        np.ones_like(x),
        np.cos(omega * x),
        np.sin(omega * x),
    ])

    X0w = X0 * sw[:, None]
    X1w = X1 * sw[:, None]
    yw = y * sw

    beta0, *_ = np.linalg.lstsq(X0w, yw, rcond=None)
    beta1, *_ = np.linalg.lstsq(X1w, yw, rcond=None)

    r0 = yw - X0w @ beta0
    r1 = yw - X1w @ beta1

    chi0 = float(r0 @ r0)
    chi1 = float(r1 @ r1)
    delta = max(0.0, chi0 - chi1)

    c0, a, b = map(float, beta1)
    amp = math.hypot(a, b)

    # a cos(theta) + b sin(theta) = A cos(theta - phi)
    phase = math.atan2(b, a)

    return ScanResult(
        omega=float(omega),
        amplitude=amp,
        phase=phase,
        delta_chi2=delta,
        local_p_chi2_2dof=float(chi2.sf(delta, df=2)),
        c0=c0,
        cos_coeff=a,
        sin_coeff=b,
    )


def scan_omegas(
    x: np.ndarray,
    residuals: np.ndarray,
    omegas: np.ndarray,
) -> tuple[np.ndarray, ScanResult]:
    scores = np.empty_like(omegas, dtype=float)
    best: Optional[ScanResult] = None

    for i, omega in enumerate(omegas):
        result = weighted_linear_sinusoid(x, residuals, float(omega))
        scores[i] = result.delta_chi2
        if best is None or result.delta_chi2 > best.delta_chi2:
            best = result

    assert best is not None
    return scores, best


def permutation_null(
    x: np.ndarray,
    residuals: np.ndarray,
    omegas: np.ndarray,
    observed_best_score: float,
    frozen_omega: Optional[float],
    observed_frozen_score: Optional[float],
    n_perm: int,
    seed: int,
) -> tuple[Optional[float], Optional[float], np.ndarray, Optional[np.ndarray]]:
    """
    Shuffle residuals among mass bins.

    Global p:
      fraction of permutations whose *maximum scan Delta-chi2* is >= observed.

    Frozen p:
      fraction whose Delta-chi2 at the single predeclared omega is >= observed.

    Add-one correction avoids reporting exactly zero:
      p = (k + 1) / (N + 1)
    """
    if n_perm <= 0:
        return None, None, np.array([]), None

    rng = np.random.default_rng(seed)
    max_scores = np.empty(n_perm, dtype=float)
    frozen_scores = (
        np.empty(n_perm, dtype=float) if frozen_omega is not None else None
    )

    for i in range(n_perm):
        yp = rng.permutation(residuals)

        perm_scores, _ = scan_omegas(x, yp, omegas)
        max_scores[i] = float(np.max(perm_scores))

        if frozen_omega is not None and frozen_scores is not None:
            frozen_scores[i] = weighted_linear_sinusoid(
                x, yp, frozen_omega
            ).delta_chi2

        if (i + 1) % max(1, n_perm // 20) == 0 or i == n_perm - 1:
            print(
                f"\rNull permutations: {i + 1:,}/{n_perm:,}",
                end="",
                flush=True,
            )
    print()

    k_global = int(np.count_nonzero(max_scores >= observed_best_score))
    p_global = (k_global + 1.0) / (n_perm + 1.0)

    p_frozen = None
    if (
        frozen_omega is not None
        and frozen_scores is not None
        and observed_frozen_score is not None
    ):
        k_frozen = int(np.count_nonzero(
            frozen_scores >= observed_frozen_score
        ))
        p_frozen = (k_frozen + 1.0) / (n_perm + 1.0)

    return p_global, p_frozen, max_scores, frozen_scores


# ----------------------------
# Output
# ----------------------------

def save_spectrum_csv(
    outdir: Path,
    edges: np.ndarray,
    centers: np.ndarray,
    counts: np.ndarray,
    model: np.ndarray,
    residuals: np.ndarray,
    analysis_mask: np.ndarray,
    background_fit_mask: np.ndarray,
) -> None:
    arr = np.column_stack([
        edges[:-1],
        edges[1:],
        centers,
        counts,
        model,
        residuals,
        analysis_mask.astype(int),
        background_fit_mask.astype(int),
    ])
    header = (
        "mass_low_GeV,mass_high_GeV,mass_center_GeV,"
        "count,background_model,pearson_residual,"
        "analysis_mask,background_fit_mask"
    )
    np.savetxt(
        outdir / "cms_dimuon_spectrum.csv",
        arr,
        delimiter=",",
        header=header,
        comments="",
    )


def save_scan_csv(
    outdir: Path,
    omegas: np.ndarray,
    scores: np.ndarray,
) -> None:
    np.savetxt(
        outdir / "omega_scan.csv",
        np.column_stack([omegas, scores]),
        delimiter=",",
        header="omega,delta_chi2",
        comments="",
    )


def make_plots(
    outdir: Path,
    centers: np.ndarray,
    counts: np.ndarray,
    model: np.ndarray,
    residuals: np.ndarray,
    analysis_mask: np.ndarray,
    omegas: np.ndarray,
    scores: np.ndarray,
    best: ScanResult,
    frozen: Optional[ScanResult],
    x: np.ndarray,
    y: np.ndarray,
    null_max: np.ndarray,
    args: argparse.Namespace,
) -> None:
    # Spectrum
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.step(centers, counts, where="mid", label="CMS dimuon data")
    ax.plot(centers, model, label="Smooth background")
    ax.set_xlabel("Dimuon invariant mass [GeV]")
    ax.set_ylabel("Counts / bin")
    ax.set_yscale("log")
    if args.log_bins:
        ax.set_xscale("log")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "spectrum_background.png", dpi=180)
    plt.close(fig)

    # Residuals
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.axhline(0.0, linewidth=1)
    ax.plot(
        centers[analysis_mask],
        residuals[analysis_mask],
        marker=".",
        linestyle="none",
    )
    ax.set_xlabel("Dimuon invariant mass [GeV]")
    ax.set_ylabel("(data - model) / sqrt(model)")
    if args.log_bins:
        ax.set_xscale("log")
    fig.tight_layout()
    fig.savefig(outdir / "residuals.png", dpi=180)
    plt.close(fig)

    # Frequency scan
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(omegas, scores)
    ax.axvline(best.omega, linestyle="--", label=f"best omega={best.omega:.5g}")
    if frozen is not None:
        ax.axvline(
            frozen.omega,
            linestyle=":",
            label=f"frozen omega={frozen.omega:.5g}",
        )
    ax.set_xlabel("omega in log(m/m0)")
    ax.set_ylabel("Delta chi-square")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "omega_scan.png", dpi=180)
    plt.close(fig)

    # Best oscillatory fit in log mass
    fit_best = (
        best.c0
        + best.cos_coeff * np.cos(best.omega * x)
        + best.sin_coeff * np.sin(best.omega * x)
    )
    order = np.argsort(x)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(x, y, marker=".", linestyle="none", label="Residuals")
    ax.plot(x[order], fit_best[order], label="Best scan sinusoid")
    ax.set_xlabel("log(m/m0)")
    ax.set_ylabel("Residual")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "best_log_periodic_fit.png", dpi=180)
    plt.close(fig)

    # Null max distribution
    if null_max.size:
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.hist(null_max, bins=40)
        ax.axvline(
            best.delta_chi2,
            linestyle="--",
            label="Observed max",
        )
        ax.set_xlabel("Permutation maximum Delta chi-square")
        ax.set_ylabel("Permutations")
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / "global_null_distribution.png", dpi=180)
        plt.close(fig)


# ----------------------------
# Main analysis
# ----------------------------

def main() -> int:
    args = parse_args()

    if args.mass_min <= 0:
        raise ValueError("mass-min must be positive.")
    if args.mass_max <= args.mass_min:
        raise ValueError("mass-max must exceed mass-min.")
    if args.bins < 20:
        raise ValueError("Use at least 20 bins.")
    if args.fit_degree < 1:
        raise ValueError("fit-degree must be >= 1.")
    if args.omega_steps < 2:
        raise ValueError("omega-steps must be >= 2.")
    if args.omega_max <= args.omega_min:
        raise ValueError("omega-max must exceed omega-min.")
    if args.m0 <= 0:
        raise ValueError("m0 must be positive.")
    if args.tight_id and args.medium_id:
        raise ValueError("Choose at most one of --tight-id and --medium-id.")

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    inputs = expand_inputs(args.input)
    golden = load_golden_json(args.golden_json)
    excluded_windows = parse_mask_windows(args.mask_window)

    # Freeze exact invocation for reproducibility.
    (outdir / "command.txt").write_text(
        " ".join(sys.argv) + "\n",
        encoding="utf-8",
    )
    (outdir / "inputs.txt").write_text(
        "\n".join(inputs) + "\n",
        encoding="utf-8",
    )
    (outdir / "config.json").write_text(
        json.dumps(vars(args), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(f"Reading {len(inputs)} NanoAOD input file(s)...")
    masses, counters = extract_dimuon_masses(inputs, args, golden)

    counts, edges, centers = make_histogram(masses, args)
    masses_in_range = int(np.sum(counts))

    print(f"Masses in requested range: {masses_in_range:,}")

    model, bg_fit_mask = robust_smooth_background(
        centers=centers,
        counts=counts,
        degree=args.fit_degree,
        iterations=args.fit_iterations,
        clip_sigma=args.clip_sigma,
        excluded_windows=excluded_windows,
    )

    residuals = (counts - model) / np.sqrt(np.maximum(model, 1.0))

    analysis_mask = (
        np.isfinite(residuals)
        & np.isfinite(model)
        & (model >= args.min_model_count)
        & (centers > 0)
    )

    if np.count_nonzero(analysis_mask) < 20:
        raise RuntimeError(
            "Too few analysis bins after expected-count threshold."
        )

    x = np.log(centers[analysis_mask] / args.m0)
    y = residuals[analysis_mask]

    # Center residuals. The regression also includes an intercept, but explicit
    # centering makes diagnostic output easier to interpret.
    y = y - np.mean(y)

    omegas = np.linspace(
        args.omega_min,
        args.omega_max,
        args.omega_steps,
    )

    print("Scanning log-periodic frequency...")
    scores, best = scan_omegas(x, y, omegas)

    frozen = None
    if args.frozen_omega is not None:
        frozen = weighted_linear_sinusoid(
            x, y, float(args.frozen_omega)
        )

    print(
        f"Best omega = {best.omega:.8g}, "
        f"DeltaChi2 = {best.delta_chi2:.6g}, "
        f"A = {best.amplitude:.6g}"
    )

    if frozen is not None:
        print(
            f"Frozen omega = {frozen.omega:.8g}, "
            f"DeltaChi2 = {frozen.delta_chi2:.6g}, "
            f"A = {frozen.amplitude:.6g}, "
            f"asymptotic local p = {frozen.local_p_chi2_2dof:.6g}"
        )

    p_global, p_frozen, null_max, null_frozen = permutation_null(
        x=x,
        residuals=y,
        omegas=omegas,
        observed_best_score=best.delta_chi2,
        frozen_omega=args.frozen_omega,
        observed_frozen_score=(
            frozen.delta_chi2 if frozen is not None else None
        ),
        n_perm=args.permutations,
        seed=args.seed,
    )

    if p_global is not None:
        print(f"Global scan permutation p = {p_global:.8g}")
    if p_frozen is not None:
        print(f"Frozen-frequency permutation p = {p_frozen:.8g}")

    save_spectrum_csv(
        outdir,
        edges,
        centers,
        counts,
        model,
        residuals,
        analysis_mask,
        bg_fit_mask,
    )
    save_scan_csv(outdir, omegas, scores)

    if null_max.size:
        np.savetxt(
            outdir / "permutation_global_max_scores.csv",
            null_max,
            delimiter=",",
            header="max_delta_chi2",
            comments="",
        )
    if null_frozen is not None:
        np.savetxt(
            outdir / "permutation_frozen_scores.csv",
            null_frozen,
            delimiter=",",
            header="frozen_delta_chi2",
            comments="",
        )

    summary = AnalysisSummary(
        events_read=counters["events_read"],
        events_after_json=counters["events_after_json"],
        selected_muons=counters["selected_muons"],
        opposite_sign_pairs=counters["opposite_sign_pairs"],
        masses_in_range=masses_in_range,
        best_scan=asdict(best),
        frozen_scan=asdict(frozen) if frozen is not None else None,
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

    (outdir / "summary.json").write_text(
        json.dumps(asdict(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    make_plots(
        outdir=outdir,
        centers=centers,
        counts=counts,
        model=model,
        residuals=residuals,
        analysis_mask=analysis_mask,
        omegas=omegas,
        scores=scores,
        best=best,
        frozen=frozen,
        x=x,
        y=y,
        null_max=null_max,
        args=args,
    )

    print(f"\nDone. Results written to: {outdir.resolve()}")
    print("Primary blind result:")
    if frozen is not None:
        print(
            f"  frozen omega={frozen.omega:.8g}, "
            f"DeltaChi2={frozen.delta_chi2:.6g}, "
            f"permutation p={p_frozen}"
        )
    else:
        print(
            "  No --frozen-omega supplied. The best-frequency scan is "
            "exploratory, not a blind replication."
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
