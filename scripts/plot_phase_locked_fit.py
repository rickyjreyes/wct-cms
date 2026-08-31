#!/usr/bin/env python3
"""Render the publication figure for the prospective CMS G2 holdout.

This script is deliberately a *presentation-only* layer.  It does not rerun
selection, background fitting, frequency selection, phase fitting, or the
phase-locked hypothesis test.  It reads the frozen outputs written by
``scripts/run_phase_locked_period.py`` and renders those saved results.

Required inputs in ``--input-dir``
---------------------------------
``summary.json``
    Frozen omega/phase, fitted signed amplitude, Delta chi-square, local
    analytic p-value, empirical null p-values, residual diagnostics, and the
    run configuration.

``spectrum.csv``
    Saved histogram, background, Pearson residuals, and the exact analysis
    mask used by the prospective test.

Outputs
-------
``phase_locked_fit.png``
    High-resolution raster figure used directly by the LaTeX manuscript.

``phase_locked_fit.pdf``
    Vector version suitable for publication and archival use.

The figure keeps the distinction between the very small *local analytic*
p-value and the finite-resolution permutation/refit-bootstrap null results
visible in the annotation so the plot cannot be read as an empirical 11-sigma
claim.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Plot the saved CMS Run2016G G2 prospective phase-locked result "
            "without rerunning the statistical analysis."
        )
    )
    p.add_argument(
        "--input-dir",
        default="results/run2016g_file2_phase_locked",
        help="Directory containing summary.json and spectrum.csv",
    )
    p.add_argument(
        "--output",
        default=None,
        help=(
            "PNG output path. Default: <input-dir>/phase_locked_fit.png. "
            "A PDF with the same stem is written automatically."
        ),
    )
    p.add_argument(
        "--dpi",
        type=int,
        default=400,
        help="PNG resolution (default: 400)",
    )
    p.add_argument(
        "--title",
        default="CMS Run2016G G2 — Prospective Phase-Locked Holdout",
        help="Figure title",
    )
    return p


def load_outputs(input_dir: Path):
    summary_path = input_dir / "summary.json"
    spectrum_path = input_dir / "spectrum.csv"

    if not summary_path.exists():
        raise FileNotFoundError(f"Missing {summary_path}")
    if not spectrum_path.exists():
        raise FileNotFoundError(f"Missing {spectrum_path}")

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    spectrum = np.genfromtxt(
        spectrum_path,
        delimiter=",",
        names=True,
        dtype=float,
        encoding="utf-8",
    )
    return summary, spectrum


def _format_p(value: float | None) -> str:
    if value is None:
        return "n/a"
    return f"{float(value):.3g}"


def _empirical_label(summary: dict) -> str:
    """Compactly report the actual Monte Carlo resolution in the saved run."""
    n_perm = int(summary.get("permutations", 0) or 0)
    n_boot = int(summary.get("parametric_bootstrap", 0) or 0)
    p_perm = summary.get("permutation_p")
    p_boot = summary.get("parametric_bootstrap_p")

    pieces: list[str] = []
    if n_perm > 0 and p_perm is not None:
        floor = 1.0 / (n_perm + 1.0)
        if np.isclose(float(p_perm), floor, rtol=1e-10, atol=0.0):
            pieces.append(f"perm: 0/{n_perm} (floor 1/{n_perm + 1})")
        else:
            pieces.append(f"perm: p={_format_p(p_perm)}")

    if n_boot > 0 and p_boot is not None:
        floor = 1.0 / (n_boot + 1.0)
        if np.isclose(float(p_boot), floor, rtol=1e-10, atol=0.0):
            pieces.append(f"refit: 0/{n_boot} (floor 1/{n_boot + 1})")
        else:
            pieces.append(f"refit: p={_format_p(p_boot)}")

    return "; ".join(pieces) if pieces else "empirical null: n/a"


def _mass_to_x(mass, m0: float):
    arr = np.asarray(mass, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log(arr / m0)


def _x_to_mass(x, m0: float):
    return m0 * np.exp(np.asarray(x, dtype=float))


def main() -> None:
    args = parser().parse_args()
    input_dir = Path(args.input_dir)
    summary, spectrum = load_outputs(input_dir)

    if summary.get("test") != "fixed_omega_fixed_phase_positive_amplitude":
        raise RuntimeError(
            "summary.json is not a fixed-omega, fixed-phase, positive-amplitude run"
        )

    observed = summary["observed"]
    omega = float(summary["omega"])
    phase = float(summary["phase"])
    amplitude = float(observed["signed_amplitude"])
    c0 = float(observed["c0"])
    delta_chi2 = float(observed["delta_chi2"])
    local_p = float(observed["local_p_one_sided_chi2_1dof"])

    config = summary.get("config", {})
    m0 = float(config.get("m0", 1.0))

    mass = np.asarray(spectrum["mass_center"], dtype=float)
    residual = np.asarray(spectrum["residual"], dtype=float)
    analysis_mask = np.asarray(spectrum["analysis_mask"], dtype=float) > 0.5

    valid = (
        analysis_mask
        & np.isfinite(mass)
        & np.isfinite(residual)
        & (mass > 0.0)
    )
    if np.count_nonzero(valid) < 3:
        raise RuntimeError("Too few valid analysis bins in spectrum.csv")

    # The prospective runner centers the selected Pearson residual field before
    # fitting c + A cos(omega*x - phase).  Reproduce only that saved plotting
    # coordinate here; do not refit any parameter.
    x = np.log(mass[valid] / m0)
    y_raw = residual[valid]
    y = y_raw - np.mean(y_raw)

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    x_dense = np.linspace(float(np.min(x)), float(np.max(x)), 2400)
    locked_curve = c0 + amplitude * np.cos(omega * x_dense - phase)

    # Single-panel publication layout.  A top secondary axis gives the physical
    # dimuon mass while the primary axis preserves the exact tested log-mass
    # coordinate.
    fig, ax = plt.subplots(figsize=(8.0, 5.1), constrained_layout=True)

    ax.scatter(
        x,
        y,
        s=16,
        alpha=0.62,
        linewidths=0,
        rasterized=True,
        label="Centered Pearson residuals",
        zorder=2,
    )
    ax.plot(
        x_dense,
        locked_curve,
        linewidth=2.1,
        label="Exact frozen waveform",
        zorder=4,
    )
    ax.axhline(0.0, linewidth=0.9, alpha=0.75, zorder=1)

    ax.set_xlim(float(np.min(x)), float(np.max(x)))
    ax.set_xlabel(r"$x=\ln(m_{\mu\mu}/1\,\mathrm{GeV})$")
    ax.set_ylabel("Centered Pearson residual")
    ax.set_title(args.title, pad=12)

    # Physical-mass reference axis.  These ticks are deliberately sparse so the
    # log-coordinate remains the visually dominant tested variable.
    mass_ticks = np.asarray([2, 3, 5, 10, 20, 50, 100, 120], dtype=float)
    mass_ticks = mass_ticks[(mass_ticks >= np.min(mass[valid])) & (mass_ticks <= np.max(mass[valid]))]
    secax = ax.secondary_xaxis(
        "top",
        functions=(
            lambda xx: _x_to_mass(xx, m0),
            lambda mm: _mass_to_x(mm, m0),
        ),
    )
    secax.set_xlabel(r"$m_{\mu\mu}$ [GeV]", labelpad=7)
    if mass_ticks.size:
        secax.set_xticks(mass_ticks)
        secax.set_xticklabels([f"{v:g}" for v in mass_ticks])

    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which="both", direction="in", right=True)
    ax.tick_params(which="major", length=5.0)
    ax.tick_params(which="minor", length=2.5)
    secax.tick_params(which="both", direction="in")

    # Light horizontal guide lines aid reading while avoiding a grid that could
    # visually mimic an oscillatory pattern.
    ax.grid(axis="y", which="major", alpha=0.14, linewidth=0.7)

    empirical = _empirical_label(summary)
    annotation = (
        rf"Frozen: $\omega={omega:.10f}$, $\phi={phase:.10f}$ rad, $A>0$" "\n"
        rf"Observed: $A={amplitude:.4f}$, $\Delta\chi^2={delta_chi2:.2f}$" "\n"
        rf"Local analytic: $p={local_p:.2e}$" "\n"
        rf"Empirical null: {empirical}" "\n"
        rf"Residual RMS={float(summary['residual_rms']):.3f}; "
        rf"max $|r|$={float(summary['residual_max_abs']):.3f}"
    )
    ax.text(
        0.015,
        0.975,
        annotation,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8.7,
        linespacing=1.25,
        bbox={"boxstyle": "round,pad=0.4", "alpha": 0.88},
        zorder=6,
    )

    ax.legend(
        loc="lower right",
        frameon=True,
        fontsize=9,
        handlelength=2.6,
    )

    # Small on-figure provenance note: the plot is rendered from frozen saved
    # outputs and is not a new optimization or statistical fit.
    fig.text(
        0.995,
        0.005,
        "Rendered from frozen summary.json + spectrum.csv; no parameter refit",
        ha="right",
        va="bottom",
        fontsize=7.2,
        alpha=0.7,
    )

    if args.output is None:
        png_path = input_dir / "phase_locked_fit.png"
    else:
        png_path = Path(args.output)

    png_path.parent.mkdir(parents=True, exist_ok=True)
    pdf_path = png_path.with_suffix(".pdf")

    metadata = {
        "source_summary": str(input_dir / "summary.json"),
        "source_spectrum": str(input_dir / "spectrum.csv"),
        "omega_fixed": omega,
        "phase_fixed": phase,
        "positive_amplitude_required": True,
        "signed_amplitude": amplitude,
        "delta_chi2": delta_chi2,
        "local_analytic_p": local_p,
        "permutation_p": summary.get("permutation_p"),
        "parametric_bootstrap_p": summary.get("parametric_bootstrap_p"),
        "analysis_bins": int(summary["analysis_bins"]),
        "plot_note": "Presentation-only rendering; no signal parameter was refit.",
    }

    fig.savefig(png_path, dpi=args.dpi, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    metadata_path = png_path.with_name(png_path.stem + "_metadata.json")
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"Wrote: {png_path}")
    print(f"Wrote: {pdf_path}")
    print(f"Wrote: {metadata_path}")
    print(f"Locked omega: {omega}")
    print(f"Locked phase: {phase}")
    print(f"Signed amplitude: {amplitude}")
    print(f"Delta chi-square: {delta_chi2}")
    print(f"Local one-sided p: {local_p}")
    print(f"Empirical null: {empirical}")


if __name__ == "__main__":
    main()
