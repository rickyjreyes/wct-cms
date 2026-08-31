#!/usr/bin/env python3
"""Plot the prospective CMS G2 phase-locked holdout result.

Reads the output files produced by scripts/run_phase_locked_period.py:
  - summary.json
  - spectrum.csv

Writes:
  - phase_locked_fit.png
  - phase_locked_fit.pdf

The plotted residuals use exactly the bins included in the frozen primary test.
The residual field is centered in the same way as the phase-locked fit before
the frozen waveform is overlaid.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Plot frozen-frequency, frozen-phase CMS residual fit."
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
            "A PDF with the same stem is also written."
        ),
    )
    p.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="PNG resolution (default: 300)",
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


def main() -> None:
    args = parser().parse_args()
    input_dir = Path(args.input_dir)

    summary, spectrum = load_outputs(input_dir)

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

    x = np.log(mass[valid] / m0)

    # run_phase_locked_period.py centers the selected residual field
    # before fitting the frozen waveform.
    y_raw = residual[valid]
    y = y_raw - np.mean(y_raw)

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    x_dense = np.linspace(float(np.min(x)), float(np.max(x)), 2000)
    locked_curve = c0 + amplitude * np.cos(omega * x_dense - phase)

    fig, ax = plt.subplots(figsize=(10.5, 5.8))

    ax.scatter(
        x,
        y,
        s=18,
        alpha=0.72,
        label="G2 centered Pearson residuals",
    )
    ax.plot(
        x_dense,
        locked_curve,
        linewidth=2.0,
        label="Frozen phase-locked waveform",
    )
    ax.axhline(0.0, linewidth=0.8)

    ax.set_xlabel(r"$x=\ln(m_{\mu\mu}/1\,\mathrm{GeV})$")
    ax.set_ylabel("Centered Pearson residual")
    ax.set_title(
        "CMS Run2016G G2 Prospective Phase-Locked Holdout"
    )

    annotation = (
        rf"$\omega={omega:.10f}$ (fixed)" "\n"
        rf"$\phi={phase:.10f}$ rad (fixed)" "\n"
        rf"$A={amplitude:.4f}$" "\n"
        rf"$\Delta\chi^2={delta_chi2:.2f}$" "\n"
        rf"$p_{{\rm local}}={local_p:.2e}$" "\n"
        rf"$N_{{\rm bins}}={int(summary['analysis_bins'])}$"
    )
    ax.text(
        0.015,
        0.975,
        annotation,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round", "alpha": 0.85},
    )

    ax.legend(loc="lower right")
    ax.grid(alpha=0.18)
    fig.tight_layout()

    if args.output is None:
        png_path = input_dir / "phase_locked_fit.png"
    else:
        png_path = Path(args.output)

    png_path.parent.mkdir(parents=True, exist_ok=True)
    pdf_path = png_path.with_suffix(".pdf")

    fig.savefig(png_path, dpi=args.dpi, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote: {png_path}")
    print(f"Wrote: {pdf_path}")
    print(f"Locked omega: {omega}")
    print(f"Locked phase: {phase}")
    print(f"Signed amplitude: {amplitude}")
    print(f"Delta chi-square: {delta_chi2}")
    print(f"Local one-sided p: {local_p}")


if __name__ == "__main__":
    main()
