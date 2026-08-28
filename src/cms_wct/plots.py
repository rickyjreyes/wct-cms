from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np


def write_plots(outdir, centers, counts, model, residuals, mask, omegas, scores, best, frozen, x, y, null_max, log_bins=False):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.step(centers, counts, where="mid", label="CMS dimuon data")
    ax.plot(centers, model, label="Smooth background")
    ax.set_xlabel("Dimuon invariant mass [GeV]")
    ax.set_ylabel("Counts / bin")
    ax.set_yscale("log")
    if log_bins:
        ax.set_xscale("log")
    ax.legend()
    fig.tight_layout(); fig.savefig(outdir / "spectrum_background.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.axhline(0.0, linewidth=1)
    ax.plot(centers[mask], residuals[mask], marker=".", linestyle="none")
    ax.set_xlabel("Dimuon invariant mass [GeV]")
    ax.set_ylabel("Pearson residual")
    if log_bins:
        ax.set_xscale("log")
    fig.tight_layout(); fig.savefig(outdir / "residuals.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(omegas, scores)
    ax.axvline(best.omega, linestyle="--", label=f"best={best.omega:.5g}")
    if frozen is not None:
        ax.axvline(frozen.omega, linestyle=":", label=f"frozen={frozen.omega:.5g}")
    ax.set_xlabel("omega in log(m/m0)"); ax.set_ylabel("Delta chi-square"); ax.legend()
    fig.tight_layout(); fig.savefig(outdir / "omega_scan.png", dpi=180); plt.close(fig)

    order = np.argsort(x)
    fit = best.c0 + best.cos_coeff * np.cos(best.omega * x) + best.sin_coeff * np.sin(best.omega * x)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(x, y, marker=".", linestyle="none", label="Residuals")
    ax.plot(x[order], fit[order], label="Best sinusoid")
    ax.set_xlabel("log(m/m0)"); ax.set_ylabel("Residual"); ax.legend()
    fig.tight_layout(); fig.savefig(outdir / "best_log_periodic_fit.png", dpi=180); plt.close(fig)

    if len(null_max):
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.hist(null_max, bins=40)
        ax.axvline(best.delta_chi2, linestyle="--", label="Observed max")
        ax.set_xlabel("Permutation maximum Delta chi-square"); ax.set_ylabel("Permutations"); ax.legend()
        fig.tight_layout(); fig.savefig(outdir / "global_null_distribution.png", dpi=180); plt.close(fig)
