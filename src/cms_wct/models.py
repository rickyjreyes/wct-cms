from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


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
    analysis_bins: int
    residual_rms: float
    residual_max_abs: float
    log_mass_span: float
    best_cycles_across_span: float
    best_scan_at_boundary: bool
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
