from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import awkward as ak
import numpy as np
import uproot


def expand_inputs(raw_inputs):
    if len(raw_inputs) == 1:
        p = Path(raw_inputs[0])
        if p.exists() and p.suffix.lower() in {".txt", ".list", ".manifest"}:
            rows = [x.strip() for x in p.read_text(encoding="utf-8").splitlines()]
            return [x for x in rows if x and not x.startswith("#")]
    return list(raw_inputs)


def load_golden_json(path: Optional[str]):
    if path is None:
        return None
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    return {
        int(run): [(int(lo), int(hi)) for lo, hi in ranges]
        for run, ranges in raw.items()
    }


def golden_mask(runs, lumis, golden):
    if golden is None:
        return np.ones(len(runs), dtype=bool)
    out = np.zeros(len(runs), dtype=bool)
    for run in np.unique(runs):
        ranges = golden.get(int(run))
        if not ranges:
            continue
        idx = runs == run
        local_lumis = lumis[idx]
        keep = np.zeros(len(local_lumis), dtype=bool)
        for lo, hi in ranges:
            keep |= (local_lumis >= lo) & (local_lumis <= hi)
        out[idx] = keep
    return out


def invariant_mass_from_pairs(pairs):
    def comp(mu):
        pt, eta, phi, mass = mu.pt, mu.eta, mu.phi, mu.mass
        px = pt * np.cos(phi)
        py = pt * np.sin(phi)
        pz = pt * np.sinh(eta)
        e = np.sqrt(np.maximum(px * px + py * py + pz * pz + mass * mass, 0.0))
        return px, py, pz, e

    px1, py1, pz1, e1 = comp(pairs.a)
    px2, py2, pz2, e2 = comp(pairs.b)
    m2 = (e1 + e2) ** 2 - (px1 + px2) ** 2 - (py1 + py2) ** 2 - (pz1 + pz2) ** 2
    return np.sqrt(np.maximum(m2, 0.0))


def extract_dimuon_masses(inputs, args, golden=None):
    branches = [
        "run", "luminosityBlock", "Muon_pt", "Muon_eta", "Muon_phi",
        "Muon_mass", "Muon_charge",
    ]
    if args.tight_id:
        branches.append("Muon_tightId")
    if args.medium_id:
        branches.append("Muon_mediumId")

    counters = {
        "events_read": 0,
        "events_after_json": 0,
        "selected_muons": 0,
        "opposite_sign_pairs": 0,
    }
    masses = []
    sources = [f"{src}:{args.tree}" for src in inputs]

    for arrays in uproot.iterate(
        sources, expressions=branches, step_size=args.step_size, library="ak", how=dict
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
        gm = golden_mask(
            ak.to_numpy(arrays["run"]), ak.to_numpy(arrays["luminosityBlock"]), golden
        )
        counters["events_after_json"] += int(np.count_nonzero(gm))
        gm = ak.Array(gm)

        pt, eta = arrays["Muon_pt"][gm], arrays["Muon_eta"][gm]
        keep = (pt >= args.muon_pt_min) & (abs(eta) <= args.muon_eta_max)
        if args.tight_id:
            keep &= arrays["Muon_tightId"][gm]
        if args.medium_id:
            keep &= arrays["Muon_mediumId"][gm]

        mu = ak.zip({
            "pt": pt[keep],
            "eta": eta[keep],
            "phi": arrays["Muon_phi"][gm][keep],
            "mass": arrays["Muon_mass"][gm][keep],
            "charge": arrays["Muon_charge"][gm][keep],
        })
        counters["selected_muons"] += int(ak.sum(ak.num(mu, axis=1)))
        pairs = ak.combinations(mu, 2, fields=["a", "b"])
        pairs = pairs[(pairs.a.charge * pairs.b.charge) < 0]
        counters["opposite_sign_pairs"] += int(ak.sum(ak.num(pairs, axis=1)))
        values = ak.to_numpy(ak.flatten(invariant_mass_from_pairs(pairs), axis=None))
        values = values[np.isfinite(values)]
        if values.size:
            masses.append(values.astype(float, copy=False))

    if not masses:
        raise RuntimeError("No opposite-sign dimuon masses survived selection")
    return np.concatenate(masses), counters
