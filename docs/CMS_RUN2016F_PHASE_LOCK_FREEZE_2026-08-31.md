# CMS Run2016F phase-locked prospective replication freeze — 2026-08-31

This document freezes the next CMS dimuon replication **before any Run2016F WCT spectral output is inspected**.

## Purpose

Run2016H file #2 and Run2016G file #1 both reproduced the previously frozen frequency with closely matching positive amplitudes and phases. The next test is deliberately sharper: frequency, phase, and amplitude sign are fixed before the new period is analyzed.

This test is prospective only. The phase below is derived from already-observed H2/G replications and therefore does not increase their retrospective significance.

## Target data

Target dataset family:

```text
/DoubleMuon/Run2016F-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD
```

Use the **post-VFP Run2016F** dataset, not the HIPM/pre-VFP Run2016F dataset, so the target is in the same post-VFP processing family as Run2016G/H.

Frozen input-selection rule:

1. Resolve the official CERN Open Data record corresponding to the exact dataset string above.
2. Use the first ROOT file in the official record's file ordering.
3. Record the record ID, ROOT filename, and checksum in the repository before running `scripts/run_phase_locked_period.py`.
4. Do not inspect a WCT/log-periodic spectrum from Run2016F before those identifiers are committed.

Golden JSON remains:

```text
data/cert/14221/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt
```

## Frozen waveform

Frequency:

```text
omega = 7.025825825825827
```

The phase prediction is the circular mean of the two genuine frozen replications:

```text
Run2016H file #2 phase = -0.30599105090079975 rad
Run2016G file #1 phase = -0.15679231968563603 rad
```

which gives

```text
phi_pred = -0.2313916852932179 rad
```

The exact prospective waveform is

```text
g(m) = cos[7.025825825825827 * ln(m / 1 GeV) - (-0.2313916852932179)]
```

or equivalently

```text
g(m) = cos[7.025825825825827 * ln(m / 1 GeV) + 0.2313916852932179].
```

The tested model is only

```text
r(m) = c + A g(m)
```

with the directional alternative

```text
A > 0.
```

There is **no phase refit**, no sine/cosine two-parameter fit, and no frequency scan in the primary prospective test.

## Frozen analysis configuration

- mass range: `2--120 GeV`
- bins: `350`
- logarithmic mass bins: yes
- muon pT minimum: `4 GeV`
- absolute muon eta maximum: `2.4`
- muon ID: tight
- event cap: `100000` events read before certification
- background family: robust Chebyshev
- background degree: `7`
- fit iterations: repository default `6`
- clip sigma: repository default `3.5`
- minimum model count: repository default `5`
- `m0 = 1 GeV`
- excluded windows:
  - `2.9--3.3 GeV`
  - `3.55--3.85 GeV`
  - `8.5--11.5 GeV`
  - `80--100 GeV`
- residual permutations: `1000` for the first frozen pass
- parametric refit bootstraps: `500` for the first frozen pass
- seed: `20260831`

The empirical null counts are intentionally the same order as the prior frozen replications for the first pass. If the result reaches a Monte Carlo floor, deeper tail work must reuse the same frozen statistic rather than change the waveform or cuts.

## Primary result fields

The primary result is:

- `observed.signed_amplitude`
- `observed.delta_chi2`
- `observed.local_p_one_sided_chi2_1dof`
- `permutation_p`
- `parametric_bootstrap_p`
- `residual_rms`
- `residual_max_abs`

A negative fitted amplitude is a prospective failure: the positive-only statistic is set to zero and the directional local p-value is 1.

## Canonical command template

After the Run2016F manifest has been committed, run:

```bash
python scripts/run_phase_locked_period.py \
  --input data/files_run2016f.txt \
  --golden-json data/cert/14221/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt \
  --output-dir results/run2016f_phase_locked \
  --omega 7.025825825825827 \
  --phase -0.2313916852932179 \
  --mass-min 2 \
  --mass-max 120 \
  --bins 350 \
  --log-bins \
  --muon-pt-min 4 \
  --muon-eta-max 2.4 \
  --tight-id \
  --max-events 100000 \
  --background-family chebyshev \
  --fit-degree 7 \
  --mask-window 2.9:3.3 \
  --mask-window 3.55:3.85 \
  --mask-window 8.5:11.5 \
  --mask-window 80:100 \
  --permutations 1000 \
  --parametric-bootstrap 500 \
  --seed 20260831
```

## Interpretation lock

- Positive amplitude with a small calibrated p-value is a successful prospective waveform replication.
- Negative amplitude is a failure of the signed prediction.
- A newly preferred frequency or phase from an exploratory Run2016F scan may be reported separately, but it cannot replace this primary result.
- Background/systematics robustness is evaluated separately and may reduce the scientific strength of a nominally successful frozen test.
- This remains CMS-internal replication and does not by itself establish WCT causation.
