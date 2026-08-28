# CMS exploratory candidate freeze — 2026-08-28

This document freezes the first CMS-discovered candidate before inspecting an independent CMS file. It is a preregistration for a CMS-to-CMS replication test, not evidence that WCT predicted the frequency in advance.

## Discovery sample

- Dataset: CMS Run2016H DoubleMuon NanoAOD (`/DoubleMuon/Run2016H-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD`)
- Discovery ROOT file: `127C2975-1B1C-A046-AABF-62B77E757A86.root`
- Golden JSON: `Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt`
- Events read: 100000
- Events surviving certification: 50951
- Opposite-sign dimuon pairs: 25275
- Masses in 2--120 GeV: 21930

## Frozen analysis configuration for independent CMS replication

- Observable: opposite-sign dimuon invariant mass `m_mumu`
- Spectral coordinate: `x = ln(m_mumu / 1 GeV)`
- Mass range: 2--120 GeV
- Binning: 350 logarithmic bins
- Muon selection: `pT >= 4 GeV`, `|eta| <= 2.4`, tight ID
- Background: robust Chebyshev degree 7 in log-count vs log-mass space
- Fit iterations: repository default
- Clip sigma: repository default
- Explicit excluded windows:
  - 2.9--3.3 GeV
  - 3.55--3.85 GeV
  - 8.5--11.5 GeV
  - 80--100 GeV
- Exploratory scan range used in discovery: omega = 3.1--80 with 1000 steps
- Event cap for first replication pass: 100000 events read before certification
- Residual permutation count for replication: 1000
- Parametric refit-bootstrap count for replication: 500
- Seed: 20260827

## Frozen candidate

The discovery scan produced:

- `omega_m = 7.025825825825827`
- amplitude = `0.7542594046882397`
- delta chi-square = `75.76163387697034`
- phase = `-0.18895382226307825`
- cycles across analyzed log-mass span = `4.565194460725239`
- residual RMS = `1.2864677547461774`
- maximum absolute residual = `6.398643818737002`
- residual-permutation global p = `1/1001 = 0.000999000999000999` (0 exceedances in 1000 trials; Monte Carlo floor)
- parametric refit-bootstrap global p = `1/501 = 0.001996007984031936` (0 exceedances in 500 trials; Monte Carlo floor)

For the independent CMS replication, the primary fixed-frequency test is therefore:

```text
--frozen-omega 7.025825825825827
```

No change to the frozen frequency, mass range, binning, muon cuts, background degree, exclusion windows, event cap, or null counts may be made after inspecting the replication result and still be called the preregistered replication test. Any altered configuration must be labeled a separate sensitivity or exploratory analysis.

## Interpretation rule

The independent CMS file tests whether the fixed candidate reproduces at the predeclared frequency. The exploratory best frequency from the replication file is secondary. A small fixed-frequency p-value is relevant to CMS-internal replication; it is not by itself a WCT confirmation.

## Coordinate-convention warning

The LHCb repository uses `ell = ln(q^2)` while this CMS analysis uses `ln(m_mumu)`. Because `q^2 = m_mumu^2`, the frequency conventions satisfy

```text
omega_CMS = 2 * k_LHCb
k_LHCb = omega_CMS / 2
```

Thus the frozen CMS candidate corresponds to `k_LHCb-equivalent = 3.5129129129129135`; a raw numerical comparison of 7.0258 to an LHCb `k = 7.0258` is not valid.
