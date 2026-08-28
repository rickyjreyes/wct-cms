# CMS Run2016G frozen replication protocol — 2026-08-28

This document freezes the next cross-period CMS replication before any Run2016G DoubleMuon event data are inspected with the WCT spectral pipeline.

## Purpose

The Run2016H discovery and file-2 replication share the same detector period, trigger/acceptance environment, reconstruction, and dataset-level systematics. This test moves to a different 2016 data-taking period while preserving the previously frozen spectral candidate and analysis choices.

## Target dataset

- CERN Open Data record: `30522`
- Dataset: `/DoubleMuon/Run2016G-UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD`
- Run range: `278820--280385`
- Input selection for the first pass: the first file returned by `cernopendata-client download-files --recid 30522 --filter-range 1-1`
- Golden JSON: `Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt`

No Run2016G WCT spectral output was inspected before this protocol was committed.

## Frozen candidate

Primary fixed frequency:

```text
omega_m = 7.025825825825827
```

The coordinate remains:

```text
x = ln(m_mumu / 1 GeV)
```

## Frozen analysis configuration

- mass range: `2--120 GeV`
- bins: `350`
- logarithmic mass bins: yes
- muon pT minimum: `4 GeV`
- absolute muon eta maximum: `2.4`
- muon ID: tight
- background: robust Chebyshev degree `7` in log-count vs log-mass
- fit iterations: repository default
- clip sigma: repository default
- excluded windows:
  - `2.9--3.3 GeV`
  - `3.55--3.85 GeV`
  - `8.5--11.5 GeV`
  - `80--100 GeV`
- exploratory scan: `omega = 3.1--80`, `1000` steps
- primary fixed frequency: `7.025825825825827`
- event cap: `100000` events read before certification
- residual permutations: `1000`
- parametric refit bootstraps: `500`
- seed: `20260827`

## Primary interpretation rule

The primary Run2016G result is the fixed-frequency test:

- `frozen_scan.amplitude`
- `frozen_scan.phase`
- `frozen_scan.delta_chi2`
- `frozen_permutation_p`
- `frozen_parametric_bootstrap_p`

The unrestricted Run2016G best frequency is secondary and may not be substituted for the frozen frequency.

A small fixed-frequency result in Run2016G would strengthen CMS cross-period persistence, but would still not establish WCT causation because both periods share the CMS detector, reconstruction family, physics channel, and broad selection mechanisms.

Failure at the frozen frequency must be retained as a failed replication; the frequency may not be moved to a new Run2016G optimum and still be called this preregistered test.

## Required model-adequacy diagnostics

Interpretation must also report:

- `residual_rms`
- `residual_max_abs`
- `best_scan_at_boundary`
- `best_cycles_across_span`

Small null p-values are conditional on the implemented background/null families and must not override clear evidence of continuum-model inadequacy.