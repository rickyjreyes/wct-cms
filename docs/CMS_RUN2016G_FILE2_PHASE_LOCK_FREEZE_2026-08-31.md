# CMS Run2016G file-2 prospective phase-locked holdout freeze — 2026-08-31

## Provenance

The previously frozen Run2016F post-VFP target could not be resolved as a public CERN Open Data Portal record. No Run2016F WCT/log-periodic spectrum was inspected.

Before inspecting any WCT result from the replacement target, the replacement rule is frozen here.

## Replacement target

Dataset:

/DoubleMuon/Run2016G-UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD

CERN Open Data record:

30522

The first file in the official record ordering,

05DD095C-F6C3-9A4F-9FB3-348A5A6403D5.root

was already used for the earlier Run2016G frozen-frequency replication and is therefore excluded.

The replacement prospective holdout is the first unused file in the official CERN ordering:

209D94D9-B6D5-A34B-A2A3-CBB7E4EA8ADF.root

Official XRootD URL:

root://eospublic.cern.ch//eos/opendata/cms/Run2016G/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v2/2430000/209D94D9-B6D5-A34B-A2A3-CBB7E4EA8ADF.root

File size:

2284233920 bytes

Adler-32:

b18e8b13

## Frozen waveform

omega = 7.025825825825827
phi   = -0.2313916852932179 rad
A > 0

The primary model remains:

r(m) = c + A cos[omega ln(m / 1 GeV) - phi]

There is no frequency refit, no phase refit, and no sign change in the primary prospective test.

## Frozen analysis

- mass range: 2--120 GeV
- bins: 350 logarithmic
- muon pT minimum: 4 GeV
- |eta| maximum: 2.4
- muon ID: tight
- maximum events read: 100000
- background: robust Chebyshev
- degree: 7
- fit iterations: 6
- clip sigma: 3.5
- minimum model count: 5
- m0: 1 GeV
- masks:
  - 2.9--3.3 GeV
  - 3.55--3.85 GeV
  - 8.5--11.5 GeV
  - 80--100 GeV
- residual permutations: 1000
- parametric refit bootstraps: 500
- seed: 20260831

## Interpretation

This is an untouched file-level prospective holdout, but it is not an independent detector or independent data-taking period. It shares the Run2016G detector, trigger, reconstruction, and period-level systematics with the previously inspected Run2016G file.

A negative signed amplitude is a prospective failure. A new optimized frequency or phase cannot replace the frozen primary result.
