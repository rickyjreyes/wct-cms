# CMS Run2016G file-2 prospective phase-locked holdout freeze — 2026-08-31

The previously frozen Run2016F post-VFP target could not be resolved as a public CERN Open Data Portal record. No Run2016F WCT/log-periodic spectrum was inspected.

Replacement dataset:
/DoubleMuon/Run2016G-UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD

CERN Open Data record:
30522

Previously used file:
05DD095C-F6C3-9A4F-9FB3-348A5A6403D5.root

Prospective holdout file:
209D94D9-B6D5-A34B-A2A3-CBB7E4EA8ADF.root

Size:
2284233920 bytes

Adler-32:
b18e8b13

Frozen waveform:

omega = 7.025825825825827
phi   = -0.2313916852932179 rad
A > 0

Primary model:

r(m) = c + A cos[omega ln(m / 1 GeV) - phi]

No frequency refit, phase refit, or sign change is allowed in the primary test.

Frozen analysis:
- 2--120 GeV
- 350 logarithmic bins
- muon pT >= 4 GeV
- |eta| <= 2.4
- tight ID
- maximum 100000 events
- robust Chebyshev degree 7
- masks 2.9--3.3, 3.55--3.85, 8.5--11.5, 80--100 GeV
- 1000 permutations
- 500 parametric refit bootstraps
- seed 20260831

Execution note:
An attempted execution before this document was committed terminated during XRootD filesystem initialization because fsspec_xrootd was not installed. No ROOT input file was opened and no WCT result was produced or inspected.

This is a file-level holdout within Run2016G, not an independent detector or independent run-period replication.
