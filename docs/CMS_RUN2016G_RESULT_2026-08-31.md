# CMS Run2016G frozen cross-period replication result — 2026-08-31

This document records the result of the preregistered Run2016G cross-period replication defined in `CMS_RUN2016G_REPLICATION_FREEZE_2026-08-28.md` and executed with `run_run2016g_frozen.sh`.

No Run2016G WCT spectral output was inspected before the protocol was frozen. The primary frequency and analysis configuration were not changed after inspection.

## Sample and frozen configuration

- Dataset: `/DoubleMuon/Run2016G-UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD`
- CERN Open Data record: `30522`
- ROOT file: `05DD095C-F6C3-9A4F-9FB3-348A5A6403D5.root`
- Events read: `100000`
- Events after golden JSON: `100000`
- Selected muons: `146551`
- Opposite-sign pairs: `46816`
- Masses in 2--120 GeV: `39633`
- Analysis bins: `287`
- Frozen frequency: `omega_m = 7.025825825825827`
- Mass bins: `350`, logarithmic
- Fit degree: `7`
- Residual permutations: `1000`
- Parametric refit bootstraps: `500`
- Seed: `20260827`

## Primary frozen-frequency result

At the preregistered frequency:

```text
omega_m = 7.025825825825827
amplitude = 0.9348796604929237
phase = -0.15679231968563603 rad
delta chi-square = 115.8921026429847
local chi-square p (2 dof) = 6.828882858892e-26
one-sided Gaussian-equivalent local Z = 10.4566656033 sigma
frozen residual-permutation p = 1/1001 = 0.000999000999000999
frozen parametric refit-bootstrap p = 1/501 = 0.001996007984031936
```

Both empirical null ensembles had zero exceedances and therefore reached their finite Monte Carlo floors. These floor values are not estimates that the true tail probability equals approximately `1e-3`; they only show that the implemented finite ensembles did not resolve the tail below their available resolution.

The analytic `10.4567 sigma` value is a local fixed-frequency diagnostic conditional on the implemented residual/background model. It is not, by itself, an empirical 10.46-sigma significance or a claim of 10.46-sigma physical discovery.

## Comparison with the earlier CMS samples

The frozen candidate was first identified in Run2016H file 1 and then frozen before Run2016H file 2 and before Run2016G were inspected.

| sample | frozen amplitude | frozen phase (rad) | frozen Delta chi-square | local p |
|---|---:|---:|---:|---:|
| Run2016H discovery file | `0.7542594047` | `-0.1889538223` | `75.76163388` | discovery statistic |
| Run2016H frozen file-2 replication | `0.9367120932` | `-0.3059910509` | `118.91483403` | `1.5065e-26` |
| Run2016G frozen cross-period replication | `0.9348796605` | `-0.1567923197` | `115.89210264` | `6.8289e-26` |

Notable cross-sample consistency:

- Run2016G amplitude differs from the Run2016H replication amplitude by only about `-0.20%`.
- Run2016G Delta chi-square differs from the Run2016H replication by about `-2.54%`.
- Run2016G phase differs from the original Run2016H discovery phase by `0.03216 rad` = `1.84 degrees`.
- Run2016G phase differs from the Run2016H file-2 replication phase by `0.14920 rad` = `8.55 degrees`.
- The Run2016H file-2 replication differed from the discovery phase by about `6.71 degrees`.

The phase values use the same `m0 = 1 GeV` convention, so these comparisons require no phase-origin correction.

## Secondary exploratory scan

The unrestricted Run2016G scan found:

```text
best omega = 7.4107107107107115
amplitude = 1.0316189456065652
delta chi-square = 144.1300157204339
local p = 5.041557769575671e-32
phase = 1.1455436706031235 rad
cycles across analyzed log-mass span = 4.81528240597925
```

The exploratory optimum is secondary. It may not replace the preregistered fixed-frequency result as the replication statistic.

## Model adequacy

Run2016G has:

```text
residual RMS = 1.7240310813893227
maximum absolute residual = 8.472512187270679
best scan at boundary = false
```

These diagnostics are improved relative to the Run2016H file-2 replication (`RMS = 1.8187`, `max |r| = 12.47`) but remain larger than an ideal unit-scale Pearson residual field. Therefore the smooth degree-7 continuum is still not an adequate complete physical/statistical model of the dimuon spectrum.

This is now the central limitation on converting the extremely small analytic local p-values into a discovery-grade physical significance. A smooth continuum misspecification, detector/trigger/acceptance structure, reconstruction effects, or Standard-Model spectral structure not represented in the null can recur across CMS periods.

## Status

**Primary replication verdict:** PASS under the frozen implemented tests.

**Scientific classification:** `CMS cross-period frozen-frequency replication under shared detector/reconstruction family`.

Run2016G independently reproduces the previously frozen frequency with a nearly identical fitted amplitude and a phase close to the original discovery phase. This substantially strengthens evidence that the Run2016H file-2 replication was not a one-file statistical accident.

It does **not** yet establish a globally calibrated greater-than-5-sigma physical anomaly or WCT causation. The next high-value work is a prespecified background/systematics sensitivity suite at the same frozen frequency, followed by deeper empirical-tail calibration only after adequate null models are established.
