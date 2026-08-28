# CMS frozen file-2 replication result — 2026-08-28

This document records the result of the preregistered CMS-to-CMS file replication defined in `CMS_CANDIDATE_FREEZE_2026-08-28.md`.

## Replication sample

- Dataset: CMS Run2016H DoubleMuon NanoAOD (`/DoubleMuon/Run2016H-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD`)
- Replication ROOT file: `183BFB78-7B5E-734F-BBF5-174A73020F89.root`
- Golden JSON: `Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt`
- Events read: 100000
- Events surviving certification: 90868
- Selected muons: 136310
- Opposite-sign dimuon pairs: 45188
- Masses in 2--120 GeV: 39230
- Analysis bins: 287

## Frozen primary test

The preregistered frequency was not changed after inspection of file #2:

```text
omega_m = 7.025825825825827
```

Frozen-frequency result:

- amplitude = `0.9367120932011154`
- delta chi-square = `118.91483403371683`
- phase = `-0.30599105090079975` rad
- local chi-square p (2 dof) = `1.506509523217018e-26` (diagnostic only)
- frozen residual-permutation p = `1/1001 = 0.000999000999000999` (0 exceedances in 1000 trials; Monte Carlo floor)
- frozen parametric refit-bootstrap p = `1/501 = 0.001996007984031936` (0 exceedances in 500 trials; Monte Carlo floor)

For comparison, discovery-file frozen-candidate values were:

- amplitude = `0.7542594046882397`
- phase = `-0.18895382226307825` rad
- delta chi-square = `75.76163387697034`

The fixed-frequency phase difference between discovery and replication is approximately `0.11704` rad (`6.71 degrees`). This phase comparison uses the same `m0 = 1 GeV` convention and therefore does not require a phase-origin adjustment.

## Secondary exploratory scan

The replication file's unrestricted exploratory scan found:

- best omega = `6.563963963963964`
- amplitude = `0.9718974045865767`
- delta chi-square = `128.55385362603215`
- phase = `-2.128625201416986` rad
- cycles across analyzed log-mass span = `4.265088926420425`
- global residual-permutation p = `1/1001`
- global parametric refit-bootstrap p = `1/501`

The exploratory optimum is secondary; the preregistered primary question was the fixed-frequency test at `omega_m = 7.025825825825827`.

## Model-adequacy warning

The replication spectrum has:

- residual RMS = `1.8186651232744468`
- maximum absolute residual = `12.468526148433225`

These are substantially larger than the ideal Pearson-residual scale and worse than the discovery file (`RMS = 1.28647`, `max |r| = 6.39864`). Therefore the smooth degree-7 continuum is not an adequate complete description of the real mass spectrum.

The small permutation and parametric-bootstrap p-values show that the observed structure is not reproduced by the implemented nulls, but they do **not** prove that the departure is WCT. The refit bootstrap is conditional on the fitted smooth continuum family; a detector, trigger, reconstruction, acceptance, or Standard-Model spectral feature absent from that null may recur in every Run2016H file.

## Status

**Result:** the preregistered fixed-frequency CMS file-to-file replication passes its implemented statistical tests.

**Scientific classification:** `CMS-internal event-split replication under shared Run2016H systematics`.

This must not be described as independent-detector confirmation or WCT confirmation. The two ROOT files contain distinct event subsets but share the same dataset, detector, run period, trigger/acceptance environment, reconstruction, and analysis model.

The next high-value test is a different data-taking period with the frequency and analysis choices frozen before inspection.