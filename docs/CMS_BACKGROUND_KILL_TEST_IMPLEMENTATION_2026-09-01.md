# CMS Background Kill Test implementation — 2026-09-01

## Status

**Implemented and regression-tested.**

The executable is:

```text
scripts/run_cms_background_kill.py
```

The reusable implementation is:

```text
src/cms_wct/background_kill.py
```

Regression tests are in:

```text
tests/test_background_kill.py
```

The GitHub Actions test workflow passed after the implementation was committed.

## Scientific question

The test attacks the current central CMS ambiguity:

> Can signal-independent flexible-background selection and the subsequent detrending/residual pipeline manufacture the frozen waveform, or alternatively absorb a real waveform strongly enough that disappearance under a flexible fit is not informative?

The background selector is not allowed to use the WCT frequency, phase, amplitude, or Delta chi-square. It uses the existing blocked held-out Poisson-deviance cross-validation machinery.

## Frozen waveform

```text
omega = 7.025825825825827
phase = -0.2313916852932179 rad
signed alternative = A > 0
```

The default pair is:

```text
Run2016H file #2
Run2016G file #1
```

The primary pair statistic is deliberately conservative:

```text
T_primary = min(DeltaChi2_locked,H, DeltaChi2_locked,G)
```

Because each locked statistic is positive-only, this requires both periods to contribute. A large value in only one period cannot carry the primary pair score.

The secondary statistic is:

```text
T_secondary = DeltaChi2_locked,H + DeltaChi2_locked,G
```

## Components implemented

### 1. WCT-blind predictive background selection

Candidate backgrounds are ranked by combined blocked held-out Poisson deviance in H and G using the existing candidate family:

```text
Chebyshev degrees 5..12
Bernstein degrees 5, 7, 9, 12
splines with smoothing factors 0.5, 1, 2
```

Only after the predictive winner is selected are the frozen frequency and phase evaluated.

### 2. End-to-end spurious-signal / detrending null

The selected smooth background is fitted to the observed H and G spectra. Each null pseudoexperiment then:

1. draws independent Poisson pseudo-counts from the fitted smooth H and G continua;
2. reruns the WCT-blind background cross-validation unless `--no-reselect` is explicitly requested;
3. refits the selected background to each pseudo-spectrum;
4. reconstructs Pearson residuals;
5. evaluates the frozen fixed-frequency/free-phase diagnostic;
6. evaluates the fixed-frequency/fixed-phase positive-amplitude statistic;
7. records the conservative pair statistic.

The empirical add-one p-value is

```text
p = (exceedances + 1) / (trials + 1)
```

This test therefore calibrates the **complete selection + detrending + frozen-waveform pipeline**, rather than shuffling already-created residuals.

### 3. Deterministic signal-absorption matrix

The test injects the frozen waveform into the smooth background in Pearson-residual amplitude units:

```text
mu_injected = B + A * sqrt(B) * cos(omega ln(m) - phase)
```

Each prespecified background candidate is then fit to the injected expectation without Poisson noise.

For H and G it reports:

```text
recovered amplitude
retention = recovered / injected
minimum pair retention
locked pair score
```

This isolates pure continuum-fit absorption from counting noise.

### 4. End-to-end Poisson injection/recovery

For each requested injected amplitude, pseudo-spectra are generated with the frozen waveform present. Every trial reruns background selection and detrending, then records:

```text
recovered amplitude in H
recovered amplitude in G
positive sign in both
primary pair score
selected background
```

The summary reports:

```text
median recovered amplitude
median minimum retention
positive-sign-in-both fraction
power above the 95th percentile of the smooth spurious-null score
background-selection frequencies
```

This prevents a highly flexible continuum from being treated as decisive merely because it can erase the observed waveform. A continuum must also demonstrate that it would recover a waveform of the relevant size if one were actually present.

## Canonical run

From the repository root:

```bash
python scripts/run_cms_background_kill.py
```

The default run uses:

```text
100000 events per period
5 blocked CV folds
16-bin contiguous CV blocks
200 end-to-end smooth-null trials
100 injection trials per amplitude
injected amplitudes = 0.25, 0.5, 0.75, 1.0 residual units
background re-selection inside every pseudoexperiment
```

For a deeper run:

```bash
python scripts/run_cms_background_kill.py \
  --null-trials 1000 \
  --injection-trials 500 \
  --injection-amplitudes 0.25 0.5 0.75 1.0
```

`--no-reselect` is provided only as a faster diagnostic. It is not the preferred final kill test because it does not propagate background-selection variability through the pseudoexperiments.

## Outputs

The default output directory is:

```text
results/cms_background_kill/
```

It contains:

```text
selection_freeze.json
cv_scores.csv
spurious_null_trials.csv
spurious_null_summary.json
absorption_matrix.csv
injection_trials.csv
injection_summary.json
summary.json
```

`selection_freeze.json` is written before the spectrum-level kill-test calculations begin and records the candidate order, frozen waveform, pair statistic, seeds, trial counts, and claim-language lock.

## Interpretation

The Background Kill Test should be considered failed for the signal interpretation if a predictively justified smooth-background pipeline frequently produces a pair statistic at least as large as observed.

A small spurious-null p-value is not sufficient by itself. Injection/recovery must also show adequate retention and power. If the selected flexible continuum absorbs most injected waveforms at the observed scale, disappearance of the real residual under that continuum does not distinguish underfitting from overfitting.

Conversely, the CMS case is strengthened if all three conditions hold together:

1. the WCT-blind predictive selector does not trivially eliminate the frozen waveform;
2. the complete smooth-null selection/detrending pipeline rarely manufactures an equally strong locked pair response;
3. injected waveforms of the observed order are recovered with substantial retention and useful detection power.

## Scope lock

This implementation is **retrospective** because the H/G spectra have already been inspected. It is a model-identifiability and analysis-artifact test, not a new prospective replication.

It does not yet replace detector, trigger, reconstruction, acceptance, correlated-systematic, Standard Model simulation, or interference controls. Those remain separate requirements before a physical discovery claim.
