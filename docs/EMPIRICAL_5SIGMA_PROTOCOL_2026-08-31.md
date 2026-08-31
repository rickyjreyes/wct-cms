# Empirical >5 sigma validation protocol — 2026-08-31

## Purpose

This document defines what this repository will require before describing a CMS result as **empirically above 5 sigma**.

The current Run2016H file-2 and Run2016G frozen-frequency replications produce local analytic fixed-frequency diagnostics near 10.5--10.6 sigma. Those values are conditional on the implemented residual/background model. The existing permutation and refit-bootstrap runs reached finite Monte Carlo floors near `10^-3`, so they do **not** yet resolve a 5-sigma empirical tail.

The next prospective CMS test is sharper: Run2016F has a separately committed freeze with frequency, phase, and positive amplitude sign fixed before inspection. If that target remains untouched under the freeze, it is the preferred prospective waveform test for deep-tail calibration.

## Statistical target

For a one-sided Gaussian convention,

\[
Z=5
\]

corresponds to

\[
p_5 = 2.866515718791946\times10^{-7}.
\]

A result is not called empirical `>5 sigma` merely because an asymptotic likelihood or chi-square conversion gives `Z > 5`. The calibrated null probability for the frozen statistic must itself be demonstrated below the threshold, subject to the controls below.

## Frozen primary statistic

For the prospective phase-locked test, the primary model is

\[
r(m)=c+A\cos[\omega\ln(m/m_0)-\phi],
\]

with all of the following fixed before the target data are inspected:

- target data selection rule;
- `omega`;
- `phi`;
- positive amplitude direction `A > 0`;
- `m0`;
- mass range;
- binning;
- muon cuts and ID;
- resonance masks;
- background family or predeclared background robustness set;
- fit degree / smoothing settings;
- clipping rule;
- minimum model-count rule;
- test statistic;
- random seeds / null-shard plan for any final deep-tail run.

The current Run2016F freeze specifies

```text
omega = 7.025825825825827
phi   = -0.2313916852932179 rad
A > 0
```

with the same `m0 = 1 GeV` convention.

A frequency or phase selected after inspecting the target data is exploratory and cannot replace this primary statistic.

## Direct Monte Carlo depth

The repository uses an add-one Monte Carlo convention,

\[
\hat p = \frac{k+1}{N+1},
\]

where `N` is the number of null pseudoexperiments and `k` is the number whose statistic equals or exceeds the observed statistic.

With zero exceedances, the smallest reportable value is `1/(N+1)`.

For the one-sided 5-sigma threshold:

| criterion | required zero-exceedance trials |
|---|---:|
| add-one numerical floor reaches `p5` | **3,488,555** |
| exact one-sided 95% binomial upper bound reaches `p5` | **10,450,778** |
| exact one-sided 99% binomial upper bound reaches `p5` | **16,065,391** |

The first number is only numerical resolution. It is not by itself a 95% confidence statement that the true null probability is below `p5`.

For a direct empirical claim, this repository uses the stronger default gate:

> **The exact one-sided 95% Clopper--Pearson upper bound on the null exceedance probability must be below `p5`.**

If zero exceedances are observed, this requires at least **10,450,778** valid null pseudoexperiments. If one or more exceedances occur, compute the exact upper bound from the observed binomial count rather than quoting the Monte Carlo floor.

Use

```bash
python scripts/plan_empirical_5sigma.py
```

to print the canonical trial counts, or evaluate a completed null run with

```bash
python scripts/plan_empirical_5sigma.py \
  --trials 10450778 \
  --exceedances 0 \
  --confidence 0.95
```

## Null hierarchy

### 1. Residual permutation

Residual permutation is retained as a fast diagnostic. It tests exchangeability of the fitted residual field and is useful for catching simple frequency structure.

It is **not sufficient by itself** for an empirical 5-sigma claim because correlated continuum mis-modeling can violate exchangeability.

### 2. End-to-end parametric refit bootstrap

The primary direct-calibration null is an end-to-end Poisson pseudoexperiment:

1. generate a pseudo-spectrum from the frozen smooth null model;
2. refit the continuum with the identical predeclared background procedure;
3. recompute residuals using the identical masks and analysis bins;
4. evaluate the exact frozen waveform statistic;
5. record whether it equals or exceeds the observed statistic.

For the phase-locked prospective test there is no need to perform an exploratory frequency scan inside the primary deep-tail toys. Frequency and phase remain fixed; the toys calibrate the one directional amplitude statistic.

### 3. Background-model robustness

A deep tail under one continuum parameterization does not prove that the continuum family is correct.

Before claiming empirical `>5 sigma`, run a predeclared robustness envelope covering reasonable continuum alternatives, including the repository's supported background families where appropriate:

- robust Chebyshev;
- robust Bernstein;
- spline background with predeclared smoothing settings;
- nearby predeclared polynomial degrees / smoothing choices;
- predeclared binning variants.

Do not select the variant that maximizes significance. The conservative interpretation is governed by the least favorable reasonable predeclared model or by another combination rule fixed before the robustness results are inspected.

## Systematic-control gate

A numerical tail below `p5` is not sufficient if the feature is unstable under basic controls. Before a discovery-grade statement, require all of the following to be reported:

1. **Prospective holdout:** the primary statistic is evaluated on data not used to choose its frequency/phase/cuts.
2. **Background stability:** the signal remains material across the predeclared continuum family/degree envelope.
3. **Binning stability:** reasonable predeclared binning changes do not create or erase the result.
4. **Era stability:** independent data-taking periods are reported separately.
5. **Simulation/control samples:** detector/reconstruction and Standard Model controls are processed through the same pipeline when suitable samples are available.
6. **Known resonance handling:** mass masks remain physics-motivated and fixed before the target result.
7. **Residual adequacy:** residual RMS, tails, and correlation structure are reported; broad non-Pearson residuals must be accounted for in the null calibration.
8. **Independent channel or detector:** a second CMS channel or independent detector is strongly preferred before a physical interpretation.
9. **Failures retained:** null and failed replications remain in the repository.
10. **No post-hoc statistic substitution:** a failed phase-locked test cannot be replaced by a newly optimized frequency, phase, background, or cut while retaining the prospective label.

## Combining replications

Do not multiply local p-values or add sigma values from Run2016H, Run2016G, and a future Run2016F result after seeing them.

If a combined significance is desired, define a **joint statistic and null procedure in advance**. The joint pseudoexperiment must reproduce the independence/dependence structure of the samples and apply the same frozen combination rule on every null draw.

Until such a joint test is preregistered and calibrated, report each replication separately.

## Tail extrapolation

A fitted exponential, generalized Pareto tail, asymptotic chi-square law, or importance-sampling estimate can be useful for computational planning. It must be labeled as modeled/extrapolated unless validated against direct Monte Carlo in an overlap region.

The cleanest claim is a direct end-to-end null count deep enough to resolve the target tail. If accelerated rare-event methods are added later, their proposal distribution, weights, diagnostics, and validation must be frozen before the final tail estimate is produced.

## Claim language

### Allowed now

```text
Run2016H file-2 and Run2016G show approximately 10.5--10.6 sigma local analytic fixed-frequency diagnostics under the implemented residual/background model. Existing empirical null runs reached Monte Carlo floors near 10^-3 and therefore do not yet establish an empirical >5 sigma significance.
```

### Allowed after the numerical floor only

If at least 3,488,555 direct null trials produce zero exceedances:

```text
The direct Monte Carlo add-one resolution reaches the one-sided 5-sigma p-value scale, with zero observed exceedances. This is a numerical-resolution statement; the stronger confidence-bound gate is not yet satisfied unless the exact upper confidence bound is also below p5.
```

### Empirical >5 sigma gate

Only when the exact one-sided 95% upper confidence bound is below `2.866515718791946e-7`, and the predeclared systematic-control gate is passed:

```text
Under the frozen statistic and stated null/systematic model, the calibrated null exceedance probability is below the one-sided 5-sigma threshold at 95% confidence.
```

This statement remains conditional on the stated experimental and background model. It is not, by itself, a claim that WCT is the physical cause of the observed structure.

## Immediate execution order

1. Preserve the existing Run2016F phase-lock freeze and commit the exact official input manifest before inspection.
2. Run the first-pass `1000` permutation / `500` refit-bootstrap prospective test exactly as frozen.
3. If the signed amplitude is negative or the prospective test fails, record the failure and stop the 5-sigma tail program for that target.
4. If the prospective result is successful and reaches the first-pass Monte Carlo floor, run the predeclared background/systematic robustness envelope.
5. Only if the result remains robust, launch deep end-to-end fixed-waveform null calibration in reproducible shards/checkpoints until either:
   - enough exceedances rule out `p < p5`, or
   - the predeclared direct-tail confidence criterion is satisfied.
6. Report the local analytic diagnostic, empirical Monte Carlo result, confidence bound, and systematics separately.
