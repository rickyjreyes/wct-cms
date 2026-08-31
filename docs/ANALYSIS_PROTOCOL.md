# Analysis protocol

## Primary question

Does a signature specified before inspection of CMS data appear in an independent CMS observable?

The primary model tested by this repository is

\[
r(m)=c+A\cos\left(\omega\ln(m/m_0)-\phi\right)+\epsilon,
\]

where `omega` should be frozen from the external WCT analysis before the CMS result is inspected.

## Primary vs exploratory result

**Primary replication test:** `--frozen-omega`.

The frequency is not fit to CMS. Only amplitude, phase, and intercept are fit. The repository reports the fixed-frequency improvement and null p-values.

For a sharper prospective replication, `scripts/run_phase_locked_period.py` supports a **fixed-frequency, fixed-phase, signed-amplitude** statistic. Under that protocol, `omega`, `phi`, and the directional alternative `A > 0` are frozen before the target data are inspected, leaving only the amplitude/intercept fit in the primary waveform test.

**Exploratory diagnostic:** omega scan.

The maximum scan statistic is corrected using null distributions of the maximum statistic. It must not be reported as though it were the pre-registered replication test.

## Cross-dataset frequency conventions

This CMS repository currently scans the coordinate

\[
x_{\rm CMS}=\ln(m_{\mu\mu}/m_0).
\]

The existing LHCb WCT analysis uses

\[
\ell_{\rm LHCb}=\ln(q^2),
\]

with `q2 = m_mumu^2`. Therefore, apart from an additive reference-scale constant absorbed by phase,

\[
\ell_{\rm LHCb}=2x_{\rm CMS}+\text{constant}.
\]

Consequently the raw frequency conventions obey

\[
\omega_{\rm CMS}=2k_{\rm LHCb},
\qquad
k_{\rm LHCb}=\omega_{\rm CMS}/2.
\]

Never compare the numerical CMS `omega` directly to an LHCb `k` without this conversion. More generally, a frozen cross-dataset claim must first reproduce the exact dimensionless coordinate and frequency convention used by the source analysis.

## Spectral-resolution safeguards

An exploratory frequency should not be interpreted as a resolved oscillation merely because it maximizes the scan statistic.

For an analyzed coordinate span

\[
\Delta x = \max\ln(m/m_0)-\min\ln(m/m_0),
\]

the number of cycles represented by a candidate frequency is

\[
N_{\rm cycles}=\frac{\omega\,\Delta x}{2\pi}.
\]

A candidate spanning less than one complete cycle is effectively a broad trend and is strongly degenerate with the background model. Fewer than two cycles should be treated cautiously. A best-fit frequency that lands on either scan boundary is unresolved and must not be reported as an interior spectral peak.

For the default 2--120 GeV interval, the full log-mass span is approximately 4.09. One complete cycle therefore requires roughly `omega >= 1.53`, and two complete cycles require roughly `omega >= 3.07` before masking changes the effective span.

## Explicit mass masks

`--mask-window low:high` excludes the specified mass interval from both the smooth-background fit and the spectral analysis. This is important for known narrow or broad Standard Model structures. Excluding a resonance only from the background fit while retaining its residual in the spectral scan would artificially create a large structured residual.

Mass windows must be selected for physics reasons before inspecting which choice maximizes the WCT statistic. Keep a record of every mask set tested.

## Null models

### Residual-permutation null

The residual-permutation null is a fast exploratory diagnostic. It shuffles the fitted residuals over the fixed log-mass coordinates and records the strongest frequency in each permutation.

With `N` permutations and the repository's add-one correction, the smallest reportable Monte Carlo p-value is

\[
p_{\min}=\frac{1}{N+1}.
\]

Thus 200 permutations can report no value below `1/201 = 0.004975...`; that value means only that zero of 200 shuffled residual sets exceeded the observed statistic. It is not a measurement of a true p-value equal to approximately 0.005.

Residual permutation assumes the residuals are exchangeable after background fitting. Correlated background mis-modeling can violate that assumption and make the permutation p-value overly optimistic.

### Refit parametric-background bootstrap

`--parametric-bootstrap N` runs a stronger null for the current continuum model:

1. Generate a Poisson pseudo-spectrum from the fitted smooth background.
2. Refit the smooth background with the identical degree, clipping, and resonance exclusions.
3. Recompute Pearson-like residuals on the fixed analysis bins.
4. Rescan the full exploratory omega range.
5. Record the maximum scan statistic, and the frozen-frequency statistic when supplied.

This propagates counting noise and background-fit re-estimation into the null. It does **not** prove that the chosen background family is correct; model-family and fit-degree stability tests remain required.

A candidate that is significant only under residual permutation but not under the refit bootstrap should be treated as a background-fitting artifact until demonstrated otherwise.

For the prospective fixed-frequency/fixed-phase test, the primary deep-tail bootstrap should evaluate the **same frozen waveform statistic on every pseudoexperiment**. There is no need to rerun an exploratory omega scan inside that primary null because the frequency and phase are already fixed.

## Empirical >5 sigma gate

A one-sided Gaussian `5 sigma` threshold corresponds to

\[
p_5 = 2.866515718791946\times10^{-7}.
\]

An asymptotic/local analytic conversion above 5 sigma is not by itself an empirical 5-sigma result. The calibrated null distribution of the frozen statistic must reach this probability scale.

Using the repository's add-one Monte Carlo convention, zero exceedances require at least:

| requirement | trials |
|---|---:|
| numerical floor `1/(N+1) <= p5` | **3,488,555** |
| exact one-sided 95% upper confidence bound `<= p5` | **10,450,778** |
| exact one-sided 99% upper confidence bound `<= p5` | **16,065,391** |

The repository's default discovery-grade direct-Monte-Carlo gate is the **95% exact upper-bound criterion**, not merely hitting the numerical floor. If one or more null trials exceed the observed statistic, use the exact binomial upper bound for that observed exceedance count.

The helper

```bash
python scripts/plan_empirical_5sigma.py
```

prints these requirements and can evaluate a completed null count.

Deep-tail calibration is only meaningful after the prospective statistic and systematic-analysis envelope are frozen. Do not change frequency, phase, cuts, masks, background model, or binning after seeing deep-tail results and retain the prospective label.

See `docs/EMPIRICAL_5SIGMA_PROTOCOL_2026-08-31.md` for the full execution and claim-language rules.

## Controls required before a physics claim

1. Repeat with multiple smooth-background families/degrees.
2. Vary binning without choosing the variation that maximizes significance.
3. Repeat across data-taking eras.
4. Run CMS simulation/background samples through the identical pipeline.
5. Test known detector/reconstruction resonances and predeclared masked windows.
6. Require an interior frequency with enough cycles to be spectrally resolvable.
7. Require consistency under the refit parametric-background bootstrap, not only residual permutation.
8. Test at least one independent CMS channel, preferably photons after dimuons.
9. Apply the same dimensionless mapping used to compare GWTC, LHC, JUNO and photodiode data.
10. Record all failed as well as successful runs.
11. For a claimed empirical `>5 sigma` result, require the predeclared robustness envelope and exact tail-calibration criterion in `docs/EMPIRICAL_5SIGMA_PROTOCOL_2026-08-31.md` to pass.

A CMS match by itself is not evidence that WCT caused the structure. The intended value is independent cross-validation under a frozen prediction and systematic controls.
