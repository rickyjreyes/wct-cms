# CMS falsification / kill-test program — 2026-09-01

This document turns the next CMS work into a ranked adversarial program. The purpose is not to accumulate additional same-model significance. The purpose is to determine whether the frozen residual survives the strongest plausible analysis-artifact and conventional-physics explanations.

## Frequency-convention lock

Two related repositories currently use different log coordinates and therefore different numerical frequency labels.

### CMS inclusive dimuon repository (`wct-cms`)

The CMS residual is analyzed in

```text
x_CMS = ln(m_mumu / m0)
omega_CMS = 7.025825825825827
```

### LHCb request-48 cross-run work (`rickyjreyes/LHC`)

The mapped LHCb analysis uses `ln(q^2)`. Since `q^2 = m^2`, the corresponding frequency is

```text
k_LHCb = 3.512912912912913
omega_CMS = 2 * k_LHCb
```

The completed LHC stage-35 deep-null record is:

```text
lhcb_wct_analysis_kit/CMS_FIXED_K_DEG6_DEEP_NULL_RESULT_2026-08-31.md
```

It reports `0/10,000` paired-null exceedances for the degree-6 cross-run phase-transfer statistic at the fixed mapped frequency, while also reporting a roughly `51.1 deg` fitted phase difference and explicitly classifying the result as post-unblinding/exploratory.

**Do not label the present CMS candidate as `k ~ 9.7`.** That number is not the frozen frequency used by the current `wct-cms` analysis or the mapped request-48 LHCb test. Any future use of a `k ~ 9.7` target must identify the source coordinate and preregistration separately.

---

## Ranked program

| Rank | Test | Core question | Current status | Repository target |
|---:|---|---|---|---|
| **1** | **CMS Background Kill Test** | Can a signal-independent flexible continuum manufacture, absorb, or select the frozen waveform? | **Central unresolved ambiguity.** The 25-configuration H/G robustness envelope showed that degree-12 polynomial and spline backgrounds can suppress/reverse the locked component while bringing residual RMS toward unity. A WCT-blind blocked-CV selector is already implemented but the spurious-signal/absorption calibration is not yet complete. | `wct-cms`: `scripts/run_hg_background_cv.py`, `src/cms_wct/background_cv.py`, plus the extension defined below. Cross-reference the LHC stage-35 degree-6 deep-null work without treating it as the same observable. |
| **2** | **CMS Freeze Test** | Does the evidence come from the genuinely frozen frequency rather than rescanning each holdout? | **Implemented, but should be made explicit as an audit.** Replication runs already report the frozen statistic separately from the exploratory scan and calibrate fixed-frequency nulls. | `wct-cms`: base pipeline, `docs/CMS_CANDIDATE_FREEZE_2026-08-28.md`, Run2016H/Run2016G result records, and the Cross-Period paper. |
| **3** | **CMS Phase Lock Test** | Does untouched data support the already-frozen frequency, phase, and positive amplitude sign? | **Prospective holdout version passed on Run2016G file-2 under the degree-7 continuum.** The locked amplitude was positive (`0.9708618`) with `Delta chi-square = 126.2832`; the finite null ensembles reached `1/1001` and `1/501` floors. This remains conditional on the continuum model and is not a flexible-background-independent significance. | `wct-cms`: `scripts/run_phase_locked_period.py`, `docs/CMS_RUN2016G_FILE2_PHASE_LOCK_FREEZE_2026-08-31.md`, `docs/CMS_RUN2016G_FILE2_PHASE_LOCK_RESULT_2026-08-31.json`, Cross-Period paper. |
| **4** | **ATLAS Precedent** | What was already done experimentally, and what is actually distinct here? | **Prior art exists and must be cited directly.** ATLAS searched periodic signals in dielectron and diphoton invariant-mass spectra with 139/fb at 13 TeV using continuous wavelet transforms and neural-network classifiers, motivated by clockwork/linear-dilaton resonance trains: JHEP 10 (2023) 079, arXiv:2305.10894, DOI `10.1007/JHEP10(2023)079`. | Cross-Period paper literature/prior-art section. |
| **5** | **Detrending Trap Test** | Can smooth fitting/subtraction itself generate a preferred log-frequency peak or phase-coherent residual? | **Not yet directly calibrated.** This should be implemented as a methodological null inside the Background Kill Test extension. | `wct-cms`: new spurious-signal null/control using the existing background-family and frozen-waveform machinery. |
| **6** | **Interference Test** | Can ordinary resonance tails, continuum interference, heavy-flavor transitions, or detector/selection structure generate the residual morphology? | **Open.** Existing resonance masks are necessary but do not model broad tails/interference or correlated detector effects. | `wct-cms` / `LHC`: collider-systematics section and dedicated control/simulation work. |

---

# 1. CMS Background Kill Test — highest priority

The present background question is no longer hypothetical. In the frozen H/G robustness envelope, most ordinary perturbations retained the positive locked component, but sufficiently flexible backgrounds did not:

```text
positive sign in both periods: 22 / 25 configurations
nonpositive rows: 4 / 50
failure configurations: cheb_d12, bernstein_d12, spline_s1
```

The spline family is especially important because it can strongly suppress/reverse the fixed waveform while producing residual RMS values close to one. Therefore a degree-7 local analytic significance cannot be treated as background-model-independent evidence.

The next extension must separate three different questions.

## 1A. Signal-independent predictive background selection

Already implemented by:

```text
scripts/run_hg_background_cv.py
src/cms_wct/background_cv.py
```

The selector uses blocked held-out Poisson deviance and does not use the WCT frequency, phase, amplitude, or `Delta chi-square` to choose the continuum.

Required output:

- winning and near-tied continuum models;
- predictive deviance in H and G separately;
- post-selection frozen-frequency and phase-locked statistics;
- explicit warning that H/G are already-observed data.

If a flexible model wins decisively on held-out prediction and removes the waveform in both periods, the degree-7 interpretation is substantially weakened.

## 1B. Spurious-signal null

For every prespecified plausible continuum family/complexity that survives predictive screening:

1. fit the smooth null without using the WCT statistic;
2. generate many pseudo-spectra from that fitted continuum;
3. rerun the **entire** detrend/refit/residual pipeline on every pseudo-spectrum;
4. evaluate both:
   - the frozen-frequency/frozen-phase directional statistic;
   - the unrestricted exploratory scan maximum;
5. record how often the pipeline creates a positive component at least as large as observed;
6. record how often the exploratory best frequency lands near the frozen frequency by chance.

This is the direct answer to: **can flexible backgrounds manufacture the signal after the same fitting/subtraction procedure used on data?**

Primary reporting should include a per-family empirical exceedance rate and a conservative envelope across the prespecified plausible background models. Do not select the null family that gives the smallest p-value.

## 1C. Signal-absorption / overfit calibration

A very flexible background can erase a real oscillation as well as correct an underfit continuum. Therefore also perform injection/recovery tests:

1. generate pseudo-spectra from each plausible smooth continuum;
2. inject the frozen waveform at several fixed amplitudes including the observed order of magnitude;
3. fit the same continuum family exactly as in the real pipeline;
4. recover the frozen amplitude/phase statistic;
5. report the retention ratio

```text
R_A = recovered amplitude / injected amplitude
```

and detection power versus injected amplitude.

A background family that wins predictive CV but systematically absorbs nearly all injected signals at this frequency cannot by itself be used to claim that disappearance of the real residual proves the signal was spurious. The identifiability/power result must be reported with the fit quality.

## 1D. Preferred final statistic

The strongest version is a **cross-fitted frozen-waveform statistic**: continuum parameters are learned from training blocks and the frozen waveform is scored only on held-out blocks. This reduces the ability of the same bins to both determine the smooth trend and create the residual being tested.

### Background-kill interpretation rule

A physical interpretation is materially weakened if a reasonable, signal-independent, predictively preferred continuum both:

- removes the frozen component in H/G; and
- produces an adequate spurious-signal calibration showing that the original degree-7 statistic is common under the improved null.

Conversely, the case strengthens if the frozen directional waveform remains unusual across the predictively competitive continuum envelope **and** injection tests show adequate sensitivity rather than systematic signal absorption.

---

# 2. CMS Freeze Test

The repository already separates the fixed-frequency replication statistic from the exploratory best-frequency scan. Make that logic impossible to miss by adding a single audit table for every holdout:

```text
sample
frozen omega
exploratory best omega
T_fixed observed
p_fixed under fixed-statistic null
T_scan observed
p_global for scan maximum
frequency chosen before holdout? yes/no
```

The scientific claim must be based on `T_fixed`, not on the holdout's rescanned maximum. The exploratory scan is only a diagnostic showing where the holdout itself would have preferred to fit.

A useful negative control is to repeat the holdout procedure for many deliberately wrong fixed frequencies selected without looking at that holdout. This estimates how exceptional the original frozen target is relative to arbitrary preselected frequencies without converting the holdout back into a scan.

---

# 3. CMS Phase Lock Test

The prospective phase-locked statistic is the sharpest existing CMS holdout design because it freezes:

```text
frequency
phase
amplitude sign A > 0
```

before opening the target file.

The Run2016G file-2 result is positive under the baseline degree-7 model. The next requirement is not another free-phase success; it is to repeat the same fixed-frequency/fixed-phase directional test under the predictively competitive background envelope from Test 1.

Do not describe phase universality as established merely because a degree-7 phase-locked holdout is strong. The separate LHC degree-6 request-48 analysis showed that fitted phase can move substantially (`~51.1 deg`) when background flexibility changes.

---

# 4. ATLAS precedent

The CMS paper must not imply that periodic invariant-mass searches are new in collider physics.

Direct precedent:

> ATLAS Collaboration, "Search for periodic signals in the dielectron and diphoton invariant mass spectra using 139 fb^-1 of pp collisions at sqrt(s)=13 TeV with the ATLAS detector," JHEP 10 (2023) 079, arXiv:2305.10894, DOI 10.1007/JHEP10(2023)079.

The correct distinction to establish is methodological/observational, not "first periodic mass-spectrum search":

- ATLAS: high-mass dielectron/diphoton spectra; semi-periodic resonance trains motivated by clockwork/linear-dilaton models; continuous wavelet transforms and neural-network classifiers.
- This repository: inclusive opposite-sign dimuon residuals over `2--120 GeV`; sinusoidal structure in `ln(m)` after smooth-continuum subtraction; fixed-frequency and fixed-phase holdout tests; end-to-end refit nulls; explicit flexible-background identifiability tests.

The paper should say that these are **different periodic-signal hypotheses and test designs**. Any stronger novelty claim requires a dedicated literature review.

---

# 5. Detrending Trap Test

This is the methodological core of Test 1B, but it should also be reported independently because it attacks the residual construction itself.

Required null classes:

- pure smooth continuum with Poisson counting noise;
- smooth continuum plus broad nonperiodic curvature/misspecification;
- alternate continuum family generating the data but a different family used for detrending;
- spectra with realistic masked gaps around resonances;
- correlated or block-structured residual perturbations when a justified model is available.

For each class, rerun the same fitting, clipping, masking, subtraction, and spectral test. Measure:

```text
P(T_fixed >= T_fixed,obs)
P(best frequency near frozen frequency)
P(positive locked amplitude and phase agreement at least as strong as observed)
```

If the detrending pipeline repeatedly concentrates power near the frozen frequency even when the generating model contains no oscillation, the residual peak is an analysis artifact until a redesigned statistic removes that behavior.

---

# 6. Interference Test

The remaining conventional-physics alternative is not simply "a resonance was missed." Broad physical components can interfere or change slope across the mass range.

The control program should test whether the frozen residual can be reproduced by reasonable combinations of:

- Drell-Yan / gamma*-Z continuum structure;
- resonance tails outside the excluded windows;
- heavy-flavor production and transitions between production regimes;
- interference terms where physically applicable;
- trigger/reconstruction/ID efficiency structure projected into dimuon mass;
- correlated detector response not represented by independent Poisson pseudo-counts.

The highest-value implementation is to pass Standard Model simulation/control samples through the **identical event selection, histogramming, continuum fit, residual construction, and frozen-waveform statistic**.

If a conventional model reproduces the frequency, phase, sign, and cross-period behavior without a new periodic component, that is a direct mundane explanation and should supersede the WCT interpretation.

---

## Claim-language lock until the kill tests are complete

The present defensible statement is:

> A frozen log-periodic residual reproduces strongly across several CMS holdouts under the declared baseline continuum, including a prospective fixed-frequency/fixed-phase/signed-amplitude holdout. However, flexible-continuum identifiability remains unresolved because some high-flexibility polynomial/spline fits suppress the waveform while improving residual whitening. The next priority is therefore a signal-independent predictive-background and spurious-signal calibration, not deeper same-model sigma.

Do not call the current local analytic `~10 sigma` diagnostics a background-model-independent or detector/systematics-calibrated discovery significance.
