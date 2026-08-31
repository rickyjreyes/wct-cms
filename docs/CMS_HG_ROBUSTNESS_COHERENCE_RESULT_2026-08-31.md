# CMS H/G robustness and coherence result — 2026-08-31

This document records the first results from the preregistered retrospective H/G robustness/coherence program in `CMS_HG_ROBUSTNESS_COHERENCE_FREEZE_2026-08-31.md`.

These diagnostics use already-observed Run2016H file #2 and Run2016G file #1. They do not constitute a new prospective replication.

## Test integrity

The expanded test suite passed before these runs, and after the vectorized coherence-null parity test was added the repository reported:

```text
17 passed in 1.52 s
```

The frozen references were:

```text
omega = 7.025825825825827
phi_pred = -0.2313916852932179 rad
signed alternative: A > 0
```

## H/G robustness envelope

The frozen envelope evaluated 25 configurations in each of the two periods, for 50 rows total. It included Chebyshev degrees 5--12, Bernstein and spline backgrounds, binning changes, pT/eta/ID changes, mass-domain changes, and resonance-mask changes.

Aggregate result:

```text
n_rows = 50
n_configs = 25
n_unique_event_selections = 10
all_positive_sign = false
min_locked_signed_amplitude = -0.04327409947458103
min_locked_delta_chi2 = 0.0
max_abs_phase_offset_rad = 2.5029746693957438
median_abs_phase_offset_rad = 0.20496274120362035
```

The detailed summarizer localized the failures:

```text
configs positive in both periods = 22 / 25
nonpositive rows = 4 / 50
|phase offset| > 0.5 rad = 9 / 50
|phase offset| > 1.0 rad = 7 / 50
```

### Configurations preserving the positive sign in both periods

```text
baseline
bernstein_d5
bernstein_d7
bernstein_d9
bins_300
bins_400
cheb_d5
cheb_d6
cheb_d8
cheb_d9
cheb_d10
cheb_d11
eta_2p1
masks_narrow
masks_wide
mass_2_to_70
mass_3_to_120
medium_id
pt_5
pt_6
spline_s0.5
spline_s2
```

Thus every prespecified binning, pT, eta, ID, mass-domain, and resonance-mask perturbation retained a positive phase-locked coefficient in both H and G. The positive sign also survived Chebyshev degrees 5--11 and Bernstein degrees 5--9.

### Configurations with a nonpositive period

Only three configuration names produced any sign failure:

```text
cheb_d12
bernstein_d12
spline_s1
```

The four failing rows were:

| period | configuration | signed A | fixed-omega Delta chi-square | phase offset (rad) | residual RMS | max |r| |
|---|---|---:|---:|---:|---:|---:|
| H2 | cheb_d12 | -0.0432741 | 1.23562 | -2.11507 | 1.05372 | 3.17268 |
| H2 | bernstein_d12 | -0.0432741 | 1.23562 | -2.11507 | 1.05372 | 3.17268 |
| H2 | spline_s1 | -0.0140454 | 0.163145 | -2.06337 | 1.01248 | 3.34666 |
| G | spline_s1 | -0.00510813 | 0.00588828 | -2.50297 | 1.01610 | 3.18954 |

`cheb_d12` and `bernstein_d12` produce numerically identical H2 results to the displayed precision. Degree-12 Chebyshev and Bernstein polynomials span the same degree-12 polynomial function space; under the same weighting/clipping solution these two rows should therefore not be treated as two independent physical failure mechanisms. They are best viewed as one high-flexibility degree-12 polynomial failure expressed in two bases.

The spline family is more consequential. `spline_s1` suppresses/reverses the fixed waveform in both periods, while `spline_s0.5` leaves only a nearly zero locked component in both periods. These flexible backgrounds also bring the residual RMS much closer to unity than the baseline degree-7 model.

### Interpretation of robustness result

The detailed result is stronger than the aggregate `all_positive_sign = false` flag alone suggests, because the candidate waveform survives all prespecified detector/selection/domain/mask perturbations and most polynomial-background choices. It is not a fragile one-cut or one-bin artifact.

However, the result also identifies the current central ambiguity: sufficiently flexible smooth continuum models can absorb most or all of the fixed-omega structure while producing residual fields near unit RMS. Therefore the present ~10.5--10.6 sigma degree-7 analytic diagnostics cannot be interpreted as background-model-independent physical significances.

Residual RMS alone cannot decide the issue, because increasing background flexibility can both correct genuine continuum misspecification and overfit/remove a real oscillatory component. The next background test must therefore use an objective predictive model-selection rule fixed independently of the signal statistic (for example held-out-bin/cross-validated predictive likelihood or a similarly prespecified complexity penalty), rather than choosing the background that leaves the largest or smallest WCT residual.

## H/G common-waveform coherence

At the frozen frequency, the observed baseline H/G coefficient vectors were:

```text
Run2016H amplitude = 0.9367120932011158
Run2016H phase     = -0.30599105090079987 rad
Run2016G amplitude = 0.9348796604929229
Run2016G phase     = -0.15679231968563617 rad
phase difference   = -0.1491987312151637 rad
```

The frozen coherence statistic returned:

```text
common amplitude = 0.9331932040821421
common phase = -0.23146485961723295 rad
common score = 233.35048747975017
separate score = 234.80693667670135
heterogeneity = 1.4564491969511835
T_coh = 231.894038282799
```

With 10,000 complete paired residual-permutation null trials:

```text
empirical p = 1 / 10001 = 9.999000099990002e-05
```

No null pair reached the observed coherence statistic, so this value is the Monte Carlo resolution floor of the run, not a resolved tail estimate below 1e-4.

Interpretation:

- Under this specific paired residual-permutation null, the **baseline H/G joint waveform agreement is highly unusual at the available 10,000-trial resolution**.
- The small heterogeneity relative to the common score quantitatively captures the close amplitude/phase match rather than multiplying the two individual analytic p-values.
- This is still retrospective: H and G were observed before `T_coh` was introduced.
- The coherence calculation uses the baseline degree-7 continuum and therefore does not override the flexible-background robustness failures.

## Combined conclusion

The two diagnostics answer different questions and must be reported together:

1. **Baseline joint coherence is strong.** The H/G common-waveform statistic had zero exceedances in 10,000 paired residual-permutation nulls.
2. **Most frozen perturbations retain the waveform sign.** Twenty-two of 25 configurations retain positive locked amplitude in both periods, including every tested binning, pT, eta, ID, mass-domain, and resonance-mask perturbation and polynomial degrees through 11.
3. **Flexible-continuum robustness is not established.** Degree-12 polynomial backgrounds suppress/reverse H2, and the spline family can strongly suppress or reverse the fixed waveform in both periods while improving residual RMS toward unity.

Therefore the correct current status is neither “the CMS signal disappeared” nor “CMS is already discovery-grade.” The baseline H/G waveform match is a strong, calibrated retrospective feature under the stated residual-permutation null and is robust to a broad set of ordinary analysis perturbations. Its physical interpretation remains limited by identifiability against sufficiently flexible smooth continuum models.

The next retrospective task is an objective, signal-independent comparison of continuum families using predictive performance/complexity control. The prospective Run2016F fixed-omega/fixed-phase/signed-amplitude protocol remains unchanged and must stay unseen until its exact first-file manifest/checksum is committed.
