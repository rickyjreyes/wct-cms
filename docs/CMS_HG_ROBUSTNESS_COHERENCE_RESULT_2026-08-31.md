# CMS H/G robustness and coherence result — 2026-08-31

This document records the first results from the preregistered retrospective H/G robustness/coherence program in `CMS_HG_ROBUSTNESS_COHERENCE_FREEZE_2026-08-31.md`.

These diagnostics use already-observed Run2016H file #2 and Run2016G file #1. They do not constitute a new prospective replication.

## Test integrity

The expanded test suite passed before these runs:

```text
16 passed in 1.54 s
```

The frozen references were:

```text
omega = 7.025825825825827
phi_pred = -0.2313916852932179 rad
signed alternative: A > 0
```

## H/G robustness envelope

The frozen envelope evaluated 25 configurations in each of the two periods, for 50 rows total. It included Chebyshev degrees 5--12, Bernstein and spline backgrounds, binning changes, pT/eta/ID changes, mass-domain changes, and resonance-mask changes.

First aggregate result:

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

Interpretation:

- The waveform is **not uniformly robust across the full prespecified envelope**: at least one period/configuration reverses the preregistered amplitude sign, producing zero positive-only improvement.
- The median absolute phase displacement is about 0.205 rad, so the central bulk of the 50 rows remains relatively close to the H/G phase prediction.
- The maximum phase displacement is 2.503 rad, so at least one perturbation substantially changes the waveform orientation.
- The current aggregate alone does not identify which configuration(s) cause the sign failure. `scripts/summarize_hg_robustness.py` was added after this result to report every nonpositive row and the largest phase excursions from the already-generated CSV without rerunning ROOT extraction. This summarizer does not change the frozen envelope or its outputs.

The robustness result therefore weakens any claim that the present ~10.5--10.6 sigma local analytic diagnostics are model-invariant physical significances. It does not erase the baseline H/G replication; it shows that its interpretation depends on analysis/background choices that must be localized and understood.

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
- The coherence calculation uses the baseline degree-7 continuum and therefore does not override the robustness-envelope failures.

## Combined conclusion

The two diagnostics answer different questions and must be reported together:

1. **Baseline joint coherence is strong.** The H/G common-waveform statistic had zero exceedances in 10,000 paired residual-permutation nulls.
2. **Uniform model robustness fails.** The full prespecified envelope contains at least one sign reversal and large phase excursions.

Therefore the correct current status is neither “the CMS signal disappeared” nor “CMS is already discovery-grade.” The baseline H/G waveform match remains a strong, calibrated retrospective feature under the stated residual-permutation null, but its physical interpretation is limited by sensitivity to some prespecified background/analysis choices.

The next immediate task is to identify exactly which frozen envelope configurations cause the sign reversal and largest phase drift. The prospective Run2016F fixed-omega/fixed-phase/signed-amplitude protocol remains unchanged and must stay unseen until its exact first-file manifest/checksum is committed.
