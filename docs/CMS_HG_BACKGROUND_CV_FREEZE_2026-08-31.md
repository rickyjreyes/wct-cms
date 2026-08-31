# CMS H/G WCT-blind background selection freeze — 2026-08-31

This document freezes an objective continuum-model selection test before its output is inspected.

Run2016H file #2 and Run2016G file #1 are already-observed data, so this is a **retrospective model-selection diagnostic**, not a new prospective replication. Its purpose is narrower: determine which background family/complexity best predicts held-out spectrum counts **without using the WCT frequency, phase, amplitude, or delta-chi-square in model selection**.

## Question

The baseline robust Chebyshev degree-7 continuum leaves a strong common H/G waveform near

```text
omega = 7.025825825825827
```

but more flexible backgrounds can absorb that structure. The relevant test is therefore:

> When background complexity is chosen only by held-out predictive performance, what continuum wins, and does the already-frozen waveform remain after that WCT-blind choice?

## Frozen data and baseline event selection

Both periods use the same established baseline event selection:

- Run2016H file #2 manifest: `data/files_replication.txt`
- Run2016G file #1 manifest: `data/files_run2016g.txt`
- golden JSON: `data/cert/14221/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt`
- event cap: `100000`
- opposite-sign dimuons
- muon `pT >= 4 GeV`
- `|eta| <= 2.4`
- tight ID
- mass range: `2--120 GeV`
- `350` logarithmic bins
- excluded windows:
  - `2.9--3.3 GeV`
  - `3.55--3.85 GeV`
  - `8.5--11.5 GeV`
  - `80--100 GeV`

No event-selection, mass-range, binning, or mask variant is selected in this test. Those remain part of the separate robustness envelope.

## Frozen continuum candidates

Candidate order is fixed as:

```text
cheb_d5
cheb_d6
cheb_d7
cheb_d8
cheb_d9
cheb_d10
cheb_d11
cheb_d12
bernstein_d5
bernstein_d7
bernstein_d9
bernstein_d12
spline_s0.5
spline_s1
spline_s2
```

Chebyshev and Bernstein parameterizations of the same polynomial degree can span the same polynomial function space; any numerical ties are retained rather than treated as independent physical evidence.

## Frozen cross-validation geometry

Cross-validation is deterministic and blocked:

```text
n_folds = 5
block_size = 16 eligible histogram bins
```

Eligible bins are the mass bins outside the frozen resonance masks. Within each contiguous eligible mass region, bins are grouped into contiguous blocks of up to 16 bins. Blocks are assigned cyclically to folds `0,1,2,3,4` and never bridge an excluded resonance gap.

For each fold:

1. fit the candidate background using training bins only;
2. robust clipping, if any, is restricted to the training mask;
3. predict the held-out validation blocks;
4. score the held-out counts with Poisson deviance.

The validation counts cannot alter the training coefficients or clipping mask.

## Frozen selection statistic

For each candidate, sum held-out Poisson deviance over all five folds in H and all five folds in G:

```text
D_combined = D_H + D_G
```

The primary selector is the candidate with the minimum `D_combined`.

For reporting, also record:

- H deviance per held-out bin;
- G deviance per held-out bin;
- combined deviance per held-out bin;
- top-five candidate ranking;
- numerical ties within the prespecified floating tolerance.

No WCT quantity is used for ranking.

## Post-selection waveform evaluation

Only after the CV winner has been selected, fit that winner on the complete baseline spectrum and evaluate the already-frozen waveform:

```text
omega = 7.025825825825827
phi_pred = -0.2313916852932179 rad
signed alternative = A > 0
```

Report separately for H and G:

- free-phase fixed-omega amplitude;
- free-phase fixed-omega phase;
- phase offset from `phi_pred`;
- fixed-omega two-parameter delta chi-square;
- phase-locked signed amplitude;
- phase-locked positive-only delta chi-square;
- one-sided local analytic diagnostic;
- residual RMS;
- maximum absolute residual.

This post-selection evaluation does not feed back into background choice.

## Interpretation lock

- If a flexible spline/high-degree model wins decisively by held-out prediction and removes the frozen waveform, that is evidence that the degree-7 continuum was underfit for this statistic. It does not prove the removed structure is nonphysical, but it prevents treating the degree-7 10-sigma diagnostic as model-robust significance.
- If a moderate-complexity model wins and the positive coherent waveform survives, that strengthens the case that the baseline feature is not merely a degree-7 underfit artifact.
- If candidate scores are nearly tied, the background family is not identified strongly enough to let one arbitrary winner determine the physical conclusion; the tied/near-tied models must be reported together.
- Because H/G have already been observed, none of these outcomes count as a new prospective replication. The frozen Run2016F phase-locked test remains the prospective evidence layer.

## Canonical command

```bash
python scripts/run_hg_background_cv.py
```

Outputs:

```text
results/hg_background_cv/selection_freeze.json
results/hg_background_cv/cv_scores.csv
results/hg_background_cv/summary.json
```
