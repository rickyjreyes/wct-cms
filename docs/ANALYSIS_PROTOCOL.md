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

The frequency is not fit to CMS. Only amplitude, phase, and intercept are fit. The repository reports the fixed-frequency improvement and a permutation p-value.

**Exploratory diagnostic:** omega scan.

The maximum scan statistic is corrected using the distribution of the maximum statistic obtained in residual permutations. It must not be reported as though it were the pre-registered replication test.

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

## Permutation p-value resolution

The residual-permutation null is an exploratory diagnostic. With `N` permutations and the repository's add-one correction, the smallest reportable Monte Carlo p-value is

\[
p_{\min}=\frac{1}{N+1}.
\]

Thus 20 permutations can report no value below `1/21 = 0.047619...`; that value means only that zero of 20 shuffled residual sets exceeded the observed statistic. It is not evidence that the true p-value has been measured to approximately 0.048.

Residual permutation also assumes the residuals are exchangeable after background fitting. A serious physics result requires stronger null controls, including parametric/bootstrap studies that repeat the background fit and tests on simulation/control regions.

## Controls required before a physics claim

1. Repeat with multiple smooth-background families/degrees.
2. Vary binning without choosing the variation that maximizes significance.
3. Repeat across data-taking eras.
4. Run CMS simulation/background samples through the identical pipeline.
5. Test known detector/reconstruction resonances and predeclared masked windows.
6. Require an interior frequency with enough cycles to be spectrally resolvable.
7. Replace or supplement residual permutation with a null that propagates background-fit uncertainty.
8. Test at least one independent CMS channel, preferably photons after dimuons.
9. Apply the same dimensionless mapping used to compare GWTC, LHC, JUNO and photodiode data.
10. Record all failed as well as successful runs.

A CMS match by itself is not evidence that WCT caused the structure. The intended value is independent cross-validation under a frozen prediction and systematic controls.
