# DTU validation power-planning guide

**Status:** author/statistician decision document. It must be completed and
approved before validation outcomes are inspected.

Record the approved inputs and exact selected grid rows in
`DTU_VALIDATION_DECISION_RECORD_TEMPLATE.md`, sign it and checksum the frozen
copy before outcome access.

## What the grid calculates

The primary estimand is the E15.5 midbrain logit-usage curvature minus the mean
corresponding forebrain and hindbrain curvature. For matched regions from one
embryo, define the regional contrast as `M - (F + H)/2`. Under an equal-variance,
compound-symmetric planning model with marginal standard deviation `sigma` and
within-embryo cross-region correlation `rho`, the standard error of the complete
three-stage contrast is

`SE(Delta) = 1.5 * sigma * sqrt((1 - rho) / n)`,

where `n` is the number of independent biological units per stage. The released
grid uses a two-sided Bonferroni threshold for six primary genes. This is more
conservative than the planned BH analysis and makes the planning threshold
fixed before outcomes are known.

`tables/validation_power_sensitivity.csv` reports 80% and 90% marginal power
over standardized effects `abs(Delta)/sigma` from 0.75 to 2.00 and within-embryo
region correlations from 0 to 0.75. The calculation uses a normal approximation
and applies an operational floor of six independent units per stage.

`tables/validation_joint_panel_power_sensitivity.csv` evaluates the complete
suggested panel rule: at least four of six genes must pass the two-sided
six-test Bonferroni threshold in the expected direction, including at least
one of the two scan-led and one of the two cross-model genes, with no gene
significant in the opposite direction. It spans 6--32 independent units per
stage, the same standardized effects and region correlations, between-gene
test-statistic correlations from 0 to 0.75, and independent per-gene assay-
failure probabilities of 0%, 10% and 20%. The calculation uses deterministic
integration of an equicorrelated normal one-factor model. It assumes the same
positive effect and variance for all six genes; this is a sensitivity model,
not a claim that the candidates are exchangeable.

`tables/validation_collection_inflation.csv` gives the minimum number of
independent units to recruit per stage to retain each analyzable target with
80%, 90% or 95% assurance under 5%, 10% or 20% whole-unit loss. A whole unit is
an independent litter for which all three matched regions remain usable.

## Decisions that still must be frozen

1. Define the smallest regional-curvature effect worth detecting on the logit
   scale. Do not use the largest discovery estimate as the target effect.
2. Supply an outcome-independent variance/covariance estimate from a blinded
   assay-feasibility run, an external dataset or a justified conservative bound.
3. Define the independent unit. Prefer independent litters with one prespecified
   embryo per litter; additional embryos from one litter are subsamples unless
   litter clustering is modelled and the sample size is inflated.
4. Select the assay-failure, whole-unit-loss and retention-assurance scenarios
   to use for recruitment. Replacing samples after viewing gene results is
   prohibited.
5. Freeze the model family, zero handling, variance estimator, blocking terms
   and permitted covariates before unblinding.
6. Replace the joint grid's equal-effect, equicorrelation and independent-
   failure assumptions with blinded feasibility estimates if those estimates
   become available; otherwise prospectively select a conservative released
   scenario.

## Interpretation boundary

The six-unit floor is an operational minimum, not a universal power
justification. For example, at 90% marginal power and within-embryo correlation
0.50, the normal approximation requires 8 independent units per stage for a
standardized effect of 1.50 and 12 for an effect of 1.25. Joint-panel power can
be lower or higher than marginal power depending on between-gene correlation,
the four-of-six rule and assay failure. Recruitment must then be inflated for
whole-unit loss. A final number cannot be selected honestly until the meaningful
effect and independent validation variance are supplied and a released joint
scenario is chosen prospectively.
