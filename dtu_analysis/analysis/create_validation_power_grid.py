"""Create a prospective power-sensitivity grid for the DTU validation design.

The calculation is deliberately expressed in standardized-effect units because
the discovery archive cannot supply an independent validation variance. It is
a planning aid, not a post hoc power claim or a substitute for statistician
approval of the final covariance model.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from statistics import NormalDist


HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
TABLES = PAPER / "tables"

FAMILY_ALPHA = 0.05
PRIMARY_TESTS = 6
POWERS = (0.80, 0.90)
STANDARDIZED_EFFECTS = (0.75, 1.00, 1.25, 1.50, 2.00)
WITHIN_EMBRYO_REGION_CORRELATIONS = (0.00, 0.25, 0.50, 0.75)
OPERATIONAL_FLOOR_PER_STAGE = 6
PLANNING_UNITS_PER_STAGE = (6, 8, 10, 12, 16, 20, 24, 32)
BETWEEN_GENE_TEST_CORRELATIONS = (0.00, 0.25, 0.50, 0.75)
PER_GENE_ASSAY_FAILURE_PROBABILITIES = (0.00, 0.10, 0.20)
PANEL_MINIMUM_SUCCESSES = 4
WHOLE_UNIT_LOSS_PROBABILITIES = (0.05, 0.10, 0.20)
COLLECTION_ASSURANCES = (0.80, 0.90, 0.95)


def required_independent_units(
    standardized_effect: float,
    within_embryo_region_correlation: float,
    power: float,
) -> tuple[float, int]:
    """Return continuous and ceiling n for the balanced three-stage design.

    With equal marginal variance sigma^2 and compound-symmetric correlation
    rho among the three matched regions, the regional contrast
    M - (F + H)/2 has variance 1.5*sigma^2*(1-rho). The temporal curvature
    weights (1, -0.5, -0.5) add a second factor of 1.5, so
    SE(Delta) = 1.5*sigma*sqrt((1-rho)/n).
    """

    if standardized_effect <= 0:
        raise ValueError("standardized_effect must be positive")
    if not 0 <= within_embryo_region_correlation < 1:
        raise ValueError("correlation must be in [0, 1)")
    if not 0 < power < 1:
        raise ValueError("power must be in (0, 1)")

    normal = NormalDist()
    per_test_alpha = FAMILY_ALPHA / PRIMARY_TESTS
    critical = normal.inv_cdf(1 - per_test_alpha / 2)
    target = normal.inv_cdf(power)
    numerator = 1.5 * math.sqrt(1 - within_embryo_region_correlation)
    continuous = ((critical + target) * numerator / standardized_effect) ** 2
    return continuous, math.ceil(continuous)


def build_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for power in POWERS:
        for correlation in WITHIN_EMBRYO_REGION_CORRELATIONS:
            for effect in STANDARDIZED_EFFECTS:
                continuous, ceiling = required_independent_units(
                    effect, correlation, power
                )
                recommended = max(ceiling, OPERATIONAL_FLOOR_PER_STAGE)
                rows.append(
                    {
                        "family_alpha": FAMILY_ALPHA,
                        "primary_tests": PRIMARY_TESTS,
                        "planning_adjustment": "two-sided Bonferroni",
                        "target_power": power,
                        "standardized_regional_curvature_effect_delta_over_sigma": effect,
                        "within_embryo_cross_region_correlation": correlation,
                        "continuous_n_independent_units_per_stage": round(continuous, 6),
                        "ceiling_n_independent_units_per_stage": ceiling,
                        "operational_floor_per_stage": OPERATIONAL_FLOOR_PER_STAGE,
                        "planning_n_independent_units_per_stage": recommended,
                        "total_independent_units_three_stages": 3 * recommended,
                        "total_region_specimens_three_regions_three_stages": 9
                        * recommended,
                    }
                )
    return rows


def binomial_tail(total: int, success_probability: float, minimum: int) -> float:
    """Return P[X >= minimum] for X ~ Binomial(total, success_probability)."""

    return sum(
        math.comb(total, successes)
        * success_probability**successes
        * (1 - success_probability) ** (total - successes)
        for successes in range(minimum, total + 1)
    )


def marginal_directional_detection_probability(
    n_per_stage: int,
    standardized_effect: float,
    within_embryo_region_correlation: float,
) -> tuple[float, float, float]:
    """Return directional power, noncentrality and the two-sided threshold.

    A successful primary result must reject the two-sided six-test Bonferroni
    test and have the prespecified positive direction. The tiny probability of
    a significant result in the wrong direction is therefore not counted.
    """

    normal = NormalDist()
    critical = normal.inv_cdf(1 - (FAMILY_ALPHA / PRIMARY_TESTS) / 2)
    noncentrality = (
        standardized_effect
        * math.sqrt(n_per_stage)
        / (1.5 * math.sqrt(1 - within_embryo_region_correlation))
    )
    probability = 1 - normal.cdf(critical - noncentrality)
    return probability, noncentrality, critical


def full_panel_rule_probability(
    positive_probability: float,
    opposite_probability: float,
) -> float:
    """Probability of the complete panel rule for three two-gene strata.

    The strata are scan-led, cross-model and the two remaining evidence axes.
    Success requires at least four positive genes overall, at least one from
    each of the first two strata, and zero opposite-direction significant
    genes. Assay failure is part of the residual outcome.
    """

    residual_probability = 1 - positive_probability - opposite_probability
    residual_probability = min(1.0, max(0.0, residual_probability))

    def stratum_probability(positive_count: int) -> float:
        return (
            math.comb(2, positive_count)
            * positive_probability**positive_count
            * residual_probability ** (2 - positive_count)
        )

    total = 0.0
    for scan_positive in range(3):
        for cross_model_positive in range(3):
            for other_positive in range(3):
                if scan_positive < 1 or cross_model_positive < 1:
                    continue
                if (
                    scan_positive + cross_model_positive + other_positive
                    < PANEL_MINIMUM_SUCCESSES
                ):
                    continue
                total += (
                    stratum_probability(scan_positive)
                    * stratum_probability(cross_model_positive)
                    * stratum_probability(other_positive)
                )
    return total


def joint_panel_success_probability(
    n_per_stage: int,
    standardized_effect: float,
    within_embryo_region_correlation: float,
    between_gene_test_correlation: float,
    per_gene_assay_failure_probability: float,
) -> tuple[float, float, float, float, float]:
    """Return the complete six-gene panel success probability.

    Test statistics follow a one-factor equicorrelated normal planning model.
    Per-gene assay failures are independent of each other and of test
    statistics. Numerical integration is deterministic and uses a trapezoid
    rule over eight standard deviations in each direction.
    """

    marginal, noncentrality, critical = marginal_directional_detection_probability(
        n_per_stage,
        standardized_effect,
        within_embryo_region_correlation,
    )
    assay_success = 1 - per_gene_assay_failure_probability
    normal = NormalDist()
    opposite_marginal = normal.cdf(-critical - noncentrality)
    if between_gene_test_correlation == 0:
        joint = full_panel_rule_probability(
            assay_success * marginal,
            assay_success * opposite_marginal,
        )
        return joint, marginal, opposite_marginal, noncentrality, critical

    lower = -8.0
    upper = 8.0
    intervals = 4000
    step = (upper - lower) / intervals
    shared_scale = math.sqrt(between_gene_test_correlation)
    residual_scale = math.sqrt(1 - between_gene_test_correlation)

    def integrand(shared_factor: float) -> float:
        conditional_positive = 1 - normal.cdf(
            (
                critical
                - noncentrality
                - shared_scale * shared_factor
            )
            / residual_scale
        )
        conditional_opposite = normal.cdf(
            (
                -critical
                - noncentrality
                - shared_scale * shared_factor
            )
            / residual_scale
        )
        panel_rule = full_panel_rule_probability(
            assay_success * conditional_positive,
            assay_success * conditional_opposite,
        )
        density = math.exp(-(shared_factor**2) / 2) / math.sqrt(2 * math.pi)
        return panel_rule * density

    total = 0.5 * (integrand(lower) + integrand(upper))
    for index in range(1, intervals):
        total += integrand(lower + index * step)
    joint = total * step
    return joint, marginal, opposite_marginal, noncentrality, critical


def build_joint_panel_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for n_per_stage in PLANNING_UNITS_PER_STAGE:
        for region_correlation in WITHIN_EMBRYO_REGION_CORRELATIONS:
            for effect in STANDARDIZED_EFFECTS:
                for gene_correlation in BETWEEN_GENE_TEST_CORRELATIONS:
                    for failure_probability in PER_GENE_ASSAY_FAILURE_PROBABILITIES:
                        joint, marginal, opposite, noncentrality, critical = (
                            joint_panel_success_probability(
                                n_per_stage,
                                effect,
                                region_correlation,
                                gene_correlation,
                                failure_probability,
                            )
                        )
                        evaluable_and_detected = (
                            (1 - failure_probability) * marginal
                        )
                        rows.append(
                            {
                                "family_alpha": FAMILY_ALPHA,
                                "primary_tests": PRIMARY_TESTS,
                                "per_test_two_sided_alpha": round(
                                    FAMILY_ALPHA / PRIMARY_TESTS, 9
                                ),
                                "directional_success_threshold_z": round(
                                    critical, 6
                                ),
                                "panel_success_rule": "at least 4 of 6 genes significant in the prespecified direction; at least 1 of 2 scan-led and 1 of 2 cross-model genes; no opposite-direction significant gene",
                                "planned_independent_units_per_stage": n_per_stage,
                                "total_independent_units_three_stages": 3 * n_per_stage,
                                "standardized_regional_curvature_effect_delta_over_sigma": effect,
                                "within_embryo_cross_region_correlation": region_correlation,
                                "between_gene_test_statistic_correlation": gene_correlation,
                                "per_gene_assay_failure_probability": failure_probability,
                                "test_statistic_noncentrality": round(
                                    noncentrality, 6
                                ),
                                "marginal_expected_direction_detection_probability_before_failure": round(
                                    marginal, 6
                                ),
                                "marginal_opposite_direction_significance_probability_before_failure": round(
                                    opposite, 9
                                ),
                                "marginal_evaluable_and_detected_probability": round(
                                    evaluable_and_detected, 6
                                ),
                                "joint_panel_success_probability": round(joint, 6),
                                "joint_panel_power_at_least_0_80": joint >= 0.80,
                                "joint_panel_power_at_least_0_90": joint >= 0.90,
                                "planning_model_boundary": "equal positive effect and variance across genes; equicorrelated normal test statistics; independent per-gene assay failure; two genes in each of scan-led, cross-model and remaining strata",
                            }
                        )
    return rows


def minimum_recruitment(
    analyzable_target: int,
    whole_unit_loss_probability: float,
    assurance: float,
) -> tuple[int, float]:
    """Find the minimum recruited units giving the requested retention assurance."""

    recruited = analyzable_target
    while True:
        probability = binomial_tail(
            recruited,
            1 - whole_unit_loss_probability,
            analyzable_target,
        )
        if probability >= assurance:
            return recruited, probability
        recruited += 1


def build_collection_inflation_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for target in PLANNING_UNITS_PER_STAGE:
        for loss_probability in WHOLE_UNIT_LOSS_PROBABILITIES:
            for assurance in COLLECTION_ASSURANCES:
                recruited, achieved = minimum_recruitment(
                    target,
                    loss_probability,
                    assurance,
                )
                rows.append(
                    {
                        "target_analyzable_independent_units_per_stage": target,
                        "whole_unit_loss_probability": loss_probability,
                        "retention_assurance_target": assurance,
                        "minimum_units_to_recruit_per_stage": recruited,
                        "total_units_to_recruit_three_stages": 3 * recruited,
                        "achieved_probability_of_retaining_target": round(
                            achieved, 6
                        ),
                        "inflation_units_per_stage": recruited - target,
                        "unit_definition": "independent litter with all three matched regions usable",
                        "boundary": "binomial independent whole-unit loss; gene-level assay failure is handled separately in the joint-panel grid",
                    }
                )
    return rows


def write_csv(rows: list[dict[str, object]], filename: str) -> None:
    path = TABLES / filename
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_guide() -> None:
    text = """# DTU validation power-planning guide

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
"""
    (PAPER / "DTU_VALIDATION_POWER_GUIDE.md").write_text(text, encoding="utf-8")


def main() -> None:
    rows = build_rows()
    joint_rows = build_joint_panel_rows()
    inflation_rows = build_collection_inflation_rows()
    write_csv(rows, "validation_power_sensitivity.csv")
    write_csv(
        joint_rows,
        "validation_joint_panel_power_sensitivity.csv",
    )
    write_csv(
        inflation_rows,
        "validation_collection_inflation.csv",
    )
    write_guide()
    print(
        f"Wrote {len(rows)} marginal-power, {len(joint_rows)} joint-panel and "
        f"{len(inflation_rows)} collection-inflation scenarios; operational "
        f"floor {OPERATIONAL_FLOOR_PER_STAGE} independent units per stage."
    )


if __name__ == "__main__":
    main()
