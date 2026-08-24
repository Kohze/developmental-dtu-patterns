"""Adversarial, output-level audit of the DTU manuscript package.

The audit assumes that headline counts, trajectory calls, validation power and
provenance claims may be wrong. It independently reconstructs their released
invariants without importing analysis functions. Any failed local check exits
non-zero after writing a machine-readable report. External evidence boundaries
are recorded separately and do not masquerade as passes.
"""

from __future__ import annotations

import csv
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
DATA = PAPER / "data"
TABLES = PAPER / "tables"
STAGES = ("10.5", "11.5", "12.5", "13.5", "14.5", "15.5", "16.5", "0")
REGIONS = ("Forebrain", "Hindbrain", "Midbrain")
LEADING_SIX = ("Scg3", "Gpm6a", "Ntrk2", "Tecr", "Armc8", "Bin1")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def number(value: str) -> float:
    if value.strip().upper() in {"NA", "NAN", ""}:
        return math.nan
    return float(value)


def close(left: float, right: float, tolerance: float = 1e-9) -> bool:
    if math.isnan(left) and math.isnan(right):
        return True
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


checks: list[dict[str, str]] = []


def check(
    check_id: str,
    adversarial_question: str,
    passed: bool,
    evidence: str,
    consequence_if_failed: str,
) -> None:
    checks.append(
        {
            "check_id": check_id,
            "adversarial_question": adversarial_question,
            "status": "PASS" if passed else "FAIL",
            "evidence": evidence,
            "consequence_if_failed": consequence_if_failed,
        }
    )


def external_boundary(
    check_id: str,
    adversarial_question: str,
    evidence: str,
    consequence: str,
) -> None:
    checks.append(
        {
            "check_id": check_id,
            "adversarial_question": adversarial_question,
            "status": "EXTERNAL_BOUNDARY",
            "evidence": evidence,
            "consequence_if_failed": consequence,
        }
    )


def binomial_tail(total: int, success_probability: float, minimum: int) -> float:
    return sum(
        math.comb(total, successes)
        * success_probability**successes
        * (1 - success_probability) ** (total - successes)
        for successes in range(minimum, total + 1)
    )


def adjusted_p_values(values: list[float], method: str) -> list[float]:
    """Independent BH/BY implementation matching R's p.adjust definition."""

    total = len(values)
    multiplier = 1.0
    if method == "BY":
        multiplier = sum(1 / index for index in range(1, total + 1))
    elif method != "BH":
        raise ValueError(f"Unsupported adjustment method: {method}")
    order = sorted(range(total), key=lambda index: values[index])
    sorted_raw = [
        values[index] * total * multiplier / rank
        for rank, index in enumerate(order, start=1)
    ]
    running = 1.0
    sorted_adjusted = [1.0] * total
    for index in range(total - 1, -1, -1):
        running = min(running, sorted_raw[index])
        sorted_adjusted[index] = max(0.0, min(1.0, running))
    answer = [1.0] * total
    for index, original_index in enumerate(order):
        answer[original_index] = sorted_adjusted[index]
    return answer


def verify_episode_rows(
    rows: list[dict[str, str]],
    stage_lookup: dict[tuple[str, str, str], dict[str, float | str]],
    q_name: str,
    gene_q_name: str,
) -> tuple[int, int, int]:
    failures = 0
    reciprocal_groups: dict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    for row in rows:
        start = STAGES.index(row["start_stage"])
        end = STAGES.index(row["end_stage"])
        if start <= 0 or end >= len(STAGES) - 1 or end < start:
            failures += 1
            continue
        stages = [
            stage_lookup.get((row["isoform_id"], row["region"], stage))
            for stage in STAGES
        ]
        if any(value is None for value in stages):
            failures += 1
            continue
        internal = stages[start : end + 1]
        direction = 1 if row["direction"] == "higher" else -1
        for value in internal:
            first = float(value["target_vs_other_1"])
            second = float(value["target_vs_other_2"])
            qualifies = (
                float(value[q_name]) < 0.05
                and float(value[gene_q_name]) < 0.05
                and abs(first) >= 0.10
                and abs(second) >= 0.10
                and abs(float(value["other_region_difference"])) < 0.10
                and (1 if first > 0 else -1) == direction
                and (1 if second > 0 else -1) == direction
            )
            if not qualifies:
                failures += 1
        flanks = (stages[start - 1], stages[end + 1])
        flank_values = [
            abs(number(str(value[field])))
            for value in flanks
            for field in ("target_vs_other_1", "target_vs_other_2")
        ]
        flank_max = max(value for value in flank_values if math.isfinite(value))
        if flank_max >= 0.10 or not close(
            flank_max, number(row["flanking_max_abs_difference"])
        ):
            failures += 1
        max_effect = max(abs(float(value["target_difference"])) for value in internal)
        mean_effect = sum(abs(float(value["target_difference"])) for value in internal) / len(internal)
        worst_q = max(float(value[q_name]) for value in internal)
        if not (
            close(max_effect, number(row["max_abs_usage_difference"]))
            and close(mean_effect, number(row["mean_abs_usage_difference"]))
            and close(worst_q, number(row["worst_pair_q"]))
            and int(row["n_stages"]) == end - start + 1
            and number(row["replicate_separation"]) > 0
            and row["replicate_consistent"].upper() == "TRUE"
        ):
            failures += 1
        group = (row["gene_id"], row["region"], row["start_stage"], row["end_stage"])
        reciprocal_groups[group].add(row["direction"])

    reciprocal_mismatches = 0
    for row in rows:
        group = (row["gene_id"], row["region"], row["start_stage"], row["end_stage"])
        expected = len(reciprocal_groups[group]) > 1
        observed = row["reciprocal_exchange"].upper() == "TRUE"
        reciprocal_mismatches += expected != observed
    return failures, reciprocal_mismatches, len(reciprocal_groups)


def main() -> None:
    diagnostics = {
        row["metric"]: int(float(row["value"]))
        for row in read_rows(TABLES / "transient_regional_scan_diagnostics.csv")
    }
    expected_tests = (
        diagnostics["isoforms_after_expression_and_multi_isoform_filtering"]
        * 3
        * diagnostics["region_stage_cells"]
        // 3
    )
    # region_stage_cells/3 is the eight-stage count; three is the region-pair count.
    check(
        "AA01",
        "Is the advertised complete test family arithmetically complete?",
        expected_tests == 300408 == diagnostics["pairwise_tests"],
        "12,517 isoforms x 3 regional pairs x 8 stages = 300,408 tests.",
        "The global adjustment family or manuscript denominator would be wrong.",
    )

    pair_lookup: dict[tuple[str, str, str, str], tuple[float, float, float, str, float, float]] = {}
    isoform_pair_counts: Counter[str] = Counter()
    pair_genes: set[str] = set()
    component_p_values: list[float] = []
    observed_bh_values: list[float] = []
    observed_by_values: list[float] = []
    gene_component_p_values: dict[str, list[float]] = defaultdict(list)
    observed_gene_values: dict[str, tuple[float, float, float, float]] = {}
    adjustment_failures = 0
    with (DATA / "transient_regional_pair_tests_all.csv").open(
        encoding="utf-8-sig", newline=""
    ) as handle:
        for row in csv.DictReader(handle):
            key = (row["isoform_id"], row["region_1"], row["region_2"], row["stage"])
            if key in pair_lookup:
                adjustment_failures += 1
            p_value = number(row["p_value"])
            q_value = number(row["isoform_q_global"])
            q_by = number(row["isoform_q_by"])
            gene_q = number(row["gene_q"])
            gene_q_by = number(row["gene_q_by"])
            gene_simes = number(row["gene_simes_p"])
            gene_bonferroni = number(row["gene_bonferroni_p"])
            values = (p_value, q_value, q_by, gene_q, gene_q_by)
            if not all(math.isfinite(value) and 0 <= value <= 1 for value in values):
                adjustment_failures += 1
            if q_value + 1e-15 < p_value or q_by + 1e-15 < q_value or gene_q_by + 1e-15 < gene_q:
                adjustment_failures += 1
            pair_lookup[key] = (
                number(row["usage_difference"]),
                q_value,
                q_by,
                row["gene_id"],
                gene_q,
                gene_q_by,
            )
            isoform_pair_counts[row["isoform_id"]] += 1
            pair_genes.add(row["gene_id"])
            component_p_values.append(p_value)
            observed_bh_values.append(q_value)
            observed_by_values.append(q_by)
            gene_component_p_values[row["gene_id"]].append(p_value)
            current_gene_values = (
                gene_simes,
                gene_q,
                gene_bonferroni,
                gene_q_by,
            )
            previous_gene_values = observed_gene_values.setdefault(
                row["gene_id"], current_gene_values
            )
            if not all(
                close(left, right)
                for left, right in zip(previous_gene_values, current_gene_values)
            ):
                adjustment_failures += 1

    recalculated_bh = adjusted_p_values(component_p_values, "BH")
    recalculated_by = adjusted_p_values(component_p_values, "BY")
    component_bh_mismatches = sum(
        not close(observed, recalculated, 2e-12)
        for observed, recalculated in zip(observed_bh_values, recalculated_bh)
    )
    component_by_mismatches = sum(
        not close(observed, recalculated, 2e-12)
        for observed, recalculated in zip(observed_by_values, recalculated_by)
    )
    adjustment_failures += component_bh_mismatches + component_by_mismatches

    gene_order = sorted(gene_component_p_values)
    recalculated_simes: list[float] = []
    recalculated_bonferroni: list[float] = []
    for gene in gene_order:
        values = sorted(gene_component_p_values[gene])
        count = len(values)
        recalculated_simes.append(
            min(1.0, min(value * count / rank for rank, value in enumerate(values, 1)))
        )
        recalculated_bonferroni.append(min(1.0, min(values) * count))
    recalculated_gene_bh = adjusted_p_values(recalculated_simes, "BH")
    recalculated_gene_by = adjusted_p_values(recalculated_bonferroni, "BY")
    gene_adjustment_mismatches = 0
    for index, gene in enumerate(gene_order):
        observed = observed_gene_values[gene]
        recalculated = (
            recalculated_simes[index],
            recalculated_gene_bh[index],
            recalculated_bonferroni[index],
            recalculated_gene_by[index],
        )
        gene_adjustment_mismatches += any(
            not close(left, right, 2e-12)
            for left, right in zip(observed, recalculated)
        )
    adjustment_failures += gene_adjustment_mismatches
    pair_complete = (
        len(pair_lookup) == 300408
        and len(isoform_pair_counts) == 12517
        and set(isoform_pair_counts.values()) == {24}
        and len(pair_genes) == 4577
    )
    check(
        "AA02",
        "Are all component tests present once, finite after fail-closed replacement, and monotonically adjusted?",
        pair_complete and adjustment_failures == 0,
        f"{len(pair_lookup):,} unique rows; {len(isoform_pair_counts):,} isoforms x 24; "
        f"{len(pair_genes):,} genes; exact BH mismatches={component_bh_mismatches}; "
        f"exact BY mismatches={component_by_mismatches}; exact gene-adjustment "
        f"mismatches={gene_adjustment_mismatches}; other violations="
        f"{adjustment_failures - component_bh_mismatches - component_by_mismatches - gene_adjustment_mismatches}.",
        "The global/BY sensitivity results could not be independently trusted.",
    )

    primary_episodes = read_rows(DATA / "transient_regional_isoform_episodes.csv")
    conservative_episodes = read_rows(
        DATA / "transient_regional_isoform_episodes_conservative.csv"
    )
    needed_groups = {
        (row["isoform_id"], row["region"])
        for row in primary_episodes + conservative_episodes
    }
    needed_isoforms = {isoform for isoform, _ in needed_groups}
    fraction_lookup: dict[str, dict[str, float]] = {}
    fraction_row_count = 0
    with (DATA / "transient_regional_filtered_isoform_fractions.csv").open(
        encoding="utf-8-sig", newline=""
    ) as handle:
        reader = csv.DictReader(handle)
        fraction_sample_ids = [
            field for field in reader.fieldnames or [] if field != "isoform_id"
        ]
        for row in reader:
            fraction_row_count += 1
            if row["isoform_id"] in needed_isoforms:
                fraction_lookup[row["isoform_id"]] = {
                    sample: number(row[sample]) for sample in fraction_sample_ids
                }
    sample_metadata: dict[str, tuple[str, str]] = {}
    for sample in fraction_sample_ids:
        region, stage_replicate = sample.split("__", 1)
        stage = re.sub(r"_[12]$", "", stage_replicate)
        sample_metadata[sample] = (region, stage)
    fraction_matrix_complete = (
        fraction_row_count == 12517
        and len(fraction_sample_ids) == 48
        and set(fraction_lookup) == needed_isoforms
    )

    def recompute_replicate_separation(row: dict[str, str]) -> float:
        values = fraction_lookup[row["isoform_id"]]
        start = STAGES.index(row["start_stage"])
        end = STAGES.index(row["end_stage"])
        stage_separations: list[float] = []
        for stage in STAGES[start : end + 1]:
            target_values = [
                value
                for sample, value in values.items()
                if sample_metadata[sample] == (row["region"], stage)
            ]
            other_values = [
                value
                for sample, value in values.items()
                if sample_metadata[sample][0] != row["region"]
                and sample_metadata[sample][1] == stage
            ]
            if len(target_values) != 2 or len(other_values) != 4:
                return math.nan
            if not all(math.isfinite(value) for value in target_values + other_values):
                return math.nan
            if row["direction"] == "higher":
                stage_separations.append(min(target_values) - max(other_values))
            else:
                stage_separations.append(min(other_values) - max(target_values))
        return min(stage_separations)

    primary_separation_mismatches = sum(
        not close(
            recompute_replicate_separation(row), number(row["replicate_separation"])
        )
        for row in primary_episodes
    )
    conservative_separation_mismatches = sum(
        not close(
            recompute_replicate_separation(row), number(row["replicate_separation"])
        )
        for row in conservative_episodes
    )
    stage_lookup: dict[tuple[str, str, str], dict[str, float | str]] = {}
    isoform_stage_counts: Counter[str] = Counter()
    stage_failures = 0
    with (DATA / "transient_regional_stage_evaluations_all.csv").open(
        encoding="utf-8-sig", newline=""
    ) as handle:
        for row in csv.DictReader(handle):
            isoform = row["isoform_id"]
            target = row["region"]
            stage = row["stage"]
            others = [region for region in REGIONS if region != target]
            oriented: list[float] = []
            q_values: list[float] = []
            q_by_values: list[float] = []
            for other in others:
                pair = sorted((target, other))
                pair_row = pair_lookup.get((isoform, pair[0], pair[1], stage))
                if pair_row is None:
                    stage_failures += 1
                    continue
                orientation = 1 if target == pair[0] else -1
                oriented.append(orientation * pair_row[0])
                q_values.append(pair_row[1])
                q_by_values.append(pair_row[2])
                if pair_row[3] != row["gene_id"]:
                    stage_failures += 1
            if len(oriented) != 2:
                stage_failures += 1
                continue
            expected = (
                oriented[0],
                oriented[1],
                max(q_values),
                max(q_by_values),
            )
            observed = (
                number(row["target_vs_other_1"]),
                number(row["target_vs_other_2"]),
                number(row["max_pair_q"]),
                number(row["max_pair_q_by"]),
            )
            if not all(close(left, right) for left, right in zip(expected, observed)):
                stage_failures += 1
            if not close(
                number(row["target_mean_if"]) - number(row["other_mean_if"]),
                number(row["target_difference"]),
            ):
                stage_failures += 1
            isoform_stage_counts[isoform] += 1
            if (isoform, target) in needed_groups:
                stage_lookup[(isoform, target, stage)] = row
    stage_complete = (
        sum(isoform_stage_counts.values()) == 300408
        and len(isoform_stage_counts) == 12517
        and set(isoform_stage_counts.values()) == {24}
    )
    check(
        "AA03",
        "Can every released target-region evaluation be reconstructed from the complete pair tests?",
        stage_complete and stage_failures == 0,
        f"300,408 stage evaluations; pair-orientation/arithmetic mismatches={stage_failures}.",
        "Episode direction, magnitude or q-value lineage could be internally inconsistent.",
    )
    del pair_lookup

    primary_failures, primary_reciprocal_failures, _ = verify_episode_rows(
        primary_episodes, stage_lookup, "max_pair_q", "gene_q"
    )
    check(
        "AA04",
        "Does every primary episode satisfy the published q, effect, geometry, flank and replicate rules?",
        len(primary_episodes) == 1348
        and fraction_matrix_complete
        and primary_failures == 0
        and primary_reciprocal_failures == 0
        and primary_separation_mismatches == 0,
        f"1,348 rows checked; fraction matrix complete={fraction_matrix_complete}; "
        f"criterion failures={primary_failures}; reciprocal-flag failures="
        f"{primary_reciprocal_failures}; replicate-separation mismatches="
        f"{primary_separation_mismatches}.",
        "The central episode count would contain calls not generated by the stated algorithm.",
    )

    conservative_failures, conservative_reciprocal_failures, _ = verify_episode_rows(
        conservative_episodes, stage_lookup, "max_pair_q_by", "gene_q_by"
    )
    episode_key_columns = (
        "isoform_id", "gene_id", "region", "start_stage", "end_stage", "direction"
    )
    primary_keys = {tuple(row[column] for column in episode_key_columns) for row in primary_episodes}
    conservative_keys = {
        tuple(row[column] for column in episode_key_columns)
        for row in conservative_episodes
    }
    check(
        "AA05",
        "Do all conservative BY/Bonferroni calls satisfy the conservative rules and remain a subset of the primary calls?",
        len(conservative_episodes) == 852
        and conservative_failures == 0
        and conservative_reciprocal_failures == 0
        and conservative_separation_mismatches == 0
        and conservative_keys <= primary_keys,
        f"852 rows checked; criterion failures={conservative_failures}; "
        f"reciprocal-flag failures={conservative_reciprocal_failures}; "
        f"replicate-separation mismatches={conservative_separation_mismatches}; "
        f"non-primary keys={len(conservative_keys - primary_keys)}.",
        "The arbitrary-dependence sensitivity would not be a reproducible conservative subset.",
    )

    primary_genes = {row["gene_id"] for row in primary_episodes}
    primary_e155 = sum(
        row["start_stage"] == "15.5" and row["end_stage"] == "15.5"
        for row in primary_episodes
    )
    conservative_genes = {row["gene_id"] for row in conservative_episodes}
    conservative_e155 = sum(
        row["start_stage"] == "15.5" and row["end_stage"] == "15.5"
        for row in conservative_episodes
    )
    all_midbrain = all(row["region"] == "Midbrain" for row in primary_episodes)
    check(
        "AA06",
        "Are headline episode totals derivable from row-level calls rather than summaries alone?",
        len(primary_genes) == 735
        and primary_e155 == 1203
        and len(primary_episodes) - primary_e155 == 145
        and len(conservative_genes) == 474
        and conservative_e155 == 817
        and len(conservative_episodes) - conservative_e155 == 35
        and all_midbrain,
        "Primary: 1,348/735, E15.5=1,203, other=145, all midbrain; "
        "conservative: 852/474, E15.5=817, other=35.",
        "The abstract's dominant quantitative result would be a summary-table artefact.",
    )

    dependence = read_rows(TABLES / "transient_regional_dependence_sensitivity.csv")
    leaders = [tuple(row["leading_six_genes"].split(";")) for row in dependence]
    top_candidates = read_rows(TABLES / "transient_regional_top_candidates.csv")
    check(
        "AA07",
        "Does the leading panel survive the arbitrary-dependence specification without manual substitution?",
        leaders == [LEADING_SIX, LEADING_SIX]
        and tuple(row["gene_id"] for row in top_candidates) == LEADING_SIX,
        "Primary and BY/Bonferroni leaders are Scg3, Gpm6a, Ntrk2, Tecr, Armc8 and Bin1.",
        "Candidate selection would be sensitive to the dependence correction or manual reranking.",
    )

    replicate_rows = read_rows(TABLES / "candidate_replicate_choice_audit.csv")
    replicate_ok = len(replicate_rows) == 6 and all(
        int(row["replicate_choice_combinations"]) == 512
        and number(row["joint_logit_expected_direction_fraction"]) == 1
        and number(row["joint_raw_expected_direction_fraction"]) == 1
        for row in replicate_rows
    )
    check(
        "AA08",
        "Can one favourable replicate suffix explain any frozen candidate direction?",
        replicate_ok,
        "All six accession pairs retain both expected directions in all 512 logit and raw-fraction selections.",
        "The panel could be an artefact of choosing one archived replicate suffix.",
    )

    assay = read_rows(TABLES / "candidate_assay_target_summary.csv")
    junctions = read_rows(TABLES / "candidate_junction_target_audit.csv")
    preferred = [row for row in junctions if row["preferred_coordinate_candidate"] == "TRUE"]
    armc8 = next(
        row for row in assay if row["gene"] == "Armc8" and row["expected_direction"] == "lower"
    )
    assay_ok = (
        len(assay) == 12
        and len(preferred) == 34
        and all(
            row["junction_kmer_length_nt"] == "40"
            and row["junction_kmer_occurrences_in_paired_accession"] == "0"
            and row["junction_kmer_occurrences_in_archived_same_gene_models"] == "1"
            for row in preferred
        )
        and armc8["preferred_junction_coordinate_candidates"] == "0"
        and armc8["longest_gene_unique_segment_nt"] == "2242"
    )
    check(
        "AA09",
        "Do the proposed archived targets actually discriminate their accession pairs within each archived gene?",
        assay_ok,
        "34 preferred 40-nt junction anchors pass within-gene checks; lower Armc8 has zero and uses a 2,242-nt terminal segment.",
        "The validation plan could specify an impossible or non-discriminating assay.",
    )

    sequences = read_rows(TABLES / "candidate_sequence_reconstruction_audit.csv")
    sequence_ok = (
        len(sequences) == 12
        and sum(row["archived_sequence_survives"] == "TRUE" for row in sequences) == 10
        and all(
            row["exact_match_to_archived_sequence"] == "TRUE"
            for row in sequences
            if row["archived_sequence_survives"] == "TRUE"
        )
        and sum(int(row["reconstructed_ambiguous_bases"]) for row in sequences) == 0
    )
    check(
        "AA10",
        "Does exon-plus-mm10 reconstruction contradict surviving archived candidate sequences?",
        sequence_ok,
        "10/10 surviving sequences match exactly; 12/12 reconstructions contain zero ambiguous bases; two Scg3 sequences lack comparators.",
        "Coordinate-to-sequence lineage would be unreliable.",
    )

    joint_rows = read_rows(TABLES / "validation_joint_panel_power_sensitivity.csv")
    joint_keys = {
        (
            row["planned_independent_units_per_stage"],
            row["standardized_regional_curvature_effect_delta_over_sigma"],
            row["within_embryo_cross_region_correlation"],
            row["between_gene_test_statistic_correlation"],
            row["per_gene_assay_failure_probability"],
        )
        for row in joint_rows
    }
    joint_ok = len(joint_rows) == 1920 and len(joint_keys) == 1920 and all(
        0 <= number(row["joint_panel_success_probability"]) <= 1
        for row in joint_rows
    )
    # A fully independent gene-statistic row has a closed-form three-stratum rule.
    closed_form_failures = 0
    for row in joint_rows:
        if number(row["between_gene_test_statistic_correlation"]) != 0:
            continue
        success = (1 - number(row["per_gene_assay_failure_probability"])) * number(
            row["marginal_expected_direction_detection_probability_before_failure"]
        )
        opposite = (1 - number(row["per_gene_assay_failure_probability"])) * number(
            row["marginal_opposite_direction_significance_probability_before_failure"]
        )
        residual = 1 - success - opposite

        def group_probability(count: int) -> float:
            return math.comb(2, count) * success**count * residual ** (2 - count)

        exact = sum(
            group_probability(scan)
            * group_probability(cross)
            * group_probability(other)
            for scan in range(3)
            for cross in range(3)
            for other in range(3)
            if scan >= 1 and cross >= 1 and scan + cross + other >= 4
        )
        closed_form_failures += abs(exact - number(row["joint_panel_success_probability"])) > 1.5e-6
    check(
        "AA11",
        "Does the joint-panel grid implement the full proposed rule rather than multiplying marginal power?",
        joint_ok and closed_form_failures == 0,
        f"1,920 unique bounded scenarios; independent-statistic closed-form failures={closed_form_failures}.",
        "The proposed recruitment could be underpowered for the actual compound success rule.",
    )

    inflation_rows = read_rows(TABLES / "validation_collection_inflation.csv")
    inflation_failures = 0
    for row in inflation_rows:
        target = int(row["target_analyzable_independent_units_per_stage"])
        loss = number(row["whole_unit_loss_probability"])
        assurance = number(row["retention_assurance_target"])
        recruited = int(row["minimum_units_to_recruit_per_stage"])
        achieved = binomial_tail(recruited, 1 - loss, target)
        if achieved + 1e-12 < assurance:
            inflation_failures += 1
        if recruited > target and binomial_tail(recruited - 1, 1 - loss, target) >= assurance:
            inflation_failures += 1
    check(
        "AA12",
        "Are recruitment inflations exact minima rather than rounded heuristics?",
        len(inflation_rows) == 72 and inflation_failures == 0,
        f"72 scenarios; non-minimal or under-assured rows={inflation_failures}.",
        "The collection could miss its analyzable-unit target more often than stated.",
    )

    claim_rows = read_rows(TABLES / "claim_evidence_ledger.csv")
    missing_claim_sources: list[str] = []
    for row in claim_rows:
        for source in row["primary_source_files"].split("; "):
            if not (PAPER / source).is_file():
                missing_claim_sources.append(f"{row['claim_id']}:{source}")
    claim_ids = [row["claim_id"] for row in claim_rows]
    check(
        "AA13",
        "Can every central claim be traced to an existing released source?",
        claim_ids == [f"DTU-C{index:02d}" for index in range(1, 13)]
        and not missing_claim_sources,
        f"12 ordered claims; missing source references={len(missing_claim_sources)}.",
        "The claim ledger would provide apparent rather than executable provenance.",
    )

    display_rows = read_rows(TABLES / "display_source_manifest.csv")
    missing_display_sources: list[str] = []
    external_display_sources = 0
    for row in display_rows:
        for field in ("rendered_asset", "source_file"):
            for source in row[field].split("; "):
                if not (PAPER / source).is_file():
                    missing_display_sources.append(f"{row['display']}:{source}")
        for source in row["regeneration_or_audit_source"].split("; "):
            if source.lower().startswith("external: "):
                external_display_sources += 1
            elif not (PAPER / source).is_file():
                missing_display_sources.append(f"{row['display']}:{source}")
    check(
        "AA14",
        "Can every displayed panel be traced to an extant source and regeneration script?",
        not missing_display_sources,
        f"{len(display_rows)} display mappings; explicit external-audit entries={external_display_sources}; "
        f"missing assets/sources={len(missing_display_sources)}.",
        "Figures or tables could not be independently regenerated from their claimed lineage.",
    )

    manuscript = (PAPER / "manuscript.tex").read_text(encoding="utf-8")
    required_language = (
        "neither specification establishes an episode-level false discovery rate",
        "not a demonstrated developmental programme",
        "not independent replication",
        "cannot distinguish a coordinated multi-cellular transition",
    )
    required_numbers = ("300,408", "1,348", "735", "1,203", "852", "474")
    check(
        "AA15",
        "Does the manuscript state the component-test/episode-FDR boundary and avoid converting discovery into mechanism?",
        all(text in manuscript for text in required_language)
        and all(text in manuscript for text in required_numbers),
        "Required no-episode-FDR, no-mechanism and non-independence language is present with all headline counts.",
        "Readers could interpret adjusted component tests as controlled episode discoveries or biological mechanism.",
    )

    release_config = json.loads((PAPER / "release_config.json").read_text(encoding="utf-8"))
    release_builder = (PAPER / "prepare_release.ps1").read_text(encoding="utf-8")
    paper_builder = (PAPER / "build.ps1").read_text(encoding="utf-8")
    required_release_files = {
        "ANTAGONISTIC_AUDIT.md",
        "DTU_VALIDATION_DECISION_RECORD_TEMPLATE.md",
        "DTU_VALIDATION_PREREGISTRATION_DRAFT.md",
        "DTU_VALIDATION_POWER_GUIDE.md",
        "release_config.json",
    }
    required_globs = {"analysis/*.py", "analysis/*.R", "data/*.csv", "tables/*.csv"}
    check(
        "AA16",
        "Does the deterministic release allow-list include the adversarial evidence needed to challenge the study?",
        required_release_files <= set(release_config["include_files"])
        and required_globs <= set(release_config["include_globs"])
        and "$releaseConfig.include_files" in release_builder
        and "$releaseConfig.include_globs" in release_builder
        and "SOURCE_DATE_EPOCH = '946684800'" in paper_builder
        and "FORCE_SOURCE_DATE = '1'" in paper_builder,
        "Audit report, decision record, protocol, power guide, all scripts and all data/table CSVs are allow-listed; the release builder consumes both configured fields and the PDF builder fixes source-date metadata.",
        "A locally passing audit could disappear from the public package.",
    )

    external_boundary(
        "AA17",
        "Has a genuinely independent cohort established biological replication?",
        "No: all current models reuse the same 48 archived sample columns.",
        "The E15.5 discontinuity remains a discovery anomaly, not a replicated developmental programme.",
    )
    external_boundary(
        "AA18",
        "Have full-length structure and cellular origin been established?",
        "No: exact primers/probes, product sequencing, long reads and matched neural/glial/vascular measurements require new material.",
        "Quantification artefact and changing cell composition remain viable explanations.",
    )
    external_boundary(
        "AA19",
        "Can the archive be reproduced from exact raw files and frozen quantifier/annotation versions?",
        "No: the original raw-file manifest, RefSeq release, Salmon/index version and complete hand-off are unresolved.",
        "Exact reproducibility begins at the archived regional objects, not FASTQ.",
    )
    external_boundary(
        "AA20",
        "Is the package legally and administratively depositable?",
        "No: licence, source-rights confirmation, final authorship metadata and persistent repository identifier are absent.",
        "The audit build is not an authorized public release.",
    )

    output_csv = TABLES / "antagonistic_audit_checks.csv"
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(checks[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(checks)

    failures = [row for row in checks if row["status"] == "FAIL"]
    passes = [row for row in checks if row["status"] == "PASS"]
    boundaries = [row for row in checks if row["status"] == "EXTERNAL_BOUNDARY"]
    verdict = "FAIL" if failures else "PASS WITH EXTERNAL EVIDENCE BOUNDARIES"
    lines = [
        "# DTU antagonistic audit",
        "",
        "Audit date: 31 July 2026  ",
        "Audit implementation: `analysis/run_antagonistic_audit.py`",
        "",
        "## Threat model",
        "",
        "This audit assumes that the headline family size, adjusted values,",
        "episode geometry, candidate selection, assay targets, prospective power",
        "and provenance may each be wrong. It reconstructs released invariants",
        "from complete derivative tables rather than trusting manuscript prose or",
        "analysis summaries. External biological and legal questions are marked",
        "as boundaries, not converted into computational passes.",
        "",
        "## Verdict",
        "",
        f"**{verdict}.** Local passes: {len(passes)}; local failures: {len(failures)}; "
        f"external boundaries: {len(boundaries)}.",
        "",
        "## Check matrix",
        "",
        "| ID | Status | Adversarial question | Evidence | Consequence if failed/unresolved |",
        "|---|---|---|---|---|",
    ]
    for row in checks:
        values = [
            row["check_id"],
            row["status"],
            row["adversarial_question"],
            row["evidence"],
            row["consequence_if_failed"],
        ]
        values = [value.replace("|", "\\|").replace("\n", " ") for value in values]
        lines.append("| " + " | ".join(values) + " |")
    lines.extend(
        [
            "",
            "## Substantiated issues corrected before this verdict",
            "",
            "- The earlier release exposed only significant component tests and",
            "  primary episodes. It now includes all 300,408 pair tests, all",
            "  300,408 target-region stage evaluations and all 852 conservative",
            "  episode calls, enabling output-level trajectory reconstruction.",
            "- Component-test adjustment was potentially readable as episode-level",
            "  FDR control. The manuscript and claim ledger now state explicitly",
            "  that no episode-level false discovery rate is established.",
            "- Replicate separation was previously present only as a stored episode",
            "  field. The complete filtered 12,517-by-48 fraction matrix is now",
            "  released, and every primary/conservative separation is recomputed.",
            "- Every component BH/BY value and every gene Simes/Bonferroni plus",
            "  BH/BY value is now independently recalculated, not merely checked",
            "  for monotonicity.",
            "- The PDF workflow previously trusted native exit codes even if a",
            "  stale artifact remained. It now requires a freshly written PDF/log",
            "  pair, an explicit successful-output record in the LaTeX log and a",
            "  fixed source-date epoch so repeated compiles are byte-identical.",
            "- The release builder previously duplicated and drifted from its",
            "  configured allow-list, omitting this report and the decision record.",
            "  It now consumes that allow-list directly, includes its own config and",
            "  uses explicit ordinal ordering for manifests and ZIP entries.",
            "",
            "## Interpretation",
            "",
            "A clean local verdict means the released calculations, tables, claims",
            "and boundaries agree with each other under the stated archived-object",
            "design. It does not authenticate the E15.5 signal biologically. That",
            "requires the independent, balanced and preregistered experiment in the",
            "validation protocol, together with finalized rights and authorship.",
            "",
        ]
    )
    (PAPER / "ANTAGONISTIC_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")

    print(
        f"Antagonistic audit: {verdict}; passes={len(passes)}, "
        f"failures={len(failures)}, external_boundaries={len(boundaries)}."
    )
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
