"""Create a machine-readable DTU claim ledger and validation-panel draft.

The outputs deliberately separate archive-supported claims from statements
that require new biological material. No new inference is performed here; all
numbers are read from the manuscript's generated audit tables.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
TABLES = PAPER / "tables"


def scalar_table(path: Path) -> dict[str, float]:
    frame = pd.read_csv(path)
    return dict(zip(frame.iloc[:, 0], frame.iloc[:, 1]))


def build_claim_ledger() -> pd.DataFrame:
    numbers = scalar_table(TABLES / "manuscript_numbers.csv")
    counter = scalar_table(TABLES / "counter_audit_summary.csv")
    scan = scalar_table(TABLES / "transient_regional_scan_diagnostics.csv")
    joint = scalar_table(TABLES / "joint_interaction_summary.csv")
    artifact = pd.read_csv(TABLES / "e155_episode_artifact_tests.csv").set_index("test")
    adjacent = pd.read_csv(TABLES / "adjacent_stage_counts.csv")
    replicate_choice_path = TABLES / "candidate_replicate_choice_audit.csv"
    replicate_choice = (
        pd.read_csv(replicate_choice_path) if replicate_choice_path.exists() else None
    )
    assay_targets = pd.read_csv(TABLES / "candidate_assay_target_summary.csv")
    if len(assay_targets) != 12:
        raise ValueError("Candidate assay-target summary must contain twelve accessions.")
    sequence_audit = pd.read_csv(TABLES / "candidate_sequence_reconstruction_audit.csv")
    if len(sequence_audit) != 12:
        raise ValueError("Candidate sequence audit must contain twelve accessions.")
    dependence = pd.read_csv(TABLES / "transient_regional_dependence_sensitivity.csv")
    conservative = dependence.loc[
        dependence["method"].eq("arbitrary-dependence sensitivity")
    ].iloc[0]

    def adjacent_count(early: str, late: str) -> int:
        row = adjacent[
            adjacent["tissue"].eq("Midbrain")
            & adjacent["early"].astype(str).eq(early)
            & adjacent["late"].astype(str).eq(late)
        ]
        return int(row.iloc[0]["genes"])

    wilcoxon = artifact.loc[
        "E15.5 higher-versus-lower isoform length (Wilcoxon rank-sum)"
    ]
    sign_test = artifact.loc[
        "Reciprocal genes with higher isoforms shorter than lower isoforms (sign test)"
    ]
    records = [
        {
            "claim_id": "DTU-C01",
            "manuscript_claim": "The surviving archive covers three brain regions at eight stages with two biological replicates per region-stage cell.",
            "evidence_type": "design inventory",
            "numerical_result": f"{int(numbers['Total libraries'])} samples; 24 region-stage cells; 2 replicates per cell",
            "primary_source_files": "tables/manuscript_numbers.csv; tables/encode_library_provenance.csv",
            "inferential_status": "descriptive and archive-supported",
            "live_alternative_or_limit": "Archived sample suffixes are not linked to portal replicate numbers.",
            "decisive_next_test": "Recover the original run manifest and raw-input-to-matrix mapping.",
        },
        {
            "claim_id": "DTU-C02",
            "manuscript_claim": "Most detected DTU genes belong to a three-region shared core.",
            "evidence_type": "primary archive reanalysis",
            "numerical_result": f"{int(counter['Primary DTU union'])} union genes; {int(counter['Primary shared core'])} shared ({100 * counter['Primary shared core'] / counter['Primary DTU union']:.1f}%)",
            "primary_source_files": "tables/dtu_membership_counts.csv; tables/counter_audit_summary.csv",
            "inferential_status": "archive-supported detection pattern",
            "live_alternative_or_limit": "Membership depends on within-contrast thresholds and is not strict tissue specificity.",
            "decisive_next_test": "Re-estimate from raw counts under one joint region-by-stage model with batch covariates.",
        },
        {
            "claim_id": "DTU-C03",
            "manuscript_claim": "The dominant temporal anomaly is an E15.5-centred midbrain peak.",
            "evidence_type": "adjacent-stage contrast inventory",
            "numerical_result": f"{adjacent_count('14.5', '15.5')} genes at E14.5-E15.5; {adjacent_count('15.5', '16.5')} at E15.5-E16.5",
            "primary_source_files": "tables/adjacent_stage_counts.csv; tables/audit_threshold_sensitivity.csv",
            "inferential_status": "robust within the archive, not independently replicated",
            "live_alternative_or_limit": "E14.5 crosses a collection-source boundary; E15.5/E16.5 are pooled tissues; cell composition and coordinated technical effects remain possible.",
            "decisive_next_test": "Independent, balanced forebrain/midbrain/hindbrain collection at E14.5/E15.5/E16.5 processed in one randomized design.",
        },
        {
            "claim_id": "DTU-C04",
            "manuscript_claim": "A globally calibrated scan resolves the peak into diverge-reconverge isoform episodes.",
            "evidence_type": f"global {int(scan['pairwise_tests']):,}-test transient scan",
            "numerical_result": f"{int(scan['primary_episodes'])} episodes in {int(scan['primary_episode_genes'])} genes; arbitrary-dependence sensitivity retains {int(conservative['replicate_separated_episodes'])} episodes in {int(conservative['genes'])} genes and the same leading six",
            "primary_source_files": "data/transient_regional_filtered_isoform_fractions.csv; data/transient_regional_pair_tests_all.csv; data/transient_regional_stage_evaluations_all.csv; data/transient_regional_isoform_episodes.csv; data/transient_regional_isoform_episodes_conservative.csv; tables/transient_regional_scan_diagnostics.csv; tables/transient_regional_dependence_sensitivity.csv",
            "inferential_status": "discovery analysis with complete-family component-test adjustment and arbitrary-dependence sensitivity; no episode-level FDR claim",
            "live_alternative_or_limit": "The complete expression-filtered family is corrected, but the trajectory rules were formulated after archive inspection.",
            "decisive_next_test": "Apply the locked scan and candidate panel to an independent cohort.",
        },
        {
            "claim_id": "DTU-C05",
            "manuscript_claim": "The calibrated episodes are concentrated in single-stage E15.5 midbrain events.",
            "evidence_type": "episode classification",
            "numerical_result": f"{int(scan['single_stage_e15_5_episodes'])}/{int(scan['primary_episodes'])} episodes ({100 * scan['single_stage_e15_5_episodes'] / scan['primary_episodes']:.1f}%); {int(scan['non_e15_5_episodes'])} other episodes",
            "primary_source_files": "tables/transient_regional_episode_summary.csv; tables/transient_regional_scan_diagnostics.csv",
            "inferential_status": "archive-supported concentration",
            "live_alternative_or_limit": "Concentration can arise from one stage-specific batch or composition shift.",
            "decisive_next_test": "Replication across independent litters with stage labels concealed during processing and quantification.",
        },
        {
            "claim_id": "DTU-C06",
            "manuscript_claim": "Simple transcript-end dominance or a systematic shift toward shorter products does not explain the E15.5 episode set.",
            "evidence_type": "post hoc structural counter-audit",
            "numerical_result": f"higher-versus-lower median length difference {wilcoxon['estimate']:.0f} nt, P={wilcoxon['p_value']:.3f}; reciprocal-gene shorter fraction {sign_test['estimate']:.3f}, P={sign_test['p_value']:.3f}",
            "primary_source_files": "tables/e155_episode_event_bias.csv; tables/e155_episode_artifact_tests.csv",
            "inferential_status": "negative control against two simple artefact classes",
            "live_alternative_or_limit": "Does not exclude mapping ambiguity, read-position bias, reference error or sequence-dependent quantification.",
            "decisive_next_test": "Junction-specific counting plus full-length long-read confirmation from new RNA.",
        },
        {
            "claim_id": "DTU-C07",
            "manuscript_claim": "A focused region-by-stage curvature model supports an E15.5-specific archive signal.",
            "evidence_type": "post hoc satuRn sensitivity model",
            "numerical_result": f"{int(joint['genes_dexseq_style_per_gene_q_lt_0_05'])}/{int(joint['multi_isoform_genes_tested'])} genes at per-gene q<0.05; {int(joint['transcripts_empirical_fdr_lt_0_05'])} transcripts at empirical FDR<0.05",
            "primary_source_files": "data/joint_interaction_genes.csv; data/joint_interaction_isoforms.csv; tables/joint_interaction_summary.csv",
            "inferential_status": "concordant sensitivity analysis, not validation",
            "live_alternative_or_limit": "Uses the same abundance-scaled archive and a contrast chosen after seeing the peak.",
            "decisive_next_test": "Fit a prespecified model to raw independent counts with batch and composition terms.",
        },
        {
            "claim_id": "DTU-C08",
            "manuscript_claim": "The six-gene panel deliberately spans distinct evidence axes rather than forming one composite ranking.",
            "evidence_type": "candidate counter-audit",
            "numerical_result": "Scg3/Gpm6a: scan rank; Armc8/Tecr: cross-model support; Bin1: effect size; Ntrk2: domain interpretation",
            "primary_source_files": "tables/candidate_counter_audit.csv; tables/candidate_mechanism_crosswalk.csv",
            "inferential_status": "prespecified validation panel derived from discovery data",
            "live_alternative_or_limit": "All candidates were selected after inspecting the same archive.",
            "decisive_next_test": "Test all six without reranking in an independent three-region by three-stage experiment.",
        },
        {
            "claim_id": "DTU-C09",
            "manuscript_claim": "The archive does not establish a biological developmental programme or mechanism.",
            "evidence_type": "claim boundary",
            "numerical_result": "two replicates per cell; bulk short-read RNA; unresolved dissection/pool batch; no independent cohort",
            "primary_source_files": "tables/encode_midbrain_provenance.csv; tables/encode_e155_pool_audit.csv; tables/input_provenance_audit.csv",
            "inferential_status": "not established",
            "live_alternative_or_limit": "Changing neural, glial or vascular composition and coordinated technical effects remain live explanations.",
            "decisive_next_test": "Independent junction and long-read replication followed by sorted or single-cell-type measurement.",
        },
    ]
    if replicate_choice is not None:
        expected = replicate_choice["joint_logit_expected_direction_fraction"]
        combinations = replicate_choice["replicate_choice_combinations"]
        records.append(
            {
                "claim_id": "DTU-C10",
                "manuscript_claim": "The frozen accession-pair directions do not depend on selecting one archived replicate suffix over the other.",
                "evidence_type": "exhaustive replicate-choice sensitivity audit",
                "numerical_result": f"all {len(replicate_choice)} pairs retain reciprocal expected directions in {int(combinations.min())}/{int(combinations.min())} one-replicate-per-cell combinations",
                "primary_source_files": "tables/candidate_replicate_choice_audit.csv",
                "inferential_status": "exhaustive archive sensitivity analysis, not independent validation",
                "live_alternative_or_limit": "All combinations reuse the same archived samples, and suffixes are not linked to portal biological-replicate numbers.",
                "decisive_next_test": "Prespecified three-region by three-stage replication in independent embryos/litters.",
            }
        )
        if not (expected.eq(1).all() and combinations.eq(512).all()):
            raise ValueError("Replicate-choice audit no longer supports the frozen claim.")
    junction_ready = assay_targets["preferred_junction_coordinate_candidates"].gt(0)
    armc8_lower = assay_targets.loc[
        assay_targets["gene"].eq("Armc8")
        & assay_targets["expected_direction"].eq("lower")
    ]
    if not (
        junction_ready.sum() == 11
        and len(armc8_lower) == 1
        and armc8_lower["preferred_junction_coordinate_candidates"].iat[0] == 0
        and armc8_lower["longest_gene_unique_segment_nt"].iat[0] == 2242
    ):
        raise ValueError("The frozen assay-target feasibility result has changed.")
    records.append(
        {
            "claim_id": "DTU-C11",
            "manuscript_claim": "The frozen panel is coordinate-auditable, but a generic two-junction assay is not feasible for every accession pair.",
            "evidence_type": "strand-aware archived exon/junction audit",
            "numerical_result": "11/12 accessions have an archived-gene-unique junction; lower Armc8 NM_001166138 has none but has a 2,242-nt gene-unique 3-prime terminal-exon segment",
            "primary_source_files": "tables/candidate_assay_target_summary.csv; tables/candidate_junction_target_audit.csv; tables/candidate_unique_exonic_segments.csv",
            "inferential_status": "provisional coordinate feasibility, not a frozen assay design",
            "live_alternative_or_limit": "The archived exon object lacks genome and annotation-version metadata; coordinates, sequence and specificity can change after annotation freeze.",
            "decisive_next_test": "Freeze the prospective reference/annotation, revalidate coordinates and off-targets, then sequence-verify every amplicon and confirm complete structures by long reads.",
        }
    )
    surviving_sequences = sequence_audit["archived_sequence_survives"].astype(bool)
    exact_surviving = sequence_audit.loc[
        surviving_sequences, "exact_match_to_archived_sequence"
    ].astype(bool)
    preferred_sequence_targets = int(
        assay_targets["preferred_junction_sequence_candidates"].sum()
    )
    if not (
        surviving_sequences.sum() == 10
        and exact_surviving.all()
        and sequence_audit["reconstructed_ambiguous_bases"].eq(0).all()
        and preferred_sequence_targets == 34
    ):
        raise ValueError("The frozen sequence-reconstruction result has changed.")
    records.append(
        {
            "claim_id": "DTU-C12",
            "manuscript_claim": "Archived exon coordinates and the named mm10 reference reproduce the surviving panel transcript sequences exactly and support provisional junction-spanning sequence candidates.",
            "evidence_type": "reference-linked transcript and junction-sequence reconstruction",
            "numerical_result": "10/10 surviving panel ntSequences match exactly; both Scg3 sequences are reconstructed but absent from the stored candidate sequence set; all 12 reconstructions contain zero ambiguous bases; 34 preferred junctions have 20+20-nt k-mers unique within archived same-gene models",
            "primary_source_files": "tables/candidate_sequence_reconstruction_audit.csv; tables/candidate_junction_target_audit.csv",
            "inferential_status": "sequence-lineage and within-gene specificity audit, not primer/probe validation",
            "live_alternative_or_limit": "No genome-wide primer-pair/off-target, secondary-structure or wet-lab efficiency test has been performed, and the original RefSeq annotation version is unresolved.",
            "decisive_next_test": "Freeze the prospective annotation/reference, design primers and probes, perform genome/transcriptome-wide in-silico specificity checks, and sequence-verify products on independent RNA.",
        }
    )
    return pd.DataFrame(records)


def build_validation_panel() -> pd.DataFrame:
    candidates = pd.read_csv(TABLES / "candidate_counter_audit.csv")
    crosswalk = pd.read_csv(TABLES / "candidate_mechanism_crosswalk.csv")
    top = pd.read_csv(TABLES / "transient_regional_top_candidates.csv")
    targets = pd.read_csv(TABLES / "candidate_assay_target_summary.csv")
    if len(targets) != 12 or targets[["gene", "expected_direction"]].duplicated().any():
        raise ValueError("Candidate assay-target summary must contain one row per accession direction.")
    higher_targets = targets.loc[
        targets["expected_direction"].eq("higher"),
        [
            "gene",
            "preferred_junction_coordinate_candidates",
            "preferred_junction_sequence_candidates",
            "longest_gene_unique_segment_nt",
        ],
    ].rename(
        columns={
            "preferred_junction_coordinate_candidates": "higher_preferred_junction_candidates",
            "preferred_junction_sequence_candidates": "higher_preferred_junction_sequence_candidates",
            "longest_gene_unique_segment_nt": "higher_longest_gene_unique_segment_nt",
        }
    )
    lower_targets = targets.loc[
        targets["expected_direction"].eq("lower"),
        [
            "gene",
            "preferred_junction_coordinate_candidates",
            "preferred_junction_sequence_candidates",
            "longest_gene_unique_segment_nt",
        ],
    ].rename(
        columns={
            "preferred_junction_coordinate_candidates": "lower_preferred_junction_candidates",
            "preferred_junction_sequence_candidates": "lower_preferred_junction_sequence_candidates",
            "longest_gene_unique_segment_nt": "lower_longest_gene_unique_segment_nt",
        }
    )
    panel = candidates.merge(crosswalk, on="gene", validate="one_to_one").merge(
        top[["gene_id", "max_worst_pair_q"]],
        left_on="gene",
        right_on="gene_id",
        validate="one_to_one",
    ).merge(higher_targets, on="gene", validate="one_to_one").merge(
        lower_targets, on="gene", validate="one_to_one"
    )
    priority = {
        "Scg3": "calibrated-scan rank",
        "Gpm6a": "calibrated-scan rank",
        "Armc8": "cross-model support",
        "Tecr": "cross-model support",
        "Bin1": "effect magnitude",
        "Ntrk2": "protein-domain interpretability",
    }
    panel["priority_axis"] = panel["gene"].map(priority)
    if not panel["higher_preferred_junction_candidates"].gt(0).all():
        raise ValueError("Every higher-accession primary endpoint requires a junction candidate.")
    if not panel["higher_preferred_junction_sequence_candidates"].gt(0).all():
        raise ValueError("Every higher-accession endpoint requires a junction sequence candidate.")
    armc8_lower = panel.loc[panel["gene"].eq("Armc8")]
    if not (
        len(armc8_lower) == 1
        and armc8_lower["lower_preferred_junction_candidates"].iat[0] == 0
        and armc8_lower["lower_longest_gene_unique_segment_nt"].iat[0] >= 20
        and panel.loc[~panel["gene"].eq("Armc8"), "lower_preferred_junction_candidates"].gt(0).all()
    ):
        raise ValueError("The frozen Armc8 reciprocal-assay exception has changed.")
    panel["primary_assay"] = panel.apply(
        lambda row: (
            "archived-gene-unique junction with a within-gene-discriminating 40-nt k-mer for the higher accession; "
            + (
                "reciprocal lower accession requires a gene-unique terminal-exon segment assay"
                if row.gene == "Armc8"
                else "archived-gene-unique junction with a within-gene-discriminating 40-nt k-mer for the reciprocal lower accession"
            )
            + "; freeze reference and verify sequence/specificity before oligo design"
        ),
        axis=1,
    )
    panel["primary_endpoint"] = (
        "midbrain E15.5 logit-usage curvature minus the mean corresponding forebrain and hindbrain curvature"
    )
    panel["expected_direction"] = panel.apply(
        lambda row: (
            f"{row.higher_midbrain_accession} positive regional-curvature contrast; "
            f"{row.lower_midbrain_accession} reciprocal negative regional-curvature contrast"
        ),
        axis=1,
    )
    panel["primary_multiplicity_family"] = (
        "six two-sided higher-accession curvature tests; BH q<0.05 plus prespecified positive direction"
    )
    panel["secondary_confirmation"] = panel["gene"].map(
        lambda gene: (
            "reciprocal gene-unique terminal-exon-segment direction, absolute copy number, "
            "total gene expression and full-length long-read structure"
            if gene == "Armc8"
            else "reciprocal junction direction, absolute copy number, total gene expression "
            "and full-length long-read structure"
        )
    )
    columns = [
        "gene",
        "higher_midbrain_accession",
        "lower_midbrain_accession",
        "priority_axis",
        "calibrated_scan_rank",
        "max_worst_pair_q",
        "max_abs_usage_difference",
        "min_replicate_separation",
        "joint_model_rank_of_4577",
        "joint_gene_q_value",
        "higher_preferred_junction_candidates",
        "lower_preferred_junction_candidates",
        "higher_preferred_junction_sequence_candidates",
        "lower_preferred_junction_sequence_candidates",
        "lower_longest_gene_unique_segment_nt",
        "primary_assay",
        "primary_endpoint",
        "expected_direction",
        "primary_multiplicity_family",
        "secondary_confirmation",
        "principal_counter_explanation",
    ]
    return panel[columns].sort_values("calibrated_scan_rank")


def write_preregistration(panel: pd.DataFrame) -> None:
    gene_lines = "\n".join(
        f"- **{row.gene}:** {row.higher_midbrain_accession} higher and "
        f"{row.lower_midbrain_accession} lower at E15.5; priority axis: "
        f"{row.priority_axis}."
        for row in panel.itertuples(index=False)
    )
    text = f"""# Independent DTU validation: preregistration draft

**Status:** analysis-ready draft; it must be approved by the authors,
experimentalist and statistician before any new outcome data are inspected.

## Confirmatory question

Does an independently collected three-region embryonic-brain cohort reproduce
the prespecified excess E15.5 midbrain curvature of the six accession pairs
selected from the archive, when all regions and stages are collected and
processed in one balanced design?

## Design that must be fixed before collection

- E14.5, E15.5 and E16.5 forebrain, midbrain and hindbrain, with matched regions
  collected from each embryo where feasible. Use at least six independent
  litters per stage with one prespecified embryo per litter where possible; any
  smaller or larger final number requires a prospective power calculation
  approved before outcome data are inspected.
- Treat independent litters, not assays, technical replicates or unmodelled
  litter-mate embryos, as the planning units. Balance litter, sex where
  measurable, dissection operator, extraction
  batch, library/plate position and processing order across region--stage cells;
  conceal region and stage labels from assay operators where practical.
- Record exclusions before molecular outcomes are inspected. Do not replace a
  failed sample selectively by gene or assay.
- Use the archived-coordinate candidates in
  `tables/candidate_junction_target_audit.csv`, but freeze the exact reference
  and annotation and revalidate every sequence and off-target before ordering
  oligos. Eleven of twelve accessions have at least one splice junction unique
  among the archived isoforms of that gene. The exception is the lower Armc8
  accession NM_001166138, whose splice junctions are all shared with its paired
  accession; its reciprocal assay must instead target the audited gene-unique
  terminal-exon segment and be confirmed by full-length reads. Measure absolute
  target copies, total gene expression and a stable reference panel; do not
  rely only on compositional percentages.
- The named mm10 reference and archived exon models exactly reconstruct all ten
  panel ntSequences that survive in the candidate object; the two Scg3
  sequences are reconstructed from the same reference but are not retained in
  that sequence set. All 34 preferred junction coordinates yield provisional
  20+20-nt junction-spanning k-mers unique within the archived same-gene
  transcript models. These k-mers are audit anchors, not primers or probes;
  genome/transcriptome-wide off-target, thermodynamic and efficiency checks
  remain mandatory after the reference freeze.
- Obtain targeted full-length long reads on the same RNA or a prespecified
  subset to confirm that the assayed junctions belong to the claimed complete
  transcripts.

## Frozen six-gene panel

{gene_lines}

The quantitative discovery values, coordinate-candidate counts and candidate-
specific counter-explanations are frozen in
`tables/preregistered_validation_panel.csv`. The coordinate audit does not
freeze primer or probe sequences because the original annotation release is
unresolved; sequence-lineage evidence is provided separately in
`tables/candidate_sequence_reconstruction_audit.csv`.

## Primary estimand and test family

For the higher-midbrain accession of gene $g$, let $L_{{g,r,s}}$ be the logit
of its absolute-copy-supported fraction in region $r$ and stage $s$. Define
$C_{{g,r}}=L_{{g,r,E15.5}}-[L_{{g,r,E14.5}}+L_{{g,r,E16.5}}]/2$. The primary
estimand is

`Delta_g = C_g,midbrain - 0.5*(C_g,forebrain + C_g,hindbrain)`.

Fit one prespecified replicate-level model per gene with fixed region--stage
cells and a litter/embryo blocking term when regions are matched. Freeze the
zero handling, model family, variance estimator, permitted design covariates
and power calculation before unblinding. Test the six higher-accession regional-
curvature contrasts two-sided and control Benjamini--Hochberg FDR across those
six tests. A gene replicates only when q<0.05 and `Delta_g > 0`. The corresponding
within-midbrain curvature and reciprocal lower-accession regional contrast are
required directional secondary checks, not additional primary tests. Report
estimates and intervals for every gene regardless of significance.

## Prospective power decision

`DTU_VALIDATION_POWER_GUIDE.md` derives the variance of the matched three-region
curvature contrast, and `tables/validation_power_sensitivity.csv` supplies 40
two-sided Bonferroni marginal-power scenarios.
`tables/validation_joint_panel_power_sensitivity.csv` supplies 1,920 scenarios
for the complete suggested decision rule: at least four of six genes in the
expected direction, including at least one scan-led and one cross-model gene,
with no opposite-direction significant gene. It crosses standardized effect,
within-embryo region correlation, between-gene test correlation and per-gene
assay failure. `tables/validation_collection_inflation.csv` gives exact binomial
recruitment inflation for whole-unit loss. Before unblinding, replace these
sensitivity grids with a signed decision record that freezes the smallest
meaningful effect, an outcome-independent covariance estimate, litter
clustering, assay-failure and whole-unit-loss assumptions, and the selected
joint-panel scenario.
Use `DTU_VALIDATION_DECISION_RECORD_TEMPLATE.md` for that freeze; retain its
signed checksum with the protocol and executable analysis version.

## Suggested panel-level decision rule

Treat the discontinuity panel as independently reproduced only if at least four
of six genes replicate, including at least one scan-led candidate (Scg3 or
Gpm6a) and at least one cross-model candidate (Armc8 or Tecr), and no candidate
is significant in the opposite direction. This rule is deliberately demanding
and must be accepted or replaced before data collection; it must not be tuned
after results are known.

## Falsification and interpretation

- Fewer than four reproducing genes, opposite-direction findings, or a signal
  explained by a recorded batch covariate falsifies the panel-level replication
  criterion.
- Junction or Armc8 terminal-segment replication without full-length
  confirmation supports a local transcript feature but not the archived
  complete transcript structure.
- Midbrain-only temporal curvature does not establish the archived regional
  divergence and therefore cannot satisfy the primary replication criterion.
- Bulk regional replication does not establish cell-intrinsic regulation. The
  replicated junctions must next be measured in matched neural, glial and vascular
  populations, or in a design that estimates those proportions.
- Even complete replication establishes a developmental isoform discontinuity,
  not its chromatin mechanism. Any H3K36me3/H3K4me1 analysis is secondary until
  the RNA event is independently confirmed.
"""
    (PAPER / "DTU_VALIDATION_PREREGISTRATION_DRAFT.md").write_text(
        text, encoding="utf-8"
    )


def main() -> None:
    ledger = build_claim_ledger()
    panel = build_validation_panel()
    ledger.to_csv(TABLES / "claim_evidence_ledger.csv", index=False)
    panel.to_csv(TABLES / "preregistered_validation_panel.csv", index=False)
    write_preregistration(panel)
    print(f"Wrote {len(ledger)} claim rows and {len(panel)} validation genes.")


if __name__ == "__main__":
    main()
