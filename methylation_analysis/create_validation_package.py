"""Freeze analysis specifications and two locus follow-up hypotheses.

The outputs are generated from the final full-family audit so numerical values
cannot drift independently of the authoritative result table.  The protocol is
explicitly a draft until approved before any new outcome data are inspected.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
AUDIT = RESULTS / "replicate_robustness_all.csv"


def selected(frame: pd.DataFrame, column: str) -> int:
    return int(frame[column].lt(0.05).sum())


def format_float(value: float) -> str:
    return f"{value:.6g}"


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    audit = pd.read_csv(AUDIT)
    if len(audit) != 11002:
        raise ValueError(f"Expected 11,002 tests, found {len(audit):,}.")

    specifications = [
        {
            "specification_id": "S01",
            "status": "archived comparator",
            "methylation_summary": "archived quantile-normalized stage average with inherited strand-sensitive overlap",
            "rna_summary": "stage-level relative isoform fraction",
            "inferential_unit": "developmental stage",
            "statistic": "Pearson correlation with parametric t reference",
            "multiplicity_family": "separate BH family within each genomic-region class",
            "selection_rule": "archived q < 0.05",
            "selected_rows": selected(audit, "archived_q"),
            "minimum_q": format_float(audit["archived_q"].min()),
            "role": "historical result; not the primary reconstruction",
        },
        {
            "specification_id": "S02",
            "status": "primary effect / parametric screen",
            "methylation_summary": "strand-agnostic coverage-weighted CpG mean; replicate-retaining stage mean",
            "rna_summary": "replicate-retaining stage mean relative isoform fraction",
            "inferential_unit": "eight developmental stages",
            "statistic": "Pearson correlation with parametric t reference",
            "multiplicity_family": "one BH family across 11,002 fixed tests",
            "selection_rule": "global q < 0.05",
            "selected_rows": selected(audit, "raw_global_q"),
            "minimum_q": format_float(audit["raw_global_q"].min()),
            "role": "effect-ranking screen; parametric reference only",
        },
        {
            "specification_id": "S03",
            "status": "exchangeability-based sensitivity calibration",
            "methylation_summary": "same profiles as S02",
            "rna_summary": "same profiles as S02",
            "inferential_unit": "developmental stage",
            "statistic": "two-sided exhaustive stage-label randomisation under an exchangeability null",
            "multiplicity_family": "one BH family across all 10,986 valid randomisation tests",
            "selection_rule": "global randomisation q < 0.05",
            "selected_rows": selected(audit, "raw_weighted_exact_global_q"),
            "minimum_q": format_float(audit["raw_weighted_exact_global_q"].min()),
            "role": "sensitivity analysis; ordered developmental stages are not plausibly exchangeable",
        },
        {
            "specification_id": "S04",
            "status": "measurement-summary sensitivity",
            "methylation_summary": "strand-agnostic unweighted CpG mean",
            "rna_summary": "stage mean relative isoform fraction",
            "inferential_unit": "developmental stage",
            "statistic": "Pearson correlation with parametric t reference",
            "multiplicity_family": "one BH family across the fixed test set",
            "selection_rule": "global q < 0.05",
            "selected_rows": selected(audit, "unweighted_global_q"),
            "minimum_q": format_float(audit["unweighted_global_q"].min()),
            "role": "tests dependence on CpG coverage weighting",
        },
        {
            "specification_id": "S05",
            "status": "overlap-rule sensitivity",
            "methylation_summary": "coverage-weighted CpG mean with inherited strand-sensitive overlap",
            "rna_summary": "stage mean relative isoform fraction",
            "inferential_unit": "developmental stage",
            "statistic": "Pearson correlation with parametric t reference",
            "multiplicity_family": "one BH family across the fixed test set",
            "selection_rule": "global q < 0.05",
            "selected_rows": selected(audit, "strand_sensitive_global_q"),
            "minimum_q": format_float(audit["strand_sensitive_global_q"].min()),
            "role": "diagnoses the inappropriate gene-strand inheritance",
        },
        {
            "specification_id": "S06",
            "status": "temporal-shape diagnostic",
            "methylation_summary": "linear-trend residual of S02 stage profile",
            "rna_summary": "linear-trend residual of stage usage profile",
            "inferential_unit": "ordinal stage positions 0--7; P0 treated as next observation",
            "statistic": "Pearson correlation with parametric t reference",
            "multiplicity_family": "one BH family across the fixed test set",
            "selection_rule": "diagnostic global q < 0.05",
            "selected_rows": selected(audit, "detrended_global_q"),
            "minimum_q": format_float(audit["detrended_global_q"].min()),
            "role": "post hoc shape diagnostic; unequal elapsed time to P0 is not modelled; not independent confirmation",
        },
        {
            "specification_id": "S07",
            "status": "temporal-change diagnostic",
            "methylation_summary": "adjacent-stage first differences of S02 profile",
            "rna_summary": "adjacent-stage first differences of usage profile",
            "inferential_unit": "seven adjacent changes not divided by elapsed time",
            "statistic": "Pearson correlation with parametric t reference",
            "multiplicity_family": "one BH family across the fixed test set",
            "selection_rule": "diagnostic global q < 0.05",
            "selected_rows": selected(audit, "difference_global_q"),
            "minimum_q": format_float(audit["difference_global_q"].min()),
            "role": "post hoc change diagnostic; E16.5-to-P0 interval is not time-scaled; not independent confirmation",
        },
        {
            "specification_id": "S08",
            "status": "influence diagnostic",
            "methylation_summary": "S02 profile after omitting one stage",
            "rna_summary": "usage profile after omitting the same stage",
            "inferential_unit": "seven retained stages per omission",
            "statistic": "eight leave-one-stage-out Pearson correlations",
            "multiplicity_family": "none; descriptive diagnostic",
            "selection_rule": "sign stability and reported range",
            "selected_rows": "not applicable",
            "minimum_q": "not applicable",
            "role": "detects single-stage leverage",
        },
        {
            "specification_id": "S09",
            "status": "replicate-choice diagnostic",
            "methylation_summary": "one independently sampled WGBS replicate per stage",
            "rna_summary": "one independently sampled RNA replicate per stage",
            "inferential_unit": "stage; assays remain unpaired",
            "statistic": "5,000 replicate-selection correlation draws for selected union",
            "multiplicity_family": "none; conditional diagnostic",
            "selection_rule": "95% interval and sign probability",
            "selected_rows": "not applicable",
            "minimum_q": "not applicable",
            "role": "propagates replicate choice without creating new stages",
        },
        {
            "specification_id": "S10",
            "status": "scale/denominator counter-audit",
            "methylation_summary": "S02 profile",
            "rna_summary": "absolute transcript and summed gene expression",
            "inferential_unit": "developmental stage",
            "statistic": "Pearson correlation",
            "multiplicity_family": "candidate-level interpretation only",
            "selection_rule": "concordance required before mechanistic language",
            "selected_rows": "not applicable",
            "minimum_q": "not applicable",
            "role": "distinguishes isoform choice from total-expression composition",
        },
    ]
    for specification in specifications:
        specification["randomisation_selected_rows"] = "not applicable"
        specification["minimum_randomisation_q"] = "not applicable"
    randomisation_columns = {
        "S03": "raw_weighted_exact_global_q",
        "S06": "detrended_exact_global_q",
        "S07": "difference_exact_global_q",
    }
    for specification in specifications:
        column = randomisation_columns.get(specification["specification_id"])
        if column:
            specification["randomisation_selected_rows"] = selected(audit, column)
            specification["minimum_randomisation_q"] = format_float(audit[column].min())
    write_csv(RESULTS / "analysis_specification_table.csv", specifications)

    key = audit.set_index(["gene_id", "isoform_id", "region", "tissue"])
    gnao_b = key.loc[("Gnao1", "NM_001113384", "upstream", "ForeBrain")]
    gnao_a = key.loc[("Gnao1", "NM_010308", "upstream", "ForeBrain")]
    taok = key.loc[("Taok3", "NM_001199685", "gene_body", "ForeBrain")]
    assert abs(gnao_b.raw_weighted_r + gnao_a.raw_weighted_r) < 1e-12

    candidates = [
        {
            "locus": "Gnao1",
            "frozen_rna_event": "NM_001113384 versus NM_010308 alternative 3-prime coding-exon usage",
            "archived_interval": "6-kb upstream analysis window",
            "observed_direction": "higher methylation accompanies a higher NM_001113384:NM_010308 ratio; developmental loss accompanies NM_010308 and total-expression expansion",
            "primary_weighted_r": format_float(gnao_b.raw_weighted_r),
            "parametric_global_q": format_float(gnao_b.raw_global_q),
            "stage_randomisation_global_q": format_float(gnao_b.raw_weighted_exact_global_q),
            "unweighted_global_q": format_float(gnao_b.unweighted_global_q),
            "methylation_reliability": format_float(gnao_b.methylation_replicate_reliability),
            "absolute_scale_counterevidence": (
                f"NM_001113384 r={format_float(gnao_b.isoform_absolute_expression_r)}; "
                f"NM_010308 r={format_float(gnao_a.isoform_absolute_expression_r)}; "
                f"total r={format_float(gnao_b.gene_total_expression_r)}"
            ),
            "current_status": "credible developmental switch; methylation mechanism not identified",
            "decisive_primary_test": "bidirectional local methylation editing followed by measurement of the absolute NM_001113384:NM_010308 ratio and nascent terminal-exon choice",
            "falsifier": "editing changes only total or NM_010308 abundance or cell state while NM_001113384 and nascent terminal-exon choice remain unchanged",
        },
        {
            "locus": "Taok3",
            "frozen_rna_event": "NM_001199685 variant-1 first-exon/5-prime-UTR usage; validated variants encode the same protein",
            "archived_interval": "approximately 159-kb gene-body summary plus margins",
            "observed_direction": "higher weighted methylation accompanies higher variant-1 relative usage",
            "primary_weighted_r": format_float(taok.raw_weighted_r),
            "parametric_global_q": format_float(taok.raw_global_q),
            "stage_randomisation_global_q": format_float(taok.raw_weighted_exact_global_q),
            "unweighted_global_q": format_float(taok.unweighted_global_q),
            "methylation_reliability": format_float(taok.methylation_replicate_reliability),
            "absolute_scale_counterevidence": f"variant-1 absolute expression r={format_float(taok.isoform_absolute_expression_r)}",
            "current_status": "temporally coherent but coverage-weighting-sensitive transcript-start hypothesis",
            "decisive_primary_test": "independent CpG-block localization then bidirectional promoter-focused editing with variant-1 nascent initiation",
            "falsifier": "no reproducible local block or editing changes neither variant-1 initiation nor absolute/relative usage",
        },
    ]
    write_csv(RESULTS / "candidate_evidence_ledger.csv", candidates)

    protocol = """# Methylation locus-validation preregistration draft

Status: author/experimentalist/statistician review required before any new
outcome data are inspected. Gnao1 and Taok3 are frozen hypotheses, not
discoveries: exhaustive stage-label randomisation produced zero global-FDR
selections under an exchangeability null, which is a sensitivity analysis
rather than exact time-series inference.

## Frozen questions

1. **Gnao1:** does methylation at a localized block within the archived 6-kb
   upstream window causally alter NM_001113384-versus-NM_010308 terminal
   coding-exon choice, or does it change total/NM_010308 transcription while
   the NM_001113384 fraction follows the denominator?
2. **Taok3:** does methylation at an independently localized promoter/first-
   exon block alter NM_001199685 variant-1 initiation, rather than merely
   tracking a diffuse 159-kb gene-body state?

No additional locus may replace either candidate after outcomes are examined.
Any added locus belongs to a separately declared exploratory family.

## Phase 1: independent localization and replication

- Use the same purified neural lineage or a matched differentiating model for
  methylation and RNA. Bulk mixed brain tissue alone is insufficient.
- Include at least four independent biological replicates per time point; the
  final number and exclusion rules require prospective power review.
- Resolve 5mC from 5hmC with oxidative/enzymatic or equivalent chemistry.
- Measure absolute transcript copies, relative usage, total gene expression,
  nascent RNA, accessibility and cell-state/composition markers.
- Verify complete Gnao1 NM_001113384/NM_010308 and Taok3 NM_001199685
  transcript structures by targeted long reads and sequence every
  junction/first-exon assay product.
- Localize replicate-consistent CpG blocks without using the perturbation
  outcomes. Failure to identify a local block ends the corresponding causal
  experiment and is reported as a negative localization result.

## Phase 2: bidirectional perturbation

- For each surviving block, use methylation writing and erasing, at least three
  non-overlapping guides, non-targeting guides, catalytically inactive editors,
  editor-only controls and a sham manipulation.
- Randomize cultures/animals and library preparation; blind molecular outcome
  quantification where practical. The biological replicate—not a guide, well,
  read or CpG—is the inferential unit.
- Confirm on-target 5mC and 5hmC separately and assay nearby/off-target regions,
  viability, differentiation and global transcriptional state.
- Measure the modification and proposed mediator before nascent RNA, then
  mature absolute transcripts and relative usage.

## Primary estimands and multiplicity

The two locus-level primary tests form one BH family (two-sided q < 0.05), but
directional support is additionally required:

- **Gnao1:** write-versus-erase contrast in the log absolute
  NM_001113384:NM_010308 copy ratio. Writing is predicted to increase this
  ratio; erasing is predicted to decrease it.
- **Taok3:** write-versus-erase contrast in logit NM_001199685 usage, supported
  by concordant variant-1 nascent first-exon initiation. Writing is predicted
  to increase variant-1 usage.

Guide-level effects are combined in a prespecified mixed model with biological
replicate as the unit and guide as a design factor; guides are not treated as
independent biological replication. Exact covariates, transformations,
missing-data rules and the power-based sample size must be approved before
unblinding.

## Mechanistic support rule

A locus supports a methylation-to-RNA mechanism only if all are observed:

1. both editors change the intended local modification in opposite directions;
2. the primary RNA estimand changes in the prespecified opposite directions;
3. a named local mediator/Pol-II/accessibility change precedes nascent RNA;
4. mediator blockade abolishes the RNA response while editing remains intact;
5. reverse editing or mediator rescue restores the molecular trajectory; and
6. effects are not explained by total transcription, cell state, toxicity or
   off-target modification.

For Gnao1, a change confined to total expression or NM_010308 with stable
NM_001113384 and nascent terminal-exon choice specifically falsifies the
splice-choice claim.
For Taok3, failure of variant-1 nascent initiation despite successful editing
falsifies the transcript-start claim. A null or opposite result is retained and
reported; candidates are not reranked.

## Approval fields

- Final biological model and time points: [required]
- Localized CpG coordinates and guide sequences: [required after Phase 1]
- Assay/junction/first-exon sequences: [required]
- Named mediator and temporal sampling schedule: [required]
- Power calculation, sample size and exclusion rules: [statistician approval]
- Randomization/blinding plan: [required]
- Author, experimentalist and statistician signatures/dates: [required]
"""
    (HERE / "METHYLATION_VALIDATION_PREREGISTRATION_DRAFT.md").write_text(
        protocol, encoding="utf-8", newline="\n"
    )


if __name__ == "__main__":
    main()
