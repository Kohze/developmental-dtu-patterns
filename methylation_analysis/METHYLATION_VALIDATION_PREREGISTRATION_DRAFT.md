# Methylation locus-validation preregistration draft

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
