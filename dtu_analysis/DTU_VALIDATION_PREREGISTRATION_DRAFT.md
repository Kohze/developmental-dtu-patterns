# Independent DTU validation: preregistration draft

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

- **Scg3:** NM_009130 higher and NM_001164790 lower at E15.5; priority axis: calibrated-scan rank.
- **Gpm6a:** NM_001253754 higher and NM_153581 lower at E15.5; priority axis: calibrated-scan rank.
- **Ntrk2:** NM_008745 higher and NM_001025074 lower at E15.5; priority axis: protein-domain interpretability.
- **Tecr:** NM_134118 higher and NM_027179 lower at E15.5; priority axis: cross-model support.
- **Armc8:** NM_028768 higher and NM_001166138 lower at E15.5; priority axis: cross-model support.
- **Bin1:** NM_001083334 higher and NM_009668 lower at E15.5; priority axis: effect magnitude.

The quantitative discovery values, coordinate-candidate counts and candidate-
specific counter-explanations are frozen in
`tables/preregistered_validation_panel.csv`. The coordinate audit does not
freeze primer or probe sequences because the original annotation release is
unresolved; sequence-lineage evidence is provided separately in
`tables/candidate_sequence_reconstruction_audit.csv`.

## Primary estimand and test family

For the higher-midbrain accession of gene $g$, let $L_{g,r,s}$ be the logit
of its absolute-copy-supported fraction in region $r$ and stage $s$. Define
$C_{g,r}=L_{g,r,E15.5}-[L_{g,r,E14.5}+L_{g,r,E16.5}]/2$. The primary
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
