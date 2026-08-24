# DTU validation prospective decision record

**Status:** blank freeze template. Complete, sign, timestamp and checksum this
record before inspecting any validation outcome. Do not overwrite the signed
version; amendments require a new version with a dated rationale and must state
whether any outcome data had been accessed.

## Record identity and access boundary

- Protocol version or commit: `[required]`
- Frozen record version: `[required]`
- Freeze date and time with timezone: `[required]`
- First permitted outcome-access date: `[required]`
- Person holding blinded sample labels: `[required]`
- People who have accessed validation outcomes before freeze: `[required; use
  "none" if applicable]`
- SHA-256 of the completed signed record: `[add after signing]`

## Reference and assay freeze

- Genome assembly and exact sequence release: `[required]`
- Transcript annotation provider, release and download checksum: `[required]`
- Versioned accessions for all twelve panel transcripts: `[required]`
- Final junction/feature target for each accession: `[attach versioned table]`
- Primer/probe sequences, chemistry and supplier identifiers: `[attach]`
- Genome/transcriptome-wide specificity report: `[attach]`
- Thermodynamic and secondary-structure acceptance criteria: `[required]`
- Efficiency, limit-of-detection and zero-call acceptance criteria: `[required]`
- Product-sequencing and long-read confirmation plan: `[required]`
- Lower Armc8 NM_001166138 terminal-segment assay approved: `[yes/no]`

The archived coordinate and 40-nt k-mer tables are lineage anchors only. They
must not be copied into this section as if they were validated oligos.

## Biological design freeze

- Independent unit: `[prefer independent litter]`
- Prespecified embryo-selection rule within each litter: `[required]`
- Regions: `[forebrain, midbrain, hindbrain required for primary estimand]`
- Stages: `[E14.5, E15.5, E16.5 required]`
- Sex ascertainment and planned use: `[required]`
- Dissection, extraction, plate and processing blocks: `[required]`
- Randomisation procedure: `[required]`
- Blinding procedure and unblinding trigger: `[required]`
- Permitted replacements before outcome access: `[required]`
- Prohibited replacements after outcome access: `[required]`

## Primary estimand and analysis freeze

- Primary estimand: `Delta_g = C_g,midbrain -
  0.5*(C_g,forebrain + C_g,hindbrain)` where each `C` is E15.5 minus the
  E14.5/E16.5 mean on the frozen logit-usage scale.
- Primary genes/accessions: `[attach the frozen six-row panel]`
- Model family and link: `[required]`
- Zero and boundary handling: `[required]`
- Litter/embryo blocking or random-effects structure: `[required]`
- Variance estimator and small-sample correction: `[required]`
- Permitted covariates and their coding: `[required]`
- Primary multiplicity rule: `[BH across six higher-accession tests unless
  prospectively replaced]`
- Direction requirement: `[Delta_g > 0 required]`
- Interval level and reporting policy: `[required; report every gene]`
- Software, package versions and executable analysis checksum: `[required]`

## Power and recruitment decision

Record one exact row from each released planning table or attach a justified
replacement calculation. Discovery-effect maximisation is not a valid source.

- Smallest meaningful absolute logit-curvature effect: `[required]`
- Outcome-independent marginal SD or standardized effect: `[required]`
- Source for effect/variance choice: `[external, blinded feasibility or
  conservative bound; cite/attach]`
- Within-embryo cross-region correlation: `[required]`
- Between-gene test-statistic correlation: `[required]`
- Per-gene assay-failure probability: `[required]`
- Target joint-panel power: `[required]`
- Selected row key from
  `tables/validation_joint_panel_power_sensitivity.csv`: `[required]`
- Target analyzable independent units per stage: `[required]`
- Whole-unit loss probability: `[required]`
- Retention-assurance target: `[required]`
- Selected row key from `tables/validation_collection_inflation.csv`:
  `[required]`
- Units to recruit per stage and total across stages: `[required]`
- Cluster/design-effect calculation if more than one embryo per litter is
  analysed: `[required or explicitly not applicable]`

## Panel rule and secondary analyses

- Accept the proposed panel rule unchanged: `[yes/no]`
- If no, replacement rule and rationale: `[required before data access]`
- Proposed rule: at least four of six genes significant in the expected
  direction, including at least one of the two scan-led genes and one of the
  two cross-model genes, with no opposite-direction significant gene.
- Reciprocal-accession checks: `[required directional secondary analysis]`
- Absolute-copy/total-gene-expression checks: `[required]`
- Long-read full-length confirmation: `[required]`
- Cell-composition measurements or matched cell-type plan: `[required]`
- Any chromatin analysis: `[secondary and conditional on RNA replication]`

## Missingness, exclusions and deviations

- Sample-level QC metrics and thresholds: `[required]`
- Region-level missingness rule: `[required]`
- Gene/assay-level non-evaluable rule: `[required]`
- Outlier rule independent of candidate outcomes: `[required]`
- Handling of incomplete matched-region triplets: `[required]`
- Permitted sensitivity analyses: `[required]`
- Deviation log location and custodian: `[required]`

All exclusions, failures and deviations must be reported by randomized unit,
stage, region and gene. The primary analysis must not be silently replaced by
the sensitivity analysis that gives the most favourable result.

## Approval

- Lead author, signature/date: `[required]`
- Experimental lead, signature/date: `[required]`
- Statistician, signature/date: `[required]`
- Data custodian confirming outcomes remained inaccessible, signature/date:
  `[required]`

After approval, render this record read-only, calculate its SHA-256, register it
with the protocol and analysis code, and retain the blank template separately.
