# Reproduction environment

The archived-object analyses were finalized under R 4.4.0. Exact R,
Bioconductor, and package versions are recorded in `data/sessionInfo.txt`. The
candidate assay-target audit uses `BSgenome.Mmusculus.UCSC.mm10` 1.4.3, the
reference package identifier retained in both archived workspaces. The
accession-architecture and claim-ledger scripts use Python 3.11.5 with packages
already available in the analysis environment. The manuscript was compiled
with pdfTeX/MiKTeX 25.12 using embedded Latin Modern vector fonts.
`build.ps1` fixes `SOURCE_DATE_EPOCH` to 1 January 2000 and enables
`FORCE_SOURCE_DATE`, making repeated builds with the same engine and inputs
byte-identical rather than embedding the wall-clock build time.

`build.ps1 -RunAnalysis -InputDir <directory>` runs the complete R analysis
when the two archived workspaces are available. Their immutable hashes and the
unresolved annotation/provenance boundary are recorded under `tables/`.
