# Systematic developmental DTU-pattern and methylation paper

This is the current folder for the manuscript and its reproducibility
materials. The `dtu_analysis/` and `methylation_analysis/` subdirectories
contain the analysis code and derived outputs, while `figures/` contains the
displayed graphical assets. The public repository is
<https://github.com/Kohze/developmental-dtu-patterns>.

This is a prepublication manuscript-and-code snapshot. See
[`RELEASE_STATUS.md`](RELEASE_STATUS.md) for the precise boundary between the
public materials and the rights-pending companion outputs.

## Scientific hierarchy

1. A general post-inference framework converts stage-specific DTU evidence and
   replicate-level transcript fractions into explicitly bounded
   diverge--reconverge candidates.
2. Applied to the mouse brain RNA-seq series, the framework resolves a strong
   E15.5-centred midbrain pattern into 1,348 candidate episodes in 735 genes.
3. Transcript structure and a fully programmatic ranking identify the six
   highest-ranked distinct reciprocal genes for independent testing.
4. A separate secondary audit of 11,002 methylation--isoform tests nominates
   `Gnao1` and `Taok3` for focused prospective work while showing that their
   support depends on methylation summarisation and expression scale.

## Accompanying R package

[`transientDTU`](https://github.com/Kohze/transientDTU) is the reusable R
implementation of the paper's post-inference decision layer. Version 0.99.1
accepts generic upstream pairwise-DTU evidence, detects bounded
diverge--reconverge episodes, checks replicate separation, annotates reciprocal
events and produces deterministic gene rankings. Its installed regression
recipe reproduced all 1,348 archived paper episodes and the exact ordered
six-gene panel. The manuscript cites the versioned package formally, and the
package documentation points users back to the article so the software and
scientific report remain mutually discoverable.

The package is developed as a separate, package-only repository so its release
history, issue tracking and Bioconductor review remain independent of the
article's analysis archive. The article release records the package version
and validation-input hashes rather than maintaining a second, potentially
stale copy of the package source.

## Citation

If these materials contribute to your work, cite the companion manuscript
using the repository's `CITATION.cff`. Analyses that use the reusable transient
DTU decision layer should additionally cite the exact `transientDTU` release.
The citation metadata will be updated with the article DOI after publication.

The comprehensive manuscript has nine main and 14 supplementary figure
environments containing 24 displayed graphical assets. The provenance register
contains 32 audited figure PDFs; eight superseded or redundant graphics remain
available but are not displayed. The main sequence adds a reproducible
framework/application overview while retaining the original thesis Figures
S16, S11, S10 and the upper panel of S9 alongside the temporal, candidate,
transcript-architecture and locus-level graphics. The old S8 and S9 panels and
the two original BH displays remain separate supplementary figures, with one
graphic per environment.

For journal review, `submission_main.tex` is a small wrapper around the
canonical `manuscript.tex`: it enables double spacing and omits the embedded
supplement without duplicating manuscript text. `supplementary_figures.tex`
builds the same S1--S14 provenance series as a separately uploadable
`Additional file 1`. The comprehensive PDF remains the archival reading copy.
The current verified builds from 24 August 2026 contain 40 pages in
`manuscript.pdf`, 33 pages in `submission_main.pdf` and 14 pages in
`supplementary_figures.pdf`.

Both build scripts can use a local untracked `tools/miktex/` toolchain and fall
back to a system MiKTeX installation when an explicit `-TexBin` is not
supplied. The large local toolchain is intentionally not part of Git history.

`prepare_journal_figures.py` creates 31 dimension-normalized PDFs in
`journal_upload_figures/` without modifying the byte-audited source figures.
`JOURNAL_FIGURE_AUDIT_2026-08-11.csv` records their final dimensions, embedded
raster resolution, sizes and SHA-256 hashes. The accompanying Markdown audit
summarizes the result. `BMC_GENOMICS_COVER_LETTER_DRAFT.md` is the current
journal cover letter with the corresponding-author and intended GitHub
repository details completed.

`prepare_release.ps1` now stages the combined article and both complete
companion analyses as one deterministic, checksummed release. The companion
contents are namespaced as `dtu_analysis/` and `methylation_analysis/` and are
selected through their existing fail-closed allow-lists. See
`PUBLIC_RELEASE_MANIFEST.md`. Until licences and source-rights confirmations
are supplied, the builder can produce an audit-only candidate but correctly
refuses a depositable release.

Key machine-readable claim and specification records are copied into
`tables/`. Complete analysis code is retained in the companion `dtu_analysis`
and `methylation_analysis` namespaces. Rights-pending generated outputs are
withheld from the public Git history until the authors clear the release gates;
the final archive must release the permitted outputs with the combined article
under one persistent identifier.

See `PUBLICATION_READINESS_AUDIT.md` for the audit verdict and
`SUBMISSION_CHECKLIST.md` for remaining gates. `PARAGRAPH_LEVEL_AUDIT.md`
records the earlier paragraph-level review and includes an addendum for the
nine prose units and twelve replacement/new captions introduced by the thesis
expansion. Superseded manuscript snapshots and dated release candidates are
not retained; `manuscript.tex` and its three built PDFs are authoritative.

`SENTENCE_LEVEL_PUBLICATION_AUDIT_2026-08-23.md` records the subsequent
sentence-by-sentence pass, the claim-calibration repairs applied to the live
source and the residual author, release, rights and reproducibility holds.

`REVIEWER_RISK_REPAIR_2026-08-10.md` records the subsequent adversarial review,
literature repair, statistical relabelling, main/supplementary restructuring
and residual submission risks.

Build with:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\build.ps1
```

Build the journal-facing main article and separate supplementary PDF with:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\build_submission.ps1
```

After resolving the licence and rights gates, build the combined repository
deposit with:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\prepare_release.ps1 `
  -Destination C:\path\to\combined-paper-release
```
