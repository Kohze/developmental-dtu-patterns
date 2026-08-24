# Journal figure dimension and raster audit

Date: 13 August 2026 (re-audited after main-figure revision)

All 32 audited figure PDFs were exported as separate
journal-upload derivatives. Each derivative is at most 170 mm wide and
225 mm high, with aspect ratio and vector content preserved. The source
files referenced by `figure_provenance.csv` were not modified.

- Vector-only derivatives: 24
- Derivatives containing directly placed raster components: 8
- Derivatives with a raster component below 300 effective dpi: 0
- Exact dimensions, raster placements and SHA-256 hashes are recorded in
  `JOURNAL_FIGURE_AUDIT_2026-08-11.csv`.

Every directly placed raster component is at least 300 effective
dpi at the exported dimensions.

## Reproduction

Run `python prepare_journal_figures.py` from this directory. The script
recreates every derivative and replaces the manifest with current
dimensions and hashes.
