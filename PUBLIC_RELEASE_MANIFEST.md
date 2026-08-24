# Combined public-release manifest

The public repository is
<https://github.com/Kohze/developmental-dtu-patterns>. Its prepublication Git
history is subject to the exclusions in `RELEASE_STATUS.md`; the final archive
will contain one versioned release for the combined article after its licence
and source-rights gates are cleared. Its top level contains the canonical manuscript, journal-facing main
article, separate supplementary PDF, figure assets, machine-readable claim
ledger and provenance tables. The complete companion analyses are retained in
two explicit namespaces:

- `dtu_analysis/` contains the systematic transient-pattern analysis;
- `methylation_analysis/` contains the methylation association and sensitivity analyses.

The reusable decision layer is released separately as the versioned R package
[`transientDTU`](https://github.com/Kohze/transientDTU). The article records
the exact package version and paper-regression input hashes and cites the
software formally. Keeping the package in its package-only repository avoids a
stale vendored copy while preserving a bidirectional link between the article
and software releases.

Each namespace is populated from its companion `release_config.json`. The
combined builder therefore inherits the existing file allow-lists instead of
copying either working directory wholesale. `RELEASE_CONTENTS.tsv` records the
byte size and SHA-256 digest of every staged file except itself. ZIP entries
are sorted and assigned a fixed timestamp for deterministic rebuilds.

The release excludes archived `.RData` and `.rds` objects, raw reads, local
MiKTeX trees, caches, work directories, build logs, submission-planning files
and superseded manuscript drafts. The retained analysis outputs can be
released only after the authors resolve the rights status of the archived
inputs. If those objects cannot be redistributed, the deposit must retain the
recorded checksums and document the authorised-input route used by the
companion scripts.

Build from this directory with:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\prepare_release.ps1 `
  -Destination C:\path\to\combined-paper-release
```

The command fails closed unless a top-level `LICENSE` and each companion
package's `LICENSE` and `SOURCE_RIGHTS_CONFIRMED.md` are present. The
`-AllowMissingLicense` and `-AllowUnconfirmedSourceRights` switches are only
for local structural audits; any output produced with either switch is not
depositable.
