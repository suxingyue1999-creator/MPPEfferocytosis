# Public-copy adaptations

`Provenance.csv` records each source path relative to the original workspace,
the original source SHA256, the adapted-copy SHA256 and the adaptation categories.
Source files were not edited. No original result was recalculated or replaced.

Public-copy changes are limited to:

- Coded discovery/validation sample literals and removal of the original identity
  alias-normalization instructions. The private mapping is not included.
- Removal of identity-only metadata assignments and exports. Batch models read
  the supplied deidentified `batch_wave` categorical annotation directly.
- Replacement of the discovery severity inference from sample-name prefixes with
  an explicit lookup in coded input metadata. The model contrast is unchanged.
- Renaming a historical sample-specific legacy-annotation flag to
  `legacy_annotation_usable`, while retaining its existing logic on coded samples.
- Configuration or relative placeholders for machine-specific input roots, and a
  configurable Salmon executable.
- Coded sample sorting and an allowlist for discovery metadata exported by the
  audit. These changes do not alter the expression counts or group tests.

Original statistical contrasts, model formulas, expression filters, QC thresholds,
candidate genes, CellChat parameters, GO rules and reference versions are retained.
The seven-pathway leave-one-sample-out block, shared-universe cameraPR block,
PRJNA differential-expression design block and CellChat inference block were
compared with the source and found textually unchanged at preparation.

The repository does not contain the local preparation tool, private diffs,
clinical metadata, serialized objects or generated scientific result tables.
Those local preparation files are explicitly excluded by `.gitignore`.

Preparation checks parse R, Python, Bash and Perl source without executing the
scientific analyses. They also check hashes, known private literals, forbidden
identity fields, machine paths and common credential patterns. These are code
release checks, not numerical reproduction or a substitute for source-data review.
