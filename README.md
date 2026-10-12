# MPP Efferocytosis

Recovered analysis source code for the study of macrophage and monocyte
efferocytosis programs in pediatric Mycoplasma pneumoniae pneumonia.

The code was recovered from the existing analysis workspace and historical
delivery snapshots. Public copies use coded sample identifiers and configurable
input roots; identity-only metadata and exact clinical dates were removed.
The statistical methods have not been replaced with newly written analyses.
Original and adapted SHA256 values are listed in [Provenance.csv](Provenance.csv).
The repository contains source code and documentation only, with no clinical
records, sample metadata, sequencing files, mass-spectrometry files, results,
serialized single-cell objects, or private sample mappings.

| Directory | Analysis covered |
|---|---|
| `scripts/Discovery` | 10-sample count-matrix audit, Scrublet, sample-specific QC, SCTransform/RPCA, annotation, pseudobulk differential expression and QC sensitivity |
| `scripts/Validation` | 23-sample RNA QC, LogNormalize/RPCA, lineage analyses, pseudobulk differential expression, gene-set tests, V(D)J and sample-specific CellChat |
| `scripts/Cross cohort` | Harmonized lineages, shared-universe cameraPR, seven replicated GO processes, standardized pathway scores, fixed/random effects and myeloid-state composition |
| `scripts/Sensitivity and figures` | The actual seven-pathway leave-one-sample-out refits in three models, candidate coverage audit and historical main-figure generation |
| `scripts/Public datasets` | GSE274269, GSE179051 and PRJNA792730 analyses, FASTQ audit, Salmon quantification, tables and figures |
| `scripts/Result imports` | Extraction of the existing discovery cell-level differential-expression tables, precomputed gene-set scores and proteomic differential tables |
| `scripts/Tables and model figure` | Historical supplementary-table and candidate-model figure generation |

See [Inputs.md](Inputs.md) for required input contracts and execution order,
[Environment.md](Environment.md) for stage-specific software, and
[Adaptations.md](Adaptations.md) for the scope of the public-copy changes.

Reproducibility boundaries:

- The archived discovery/cross-cohort computations and the validation/public
  computations used different R environments. A single common version is not
  claimed; recorded stage versions are documented separately.
- The four discovery gene-set scores were imported from a precomputed
  cell-type summary. The source implementation of those scores was not recovered.
- Discovery cell-level candidate effects were imported from pre-existing
  differential-expression spreadsheets and are distinct from the newer count-based
  pseudobulk analyses.
- The reported myeloid-state ordering used pre-existing trajectory outputs.
  No recovered script in this release computes that original trajectory.
- Proteomic scripts import delivered differential tables. They do not recreate
  the raw DIA-NN processing or the original protein-level Welch tests. The original
  execution script for the subsequent full-table BH recalculation was not located;
  no substitute implementation is presented as original code.
- Figure scripts document historical result assembly. Later manuscript-layout
  edits and all final exported image versions are not reproduced by these scripts.
- Validation scripts also contain the originally analysed plastic-bronchitis and
  allergy contrasts. This manuscript uses the severity contrast and the specified
  covariate adjustments; the other contrasts are retained as source provenance.

Verification at preparation was limited to source/release checksums, syntax
parsing, code/privacy scans and comparison of selected core statistical blocks.
The public copies have not been rerun end to end, and no new scientific results
were produced. Historical freeze checks require the corresponding audited input
objects and delivery manifests, which are not included here. No synthetic
replacements for those clinical inputs are provided.

Code availability does not substitute for data deposition. Public comparison
datasets are available under their existing accessions; clinical study datasets
require separate access arrangements. No clinical-data accession or reviewer
access link is asserted by this repository. No software licence is asserted
pending the authors' determination of the applicable rights.
