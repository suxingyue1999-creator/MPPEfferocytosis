# Input contracts and execution order

All clinical inputs must already be deidentified before they are supplied. Use
`D01` to `D10` for discovery and `V01` to `V23` for validation consistently across
matrix directories, cell prefixes, metadata, AIRR files and intermediate outputs.
No private mapping, real metadata or example clinical record is included.

Discovery metadata requires `sample_id`, `severity`, `severity_en` and `cohort_id`.
`severity_en` uses `mild` and `severe`. The discovery audit reads the coded
`10samples_metadata.csv` under `BD_SCRIPT_ROOT`; the cross-cohort severity lookup
reads the equivalent file specified by `MPP_DISCOVERY_METADATA`.

Validation metadata requires coded `sample_code`/`sample_id`, `severity`,
`allergy_group`, `plastic_group`, `batch_wave`, `age` and `sex`. The original
severity/covariate coding is retained: severity `轻症`/`重症`, allergy and plastic
phenotypes `否`/`是`, and batch levels `2023_Dec`/`2024_Spring`. Supply `batch_wave`
directly; no exact clinical date is needed. Original patient labels must not be
reintroduced into any input. Legacy annotation flags use the anonymous
`legacy_annotation_usable` column.

The validation audit retains the original legacy-object/barcode checks. An
authorized rerun therefore also requires a deidentified legacy Seurat object,
current RNA/ADT matrices, their annotation tables and per-sample Scrublet calls.
The same identities must be recoded in those inputs. The main RNA analyses use
only Gene Expression counts; ADT features are not used as primary RNA counts.

| Configuration | Meaning |
|---|---|
| `MPP_WORKSPACE` | Deidentified working-data tree for validation, imports and figure/table scripts; defaults to `inputs/workspace` |
| `BD_INPUT_ROOT` | Discovery BD filtered/unfiltered MEX archives and audited companion files, organized by coded sample |
| `BD_ANALYSIS_ROOT` | Discovery output tree, kept separate from source input files |
| `BD_SCRIPT_ROOT` | Discovery source/configuration root containing coded `10samples_metadata.csv` |
| `BD_ENV_ROOT` | Existing discovery environment records used by the historical report assembly |
| `BD_DISCOVERY_ROOT` | Completed discovery analysis root required by the cross-cohort scripts |
| `BD_VALIDATION_ROOT` | Validation analysis root with the completed `results` tree |
| `BD_VALIDATION_CLINICAL` | Deidentified validation metadata for cross-cohort audit |
| `BD_CROSS_ROOT` | Cross-cohort analysis output root |
| `BD_R_LIB` | Optional existing R library directory for historical scripts |
| `MPP_DISCOVERY_METADATA` | Coded discovery metadata for severity lookup; defaults to `inputs/discovery.csv` |
| `MPP_PUBLIC_ROOT` | Public-data inputs and outputs, with the original relative layout; defaults to `inputs/public` |
| `MPP_PROTEOMICS_ROOT` | Proteomic report/data tree for the protein import script; defaults to `inputs/proteomics` |
| `SALMON` | Salmon executable for PRJNA792730; defaults to `salmon` on PATH |

Optional worker controls are preserved under `BD_SAMPLE_WORKERS`,
`BD_SCRUBLET_WORKERS`, `BD_SCT_WORKERS`, `BD_MARKER_WORKERS`,
`BD_LINEAGE_WORKERS`, `BD_CELLTYPE_WORKERS` and `BD_HASH_WORKERS`.
The original `BD_REUSE_STAGE02_DISCOVERY` and `LINEAGE_ONLY` stage controls
are also retained; they should be used only with matching audited inputs.

The original relative report/data-tree names are retained in source code. Setting
a root does not change these layouts. The file-reading statements specify the
exact required filenames; absent input should produce an explicit error. Use
fresh analysis/output roots for a rerun and retain read-only source copies.

Run the stages in dependency order:

1. Discovery: audit inputs/object, run Scrublet, QC, apply the original recorded
   QC-review policy, integrate, annotate and compare. Corresponding validation
   scripts audit the generated products. The QC policy retains the historical
   coded-sample sensitivity exclusions; it is not a general automatic policy for
   a newly recruited cohort.
2. Validation: `00_audit_23.R`, `01_rebuild_rna_qc.R`,
   `02_cluster_annotate_rna.R`, `03_pseudobulk_de.R`,
   `03b_subcluster_lineages.R`, `07_lineage_pseudobulk_de.R`,
   `08_lineage_gsea.R`. The source preserves the original phenotype contrasts.
3. V(D)J: `04_analyze_vdj.R` then `09_vdj_clinical_comparison.R`.
   Interaction analysis: `12_rebuild_cellchat_by_sample.R`, then
   `13_summarize_cellchat_results.py`.
4. Discovery report/integrity preparation, then cross-cohort audit, discovery
   freeze, harmonized comparison, pathway integration and replicated-state
   analysis. Historical freeze scripts expect archived validation/manifests;
   those artifacts must be supplied or generated through the original pipeline.
5. `Sensitivity and figures/01_prepare_derived_data.R` performs the actual
   leave-one-sample-out refits for the seven pathways and three reported models;
   `02_build_main_figures.R` plots their summaries. These scripts also require
   completed atlas and pseudobulk products.
6. Result-import and table/model-figure scripts consume existing result files.
   They do not replace the missing original score/trajectory/protein-processing
   implementations described in README.
7. Public datasets: GEO uses the published counts and series metadata;
   PRJNA792730 uses the nine predefined runs, ENA byte/MD5 manifest, GENCODE v44
   transcript reference and FASTQ audit before Salmon. Build gene counts, then
   call `analyze_prjna792730_counts.R <count_matrix> <output_dir>`, followed by
   artifact generation. `build_prjna792730_gene_counts.R` describes its required
   command-line arguments in its usage check.

This is a recovered source-code release, not a bundled data distribution or a
claimed one-command reproduction. Scripts that import older report outputs must
receive those deidentified outputs separately. Clinical-data access and any
original raw-data deposition remain separate from code availability.
